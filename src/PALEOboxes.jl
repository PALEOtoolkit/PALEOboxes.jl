"""
    PALEOboxes

PALEOboxes provides the model coupler for the PALEO model framework.

It is registered as a Julia package with public github repository
[PALEOboxes.jl](https://github.com/PALEOtoolkit/PALEOboxes.jl) 
and online 
[documentation](https://paleotoolkit.github.io/PALEOboxes.jl)

A PALEO `Model` contains `Domain`s, each of which contain Variables 
defining `Field`s containing `Data` arrays, and Reactions with `ReactionMethod`s
that operate on the Variables to calculate model time evolution.

PALEOboxes creates the model from a .yaml configuration file, and implements a coupler that provides 
a unified mechanism for:
1. ‘low-level’ coupling (e.g. linking individual redox Reactions within a Domain, on which is built
2. ‘module-level’ coupling (linking e.g. atmosphere and ocean components) based on standardising names for communicating fluxes, and which enables
3. separation of biogeochemical reaction and transport. 
"""
module PALEOboxes

import YAML
import Graphs # formerly LightGraphs
import DataFrames
using DocStringExtensions
import OrderedCollections
import Logging
import Printf
import Atomix

import PrecompileTools
import TimerOutputs: @timeit, @timeit_debug

include("utils/DocStrings.jl")

include("Types.jl")
include("CoordsDims.jl")
include("Fields.jl")

include("data/ScalarData.jl")
include("data/ArrayScalarData.jl")
include("data/IsotopeData.jl")

include("VariableAttributes.jl")
include("VariableReaction.jl")
include("VariableDomain.jl")

include("Parameter.jl")
include("ParameterAggregator.jl")

include("ReactionMethodSorting.jl")
include("Model.jl")
include("Domain.jl")
include("CellRange.jl")
include("ReactionMethod.jl")
include("Reaction.jl")
include("ReactionFactory.jl")
include("ModelData.jl")
include("Grids.jl")

include("reactionmethods/SetupInitializeUtilityMethods.jl")
include("reactionmethods/VariableStatsMethods.jl")
include("reactionmethods/RateStoich.jl")

include("utils/Interpolation.jl")
include("utils/TestUtils.jl")
include("utils/SIMDutils.jl")
include("utils/IteratorUtils.jl")
include("utils/DocUtils.jl")
include("utils/StringUtils.jl")
include("utils/ChemistryUtils.jl")
include("utils/ADUtils.jl")

include("variableaggregators/VariableAggregator.jl")

include("reactioncatalog/Reactions.jl")

# Deprecated functions
"""
    get_statevar

DEPRECATED - moved to PALEOmodel
"""
function get_statevar end

"""
    get_statevar_norm

DEPRECATED - moved to PALEOmodel
"""
function get_statevar_norm end

#####################################################
# Precompilation
# Run code to precompile
#######################################################

"""
    precompile_reaction(rdict::Dict{String, Type}, classname::AbstractString; logger=Logging.NullLogger())
    precompile_reaction(ReactionType::Type{<:AbstractReaction}; logger=Logging.NullLogger())

For use in @PrecompileTools.compile_workload: create Reaction and call register_methods
"""
function precompile_reaction(rdict::Dict{String, Type}, classname::AbstractString; logger=Logging.NullLogger())

    try
        Logging.with_logger(logger) do
            rj = create_reaction(rdict, classname, "test", Dict{String, Any}())
            rj.base.domain = Domain(name="test", ID=1, parameters=Dict{String, Any}())
            register_methods!(rj)
        end
    catch ex
        @warn "precompile_reaction(rdict, $classname) failed with exception:" ex
    end

    return nothing
end

function precompile_reaction(ReactionType::Type{<:AbstractReaction}; logger=Logging.NullLogger())

    try
        Logging.with_logger(logger) do
            rj = create_reaction(ReactionType, "test", Dict{String, Any}())
            rj.base.domain = Domain(name="test", ID=1, parameters=Dict{String, Any}())
            register_methods!(rj)
        end
    catch ex
        @warn "precompile_reaction($ReactionType) failed with exception:" ex
    end

    return nothing
end

"""
    run_model(configfile::AbstractString, configname::AbstractString; call_do_deriv=true, tforce=0.0, logger=Logging.NullLogger())
    run_model(model::Model; call_do_deriv=true, tforce=0.0, logger=Logging.NullLogger())

For use in @PrecompileTools.compile_workload: create_model_from_config and optionally call do_deriv at time tforce
"""
function run_model(configfile::AbstractString, configname::AbstractString; logger=Logging.NullLogger(), kwargs...)
    
    try
        model = Logging.with_logger(logger) do 
            create_model_from_config(configfile, configname)
        end
        run_model(model; logger, kwargs...)
    catch ex
        @warn "run_model($configfile, $configname) failed with exception:" ex
    end
    
    return nothing
end

function run_model(model::Model; call_do_deriv=true, tforce=0.0, logger=Logging.NullLogger())

    try
        Logging.with_logger(logger) do
            arrays_idx = 1
            modeldata =  create_modeldata(model)
            allocate_variables!(model, modeldata, arrays_idx)

            check_ready(model, modeldata)

            initialize_reactiondata!(model, modeldata; create_dispatchlists_all=true)

            check_configuration(model)

            dispatch_setup(model, :setup, modeldata)
            dispatch_setup(model, :norm_value, modeldata)   
            dispatch_setup(model, :initial_value, modeldata)

            # take a time step at forcing time tforce - TODO, can be model dependent on missing setup
            if call_do_deriv
                # set tforce (if present)
                hostdep_vars = VariableDomain[]
                for dom in model.domains
                    dv, _ = get_host_variables(dom, VF_Undefined)
                    append!(hostdep_vars, dv )
                end
                hostdep = VariableAggregatorNamed(hostdep_vars, modeldata, arrays_idx)
                set_values!(hostdep, Val(:global), Val(:tforce), tforce; allow_missing=true)

                dispatchlists = modeldata.dispatchlists_all
                do_deriv(dispatchlists)
            end
        end
    catch ex
        @warn "run_model($model; call_do_deriv=$call_do_deriv) failed with exception:" ex
    end

    return nothing
end


@PrecompileTools.setup_workload begin
    # create Reactions and register methods to precompile this code

    # Putting some things in `setup` can reduce the size of the
    # precompile file and potentially make loading faster.
    

    rdict = find_all_reactions()
    reactionlist = [
        "ReactionFluxTransfer", "ReactionReservoirScalar", "ReactionFluxPerturb", "ReactionReservoir",
        "ReactionReservoirForced", "ReactionSum", "ReactionFluxTarget", "ReactionForceInterp", "ReactionGrid2DNetCDF",
        "ReactionAreaVolumeValInRange", "ReactionReservoirWellMixed", "ReactionForceGrid", "ReactionConst", 
        "ReactionRestore", "ReactionScalarConst", "ReactionVectorSum", "ReactionWeightedMean",
        "ReactionReservoirTotal", "ReactionUnstructuredVectorGrid", "ReactionCartesianGrid", "ReactionReservoirConst",
    ]

    @PrecompileTools.compile_workload begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        for r in reactionlist
            precompile_reaction(rdict, r)
        end

        # Negligible difference ?
        run_model(joinpath(@__DIR__, "../test/configreservoirs.yaml"), "model1")
        run_model(joinpath(@__DIR__, "../test/configfluxes.yaml"), "model1")
    end
end

end # module
