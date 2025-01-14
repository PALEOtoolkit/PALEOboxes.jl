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

# https://discourse.julialang.org/t/is-compat-jl-worth-it-for-the-public-keyword/119041/22
# define @public for backwards compatibility with Julia < v1.11
macro public(ex)
    if VERSION >= v"1.11"
        args = ex isa Symbol ? (ex,) : Base.isexpr(ex, :tuple) ? ex.args : error("@public parse error")
        esc(Expr(:public, args...))
    else
        nothing
    end
end

@public DocStrings

@public AbstractModel, Model, create_model_from_config, check_variable_links, check_ready, initialize_reactiondata!
@public DispatchMethodLists, dispatch_setup, create_dispatch_methodlists, do_deriv, dispatch_methodlist
@public get_num_domains, show_methods_setup, show_methods_initialize, show_methods_do, show_variables, show_links, show_parameters
@public get_data, get_table

@public AbstractDomain, get_domain, AbstractSubdomain, set_data_dimension!, get_data_dimension

@public VariableDomain, VariableDomPropDep, VariableDomContribTarget, set_data!, allocate_variables!
@public VariableAggregator, VariableAggregatorNamed, ParameterAggregator

@public NamedDimension, get_dimensions, get_dimension, set_coordinates!, get_coordinates

@public AbstractMesh, AbstractMeshOrNothing, AbstractSubdomain, get_mesh, has_internal_cartesian, internal_size, cartesian_size, Grids

@public AbstractCellRange, create_default_cellrange, CellRange, CellRangeColumns

@public Reservoirs, Fluxes, FluxPerturb, Forcings, GridForcings, GridReactions, VariableStats

@public AbstractReaction, ReactionBase, get_reaction, set_model_geometry, register_methods!, register_dynamic_methods!,
    add_method_setup!, add_method_initialize!, add_method_do!, create_reaction
@public show_all_reactions, doc_reaction

@public AbstractParameter, Parameter, ParametersTuple, 
    ParBool, ParDouble, ParInt, ParString, ParEnum, ParDoubleVec, ParDoubleVecVec, ParIntVec, ParStringVec,
    setvalue!, setvalueanddefault!, setfrozen!, set_parameter_value!, get_parameter_value

@public VariableBase, show_variable, get_variable, get_variables, has_variable, get_reaction_variables
@public VariableType, VariableFunction, VariablePhase
@public Attribute, StandardAttributes, get_attribute, set_attribute!, has_attribute, set_variable_attribute!, get_variable_attribute
@public VariableReaction, 
    VarProp, VarPropScalar, VarPropStateIndep, VarPropScalarStateIndep,
    VarDep, VarDepColumn, VarDepScalar, VarDepStateIndep, VarDepColumnStateIndep, VarDepScalarStateIndep,
    VarTarget, VarTargetScalar,
    VarContrib, VarContribColumn, VarContribScalar,
    VarStateExplicit, VarStateExplicitScalar,
    VarTotal, VarTotalScalar,
    VarStateTotal, VarStateTotalScalar,
    VarDeriv, VarDerivScalar,
    VarState, VarStateScalar,
    VarConstraint, VarConstraintScalar
@public AbstractVarList, Varlist_single, VarList_namedtuple, VarList_namedtuple_fields, VarList_components, VarList_tuple, VarList_vector,
    VarList_ttuple, VarList_vvector, VarList_nothing, VarList_tuple_nothing, VarList_fields

@public AbstractReactionMethod, ReactionMethod
@public add_method_setup_initialvalue_vars_default!, add_method_initialize_zero_vars_default!, add_method_do_totals_default!
@public RateStoich, create_ratestoich_method, parse_number_name
@public LinInterp, interp, value_ad, zero_ad, smoothstepcubic

@public AbstractField, Field, get_field, add_field!

@public AbstractSpace, ScalarSpace, CellSpace, ColumnSpace

@public AbstractData, allocate_values, check_values, check_data_type, init_values!, zero_values!, dof_values!, 
    copyfieldto!, copytofield!, add_field_vec!, num_components, get_components
@public UndefinedData, ScalarData, AbstractIsotopeScalar, IsotopeLinear, ArrayScalarData
@public isotope_totaldelta, get_total, get_delta

@public AbstractModelData, ModelData, create_modeldata, add_arrays_data!

@public Constants, TestUtils, IteratorUtils, ChemistryUtils, DocStrings, SIMDutils

@public collate_markdown, precompile_reaction, run_model



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
include("variableaggregators/ParameterAggregator.jl")

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

include("deprecated.jl")


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
        "ReactionReservoirForced", "ReactionSum", "ReactionFluxTarget", "ReactionForceInterp",
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
