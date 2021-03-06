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

include("utils/DocStrings.jl")

include("Types.jl")
include("CoordsDims.jl")
include("Fields.jl")

include("data/AtomicScalar.jl")
include("data/ScalarData.jl")
include("data/ArrayScalarData.jl")
include("data/IsotopeData.jl")

include("VariableAttributes.jl")
include("VariableReaction.jl")
include("VariableDomain.jl")
include("ReactionMethodSorting.jl")
include("Model.jl")
include("Domain.jl")
include("CellRange.jl")
include("Parameter.jl")
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

# create Reactions and register methods to precompile this code
function precompile_reaction(rdict, classname)
    rj = create_reaction(rdict, classname, "test", Dict{String, Any}())
    rj.base.domain = Domain(name="test", ID=1, parameters=Dict{String, Any}())
    register_methods!(rj)

    return nothing
end

function precompile_reactions()
    rdict = find_all_reactions()
    reactionlist = [
        "ReactionFluxTransfer", "ReactionReservoirScalar", "ReactionFluxPerturb", "ReactionReservoir",
        "ReactionReservoirForced", "ReactionSum", "ReactionFluxTarget", "ReactionForceInterp", "ReactionGrid2DNetCDF",
        "ReactionAreaVolumeValInRange", "ReactionReservoirWellMixed", "ReactionForceGrid", "ReactionConst", 
        "ReactionRestore", "ReactionScalarConst", "ReactionForceInsolation", "ReactionVectorSum", "ReactionWeightedMean",
        "ReactionReservoirTotal", "ReactionUnstructuredVectorGrid", "ReactionCartesianGrid", "ReactionReservoirConst",
    ]
    for r in reactionlist
        precompile_reaction(rdict, r)
    end

    return nothing
end

precompile_reactions()

# create and take a timestep for a test configuration
# function run_model(configfile, configname)
    
#     model =  create_model_from_config(configfile, configname)

#     modeldata =  create_modeldata(model)
#     allocate_variables!(model, modeldata)

#     check_ready(model, modeldata)

#     initialize_reactiondata!(model, modeldata)

#     check_configuration(model)

#     dispatch_setup(model, :setup, modeldata)
#     dispatch_setup(model, :norm_value, modeldata)   
#     dispatch_setup(model, :initial_value, modeldata)

#     # take a time step
#     # dispatchlists = modeldata.dispatchlists_all
#     # do_deriv(dispatchlists)

#     return nothing
# end

# Negligible difference
# run_model(joinpath(@__DIR__, "../test/configreservoirs.yaml"), "model1")
# run_model(joinpath(@__DIR__, "../test/configfluxes.yaml"), "model1")

end # module
