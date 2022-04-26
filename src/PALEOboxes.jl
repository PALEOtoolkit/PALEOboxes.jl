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

include("reactionmethods/VariableStatsMethods.jl")
include("reactionmethods/RateStoich.jl")

include("utils/Interpolation.jl")
include("utils/TestUtils.jl")
include("utils/SIMDutils.jl")
include("utils/IteratorUtils.jl")

include("solverview/VariableAggregator.jl")
include("solverview/SolverView.jl")

include("reactioncatalog/Reactions.jl")

end # module
