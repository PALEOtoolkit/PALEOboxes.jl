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
