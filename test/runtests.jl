
using Test
# using NBInclude


import PALEOboxes as PB

@testset "PALEOboxes" begin

include("runfieldtests.jl")
include("rundoctests.jl")
include("runyamltests.jl")
include("runbasetests.jl")
include("rungridstests.jl")
include("runreservoirtests.jl")
include("runfluxtests.jl")
include("runratestoichtests.jl")
include("runsimdtests.jl")

end # @testset
