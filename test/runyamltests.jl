
using Test
import YAML

@testset "YAML.jl" begin
    # test YAML.jl handling of merge keys and duplicate keys
    # see comments in merge.yaml test file

    data = YAML.load_file("merge.yaml")

    @test length(data) == 11

    # CENTER
    @test data[1] == Dict("x" => 1, "y" => 2)
    # LEFT
    @test data[2] == Dict("x" => 0, "y" => 2)
    # BIG
    @test data[3] == Dict("r" => 10)
    # SMALL
    @test data[4] == Dict("r" => 1)

    # All the following maps are equal:

    # explicit keys
    @test data[5] == Dict("label" => "center/big", "x" => 1, "r" => 10, "y" => 2)
    # Merge one map
    @test data[6] == Dict("label" => "center/big", "x" => 1, "r" => 10, "y" => 2)
    # Merge multiple maps
    @test data[7] == Dict("label" => "center/big", "x" => 1, "r" => 10, "y" => 2)
    # Override
    @test data[8] == Dict("label" => "center/big", "x" => 1, "r" => 10, "y" => 2)

    # SD added: testing key yaml.jl key override and merge

    # yaml.jl merge keys don't override explicit keys (as spec)
    @test data[9] == Dict("r" => 10)
    # yaml.jl issue ? duplicates are allowed, and later keys override earlier ones
    @test data[10] == Dict("r" => 10)
    # yaml.jl issue ? allows duplicate merge keys, and later ones override earlier ones !
    @test data[11] == Dict("r" => 1)
end
