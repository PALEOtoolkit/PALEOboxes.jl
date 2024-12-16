
using Test
using Logging

import PALEOboxes as PB


@testset "Fields" begin

    @testset "ScalarData, ScalarSpace" begin

        f = PB.allocate_field(PB.ScalarData, (), Float64, PB.ScalarSpace, nothing; allocatenans=true)
        
        @test isnan(f.values[])

    end

    @testset "ScalarData, CellSpace, no grid" begin
        # check that a CellSpace Field with no grid behaves as a ScalarSpace Field

        f = PB.allocate_field(PB.ScalarData, (), Float64, PB.CellSpace, nothing; allocatenans=true)
        
        @test isnan(f.values[])
    end

    @testset "ScalarData, CellSpace, UnstructuredColumnGrid" begin
        g = PB.Grids.UnstructuredColumnGrid(ncells=5, Icolumns=[collect(1:5)])

        f = PB.allocate_field(PB.ScalarData, (), Float64, PB.CellSpace, g; allocatenans=true)
        
        @test isnan.(f.values) == [true for i in 1:5]
      
    end
  
    @testset "ArrayScalarData, ScalarSpace" begin
        f = PB.allocate_field(PB.ArrayScalarData, (PB.NamedDimension("test", 2, []), ), Float64, PB.ScalarSpace, nothing; allocatenans=true)

        @test isnan.(f.values) == [true, true]
        
    end

    @testset "ArrayScalarData, CellSpace, no grid" begin
        # TODO this is possibly inconsistent with (ScalarData, CellSpace, no grid),
        # as Field.values here is a (1, 2) Array, not a (2,) Vector
        f = PB.allocate_field(PB.ArrayScalarData, (PB.NamedDimension("test", 2, []), ), Float64, PB.CellSpace, nothing; allocatenans=true)

        @test_broken size(f.values) == (2, ) # TODO should be a Vector ?
        @test size(f.values) == (1, 2)
        @test isnan.(f.values) == [true true] # TODO 1x2 Array or Vector ?        
        
    end

    @testset "ArrayScalarData, CellSpace, UnstructuredColumnGrid" begin
        g = PB.Grids.UnstructuredColumnGrid(ncells=5, Icolumns=[collect(1:5)])

        f = PB.allocate_field(PB.ArrayScalarData, (PB.NamedDimension("test", 2, []), ), Float64, PB.CellSpace, g; allocatenans=true)

        @test size(f.values) == (5, 2)

       
    end

end
