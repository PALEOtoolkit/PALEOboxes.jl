
using Test

import PALEOboxes as PB

include("ReactionPaleoMockModule.jl")

@testset "Grids" begin

@testset "Cartesian Grid" begin

    model = PB.create_model_from_config(joinpath(@__DIR__, "configgrids.yaml"), "cartesiangrid")

    # Test domains
    @test PB.get_num_domains(model) == 1

    ocean_domain = PB.get_domain(model, "ocean")
    dsize = (2, 3, 4)
    @test PB.get_length(ocean_domain) == prod(dsize)

    hostvars = PB.get_variables(ocean_domain, hostdep=true)
    modelvars = Dict((var.name, var) for var in PB.get_variables(ocean_domain, hostdep=false))

    modeldata = PB.create_modeldata(model)

    PB.allocate_variables!(ocean_domain, modeldata, hostdep=false)
    @test length(PB.get_unallocated_variables(ocean_domain, modeldata)) == 1

    @test length(hostvars) == 1
    @test hostvars[1].name == "scalar_dep"
    @test PB.is_scalar(hostvars[1])
    hostvardata = Float64[NaN]
    PB.set_data!(hostvars[1], modeldata, hostvardata)
    @test PB.check_ready(ocean_domain, modeldata, throw_on_error=false) == true

    @test PB.check_configuration(model, throw_on_error=false) == true

    PB.initialize_reactiondata!(model, modeldata)
    
    cellrange = PB.Grids.create_default_cellrange(ocean_domain, ocean_domain.grid)
    cellranges = [cellrange]

    PB.dispatch_setup(model, :setup, modeldata, cellranges)
    PB.dispatch_setup(model, :norm_value, modeldata, cellranges)
    PB.dispatch_setup(model, :initial_value, modeldata, cellranges)

    dispatchlists = modeldata.dispatchlists_all

    hostvardata[] = 23.0
    PB.do_deriv(dispatchlists)
    @test PB.get_data(modelvars["julia_paleo_mock/scalar_prop"], modeldata) == hostvardata

    mock_phy_var = modelvars["mock_phy"]
    @test size(PB.get_data(mock_phy_var, modeldata)) == dsize
    expected_result = Array{Float64, length(dsize)}(undef, dsize)
    fill!(expected_result, 0.15)
    @test PB.get_data(mock_phy_var, modeldata) == expected_result 

    model = nothing
end

    
@testset "Cartesian Linear 2D" begin

    model = PB.create_model_from_config(joinpath(@__DIR__, "configgrids.yaml"), "cartesianlinear2D")

    # Test domains
    @test PB.get_num_domains(model) == 1

    surface_domain = PB.get_domain(model, "surface")
    dsize = (144, 90)
    @test PB.get_length(surface_domain) == prod(dsize)

    modeldata = PB.create_modeldata(model)
    PB.allocate_variables!(model, modeldata)
    @test PB.check_ready(model, modeldata, throw_on_error=false) == true

    modelcreated_vars_dict = Dict([(var.name, var) for var in PB.get_variables(surface_domain, hostdep=false)])

    PB.initialize_reactiondata!(model, modeldata)    
      
    PB.dispatch_setup(model, :setup, modeldata)
    PB.dispatch_setup(model, :norm_value, modeldata)
    PB.dispatch_setup(model, :initial_value, modeldata)

    Asurf = PB.get_data(modelcreated_vars_dict["Asurf"], modeldata)
    @test length(Asurf) == prod(dsize)
    @test size(Asurf) == (prod(dsize), ) # linear Vector
    @test isapprox(sum(Asurf), 4*pi*6.371229e6^2, rtol=1e-14)

    model = nothing
end

@testset "Cartesian Array 2D" begin

    model = PB.create_model_from_config(joinpath(@__DIR__, "configgrids.yaml"), "cartesianarray2D")

    # Test domains
    @test PB.get_num_domains(model) == 1

    surface_domain = PB.get_domain(model, "surface")
    dsize = (144, 90)
    @test PB.get_length(surface_domain) == prod(dsize)

    modeldata = PB.create_modeldata(model)
    PB.allocate_variables!(model, modeldata)
    @test PB.check_ready(model, modeldata, throw_on_error=false) == true

    modelcreated_vars_dict = Dict([(var.name, var) for var in PB.get_variables(surface_domain, hostdep=false)])

    PB.initialize_reactiondata!(model, modeldata)    
      
    PB.dispatch_setup(model, :setup, modeldata)
    PB.dispatch_setup(model, :norm_value, modeldata)
    PB.dispatch_setup(model, :initial_value, modeldata)

    Asurf = PB.get_data(modelcreated_vars_dict["Asurf"], modeldata)
    @test length(Asurf) == prod(dsize)
    @test size(Asurf) == dsize # 2D Array
    @test isapprox(sum(Asurf), 4*pi*6.371229e6^2, rtol=1e-14)

    model = nothing
end

end
