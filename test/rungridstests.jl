
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

    modeldata = PB.create_modeldata(model)

    PB.allocate_variables!(ocean_domain, modeldata, 1, hostdep=false)
    @test length(PB.get_unallocated_variables(ocean_domain, modeldata, 1)) == 1

    @test length(hostvars) == 1
    @test hostvars[1].name == "scalar_dep"
    @test PB.is_scalar(hostvars[1])
    hostvardata = Float64[NaN]
    PB.set_data!(hostvars[1], modeldata, 1, hostvardata)
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

    all_vars = PB.VariableAggregatorNamed(modeldata)
    @info "allvars: $all_vars"
    all_data = all_vars.values

    # NB: "ocean.julia_paleo_mock/scalar_prop" substitute / --> __
    @test all_data.ocean.julia_paleo_mock__scalar_prop == hostvardata

    @test all_data.ocean.mock_phy == fill(0.15, dsize)

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
    PB.allocate_variables!(model, modeldata, 1)
    @test PB.check_ready(model, modeldata; throw_on_error=false) == true

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
    PB.allocate_variables!(model, modeldata, 1)
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
