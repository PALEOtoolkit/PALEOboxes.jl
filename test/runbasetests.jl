import PALEOboxes as PB
using Test
using BenchmarkTools

include("ReactionPaleoMockModule.jl")

@testset "PB.jl base" begin
    model = PB.create_model_from_config("configbase.yaml", "model1")

    PB.show_methods_setup(model)
    PB.show_methods_initialize(model)
    PB.show_methods_do(model)

    @test length(model.domains) == 3

    ocean_domain = model.domains[1]
    @test ocean_domain.name == "ocean"
    @test ocean_domain.ID == 1

    @test length(ocean_domain.reactions) == 2
    reaction = PB.get_reaction(ocean_domain, "julia_paleo_mock")
    @test reaction.name == "julia_paleo_mock"
    println(reaction)

    @test length(PB.get_parameters(reaction)) == 4
    par = PB.get_parameter(reaction, "par1")
    @test par.name == "par1"
    @test par.v == 0.15

    
    @test length(PB.get_variables(reaction)) == 3
    var = PB.get_variable(reaction, "do_react", "phy")
    @test var.localname == "phy"
    @test var.linkreq_name == "mock_phy"

    @test length(PB.get_variables(ocean_domain)) == 5

    modelvars = Dict((var.name, var) for var in PB.get_variables(ocean_domain, hostdep=false))
    @test length(modelvars) == 4
    @test !isnothing(PB.get_variable(ocean_domain, "mock_phy"))
    @test haskey(modelvars, "mock_phy")
    @test haskey(modelvars, "julia_paleo_mock/scalar_prop")
    @test haskey(modelvars, "mock2_phy")
    @test haskey(modelvars, "julia_paleo_mock2/scalar_prop")

    hostvars = PB.get_variables(ocean_domain, hostdep=true)
    @test length(hostvars) == 1
    @test hostvars[1].name == "scalar_dep"
    @test PB.is_scalar(hostvars[1])

    domain_size = 1000
    # domain_size = 10
    ocean_domain.grid = PB.Grids.UnstructuredVectorGrid(ncells=domain_size)
    @test PB.get_length(ocean_domain) == domain_size


    modeldata = PB.create_modeldata(model)

    @test length(PB.get_unallocated_variables(ocean_domain, modeldata)) == 5
    PB.allocate_variables!(ocean_domain, modeldata, hostdep=false)
    @test length(PB.get_unallocated_variables(ocean_domain, modeldata)) == 1
    @test PB.check_ready(ocean_domain, modeldata, throw_on_error=false) == false

    hostvardata = Float64[NaN]
    PB.set_data!(hostvars[1], modeldata, hostvardata)
    @test PB.check_ready(ocean_domain, modeldata, throw_on_error=false) == true

    @test PB.check_configuration(model, throw_on_error=false) == true

    PB.initialize_reactiondata!(model, modeldata)
    
    cellrange = PB.Grids.create_default_cellrange(ocean_domain, ocean_domain.grid)
    cellranges = [cellrange]

    PB.dispatch_setup(model, :initial_value, modeldata, cellranges)

    dispatchlists = modeldata.dispatchlists_all

    hostvardata[] = 27.3
    PB.do_deriv(dispatchlists)   
    @test PB.get_data(modelvars["julia_paleo_mock/scalar_prop"], modeldata) == hostvardata
    @test PB.get_data(modelvars["julia_paleo_mock2/scalar_prop"], modeldata) == hostvardata


    mock_phy_var = modelvars["mock_phy"]
    expected_result = Vector{Float64}(undef, domain_size)
    expected_result .= 0.15
    @test PB.get_data(mock_phy_var, modeldata) == expected_result 

end


@testset "PB.jl loop" begin
    
    @test_throws ErrorException("The input graph contains at least one loop.") PB.create_model_from_config("configbase.yaml", "model_with_loop")
   
end
