import PALEOboxes as PB
using Test
using BenchmarkTools

include("ReactionPaleoMockModule.jl")

@testset "PALEOboxes base" begin
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

    @test length(PB.get_unallocated_variables(ocean_domain, modeldata, 1)) == 5
    PB.allocate_variables!(ocean_domain, modeldata, 1; hostdep=false, check_units_opt=:error)
    @test length(PB.get_unallocated_variables(ocean_domain, modeldata, 1)) == 1
    @test PB.check_ready(ocean_domain, modeldata, throw_on_error=false) == false

    hostvardata = Float64[NaN]
    PB.set_data!(hostvars[1], modeldata, 1, hostvardata)
    @test PB.check_ready(ocean_domain, modeldata, throw_on_error=false) == true

    @test PB.check_configuration(model, throw_on_error=false) == true

    PB.initialize_reactiondata!(model, modeldata; create_dispatchlists_all=true)
    
    cellrange = PB.Grids.create_default_cellrange(ocean_domain, ocean_domain.grid)
    cellranges = [cellrange]

    PB.dispatch_setup(model, :setup, modeldata, cellranges)
    PB.dispatch_setup(model, :norm_value, modeldata, cellranges)
    PB.dispatch_setup(model, :initial_value, modeldata, cellranges)

    dispatchlists = modeldata.dispatchlists_all

    hostvardata[] = 27.3
    PB.do_deriv(dispatchlists)   

    allvars = PB.VariableAggregatorNamed(modeldata)
    println(allvars)
    @test PB.num_vars(allvars) == 5
    @test length(allvars) == 3 + 2*domain_size

    @test PB.get_data(modelvars["julia_paleo_mock/scalar_prop"], modeldata) == hostvardata
    @test PB.get_data(modelvars["julia_paleo_mock2/scalar_prop"], modeldata) == hostvardata
    # NB: / in name --> __ in field, so "ocean.julia_paleo_mock/scalar_prop" -> field ocean.julia_paleo_mock__scalar_prop
    @test allvars.values.ocean.julia_paleo_mock__scalar_prop == hostvardata
    @test allvars.values.ocean.julia_paleo_mock2__scalar_prop == hostvardata
 
    @test allvars.values.ocean.mock_phy == fill(0.15, domain_size)

end


@testset "PALEOboxes base loop" begin
    
    @test_throws ErrorException("The input graph contains at least one loop.") PB.create_model_from_config("configbase.yaml", "model_with_loop")
   
end

@testset "PALEOboxes base invalid parameter name" begin
    
    @test_throws ErrorException("reaction ocean.julia_paleo_mock has no Parameter(s):\n    misspelt_par:    something\n") PB.create_model_from_config("configbase.yaml", "model_with_invalid_parameter")
   
end

@testset "PALEOboxes base empty parameters" begin
    
    model = PB.create_model_from_config("configbase.yaml", "model_with_empty_parameters")

    @test length(model.domains) == 3
   
end

@testset "PALEOboxes base empty variable link" begin
    
    @test_throws MethodError PB.create_model_from_config("configbase.yaml", "model_with_empty_variable_link")

end

@testset "PALEOboxes base invalid variable attribute" begin
    
    @test_throws ErrorException PB.create_model_from_config("configbase.yaml", "model_with_invalid_variable_attribute")

end

# test parsing concatenated yaml files
@testset "PALEOboxes base concatenate" begin
    model = PB.create_model_from_config(["configbase_pt1.yaml", "configbase_pt2.yaml"], "model1")

    @test length(model.domains) == 3

end

@testset "PALEOboxes base unlinked variable" begin

    model = PB.create_model_from_config("configbase.yaml", "model_with_unlinked_variable")

    @test PB.check_variable_links(model; throw_on_error=true, expect_hostdep_varnames=["ocean.host_supplied_dep"]) == true
    
    @test_throws ErrorException PB.check_variable_links(model) 

end