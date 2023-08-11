import PALEOboxes as PB
using Test

# import Infiltrator
include("ReactionRateStoichMock.jl")


@testset "RateStoich" begin
    model = PB.create_model_from_config("configratestoich.yaml", "model1")

    domain = PB.get_domain(model, "ocean")   
    domain.grid = PB.Grids.UnstructuredVectorGrid(ncells=10)  
    @test PB.get_length(domain) == 10
    modeldata = PB.create_modeldata(model)
    PB.allocate_variables!(model, modeldata, 1)
    PB.check_ready(model, modeldata)

    all_vars = PB.VariableAggregatorNamed(modeldata)
    all_data = all_vars.values

    PB.initialize_reactiondata!(model, modeldata; create_dispatchlists_all=true)
    
    PB.dispatch_setup(model, :setup, modeldata)
    PB.dispatch_setup(model, :norm_value, modeldata)
    PB.dispatch_setup(model, :initial_value, modeldata)

    dispatchlists = modeldata.dispatchlists_all
    PB.do_deriv(dispatchlists)

    ocean_H2S = all_data.ocean.H2S
    # @Infiltrator.infiltrate
    @test PB.get_total.(ocean_H2S) == fill(10.0, 10)
    @test PB.get_delta.(ocean_H2S) == fill(-30.0, 10)

    ocean_H2S_delta = all_data.ocean.H2S_delta
    @test maximum(ocean_H2S_delta .- fill(-30.0, 10)) < 1e-13

    ocean_myrate = all_data.ocean.myrate
    @test ocean_myrate == fill(42.0, 10)

    ocean_H2S_sms = all_data.ocean.H2S_sms
    @test PB.get_total.(ocean_H2S_sms) == fill(-1.0*42.0, 10)   
    deltaabserr = abs.(PB.get_delta.(ocean_H2S_sms) .- fill(-30.0 - 10.0, 10))
    @test maximum(deltaabserr) < 1e-13

    ocean_SO4_sms = all_data.ocean.SO4_sms
    @test PB.get_total.(ocean_SO4_sms) == fill(1.0*42.0, 10)
    deltaabserr = abs.(PB.get_delta.(ocean_SO4_sms) .- fill(-30.0 - 10.0, 10))
    @test maximum(deltaabserr) < 1e-13

    @info "stoichiometry"
    ratevars = PB.get_variables(domain, v->PB.has_attribute(v, :rate_processname))
    @test length(ratevars) == 1
    
    rv = ratevars[1]
    @test rv.name == "myrate"
    @test PB.get_attribute(rv, :rate_processname) == "redox"
    speciesstoich = Dict(zip(PB.get_attribute(rv, :rate_species), PB.get_attribute(rv, :rate_stoichiometry)))
    @test speciesstoich == Dict("SO4::Isotope" => 1.0, "H2S::Isotope" => -1.0, "O2" => -2.0)      
end
