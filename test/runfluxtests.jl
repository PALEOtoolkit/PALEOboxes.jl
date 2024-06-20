
using Test
using Logging

import PALEOboxes as PB


@testset "Fluxes" begin

    model =  PB.create_model_from_config("configfluxes.yaml", "model1")

    global_domain =  PB.get_domain(model, "global")
    @test  PB.get_length(global_domain) == 1
    ocean_domain =  PB.get_domain(model, "ocean")
    ocean_domain.grid = PB.Grids.UnstructuredVectorGrid(ncells=1000) 
    ocean_length = PB.get_length(ocean_domain)
    @test  ocean_length == 1000

    modeldata =  PB.create_modeldata(model)

    PB.allocate_variables!(model, modeldata, 1; check_units_opt=:error)
 
    @test  PB.check_ready(model, modeldata) == true    

    # get variables
    all_vars = PB.VariableAggregatorNamed(modeldata)
    @test PB.num_vars(all_vars) == 9
    @info "all_vars: $all_vars"
    all_data = all_vars.values

    PB.initialize_reactiondata!(model, modeldata; create_dispatchlists_all=true)
      
    @info "dispatch_setup"
    PB.dispatch_setup(model, :setup, modeldata)
    PB.dispatch_setup(model, :norm_value, modeldata)
    PB.dispatch_setup(model, :initial_value, modeldata)
    
    @info "const fluxes:"
    @test all_data.ocean.flux_A == fill(10.0, ocean_length)
    @test all_data.ocean.flux_B == fill(20.0, ocean_length)
    @test all_data.ocean.flux_C == 
        fill(PB.isotope_totaldelta(PB.IsotopeLinear, 2.0, -1.0), ocean_length)
   
    dispatchlists = modeldata.dispatchlists_all
    PB.do_deriv(dispatchlists)
    PB.do_deriv(dispatchlists)

    @info "transferred fluxes"
    @test all_data.ocean.tflux_A == fill(2*10.0, ocean_length)
    @test all_data.ocean.tflux_C == 
        fill(2*PB.isotope_totaldelta(PB.IsotopeLinear, 2.0, -1.0), ocean_length)

    @info "transferred flux totals"
    @test all_data.ocean.tflux_total_A[] == ocean_length*2*10.0
    @test all_data.ocean.tflux_total_C[] == 
        ocean_length*2*PB.isotope_totaldelta(PB.IsotopeLinear, 2.0, -1.0)

    @info "global transferred fluxes"
    @test all_data.global.tflux_B[] == ocean_length*-1*20.0
    @test all_data.global.tflux_C[] == 
        ocean_length*-1*PB.isotope_totaldelta(PB.IsotopeLinear, 2.0, -1.0)

    @info "test complete"
  
end
