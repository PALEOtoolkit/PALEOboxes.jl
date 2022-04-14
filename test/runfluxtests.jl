
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

    PB.allocate_variables!(model, modeldata)
 
    @test  PB.check_ready(model, modeldata) == true    

    # get variables
    global_vars_dict =  Dict(var.name=>var for var in PB.get_variables(global_domain))
    ocean_vars_dict =  Dict(var.name=>var for var in PB.get_variables(ocean_domain))

    PB.initialize_reactiondata!(model, modeldata)
      
    @info "dispatch_setup"
    PB.dispatch_setup(model, :initial_value, modeldata)
    
    @info "const fluxes:"
    @test PB.get_data(ocean_vars_dict["flux_A"], modeldata) == fill(10.0, ocean_length)
    @test PB.get_data(ocean_vars_dict["flux_B"], modeldata) == fill(20.0, ocean_length)
    @test PB.get_data(ocean_vars_dict["flux_C"], modeldata) == 
        fill(PB.isotope_totaldelta(PB.IsotopeLinear, 2.0, -1.0), ocean_length)
   
    dispatchlists = modeldata.dispatchlists_all
    PB.do_deriv(dispatchlists)
    PB.do_deriv(dispatchlists)

    @info "transferred fluxes"
    @test PB.get_data(ocean_vars_dict["tflux_A"], modeldata) == fill(2*10.0, ocean_length)
    @test PB.get_data(ocean_vars_dict["tflux_C"], modeldata) == 
        fill(2*PB.isotope_totaldelta(PB.IsotopeLinear, 2.0, -1.0), ocean_length)

    @info "transferred flux totals"
    @test PB.get_data(ocean_vars_dict["tflux_total_A"], modeldata)[] == ocean_length*2*10.0
    @test PB.get_data(ocean_vars_dict["tflux_total_C"], modeldata)[] == 
        ocean_length*2*PB.isotope_totaldelta(PB.IsotopeLinear, 2.0, -1.0)

    @info "global transferred fluxes"
    @test PB.get_data(global_vars_dict["tflux_B"], modeldata)[] == ocean_length*-1*20.0
    @test PB.get_data(global_vars_dict["tflux_C"], modeldata)[] == 
        ocean_length*-1*PB.isotope_totaldelta(PB.IsotopeLinear, 2.0, -1.0)

    @info "test complete"
  
end
