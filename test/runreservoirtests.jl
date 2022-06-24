
using Test
using Logging

import PALEOboxes as PB

import Infiltrator

@testset "Reactions" begin
    @test PB.find_reaction("ReactionReservoirScalar") == PB.Reservoirs.ReactionReservoirScalar

    PB.show_all_reactions()  # just check runs without error
end

@testset "Reservoirs" begin

    # logfile = open("Reservoirs_log.txt", "w")
    # println("writing log output to ", logfile)
    # global_logger(SimpleLogger(logfile))

    model =  PB.create_model_from_config("configreservoirs.yaml", "model1")

    # Test domains
    @test  PB.get_num_domains(model) == 2
    global_domain =  PB.get_domain(model, "global")
    @test global_domain.name == "global"
    @test  PB.get_length(global_domain) == 1
    ocean_domain =  PB.get_domain(model, "ocean")
    @test  PB.get_length(ocean_domain) == 1000
    ocean_length = PB.get_length(ocean_domain)

    # allocate internal variables
    modeldata =  PB.create_modeldata(model)
    PB.allocate_variables!(model, modeldata, hostdep=false)

    @test length( PB.get_unallocated_variables(global_domain, modeldata)) == 4
    @test  PB.check_ready(global_domain, modeldata, throw_on_error=false) == false    

    # allocate arrays for host dependencies and set data pointers
    PB.allocate_variables!(model, modeldata, hostdep=true)
    @test  PB.check_ready(global_domain, modeldata) == true

    # get modelcreated variables
    global_modelcreated_vars_dict =  Dict(var.name=>var for var in PB.get_variables(global_domain, hostdep=false))
    ocean_modelcreated_vars_dict =  Dict(var.name=>var for var in PB.get_variables(ocean_domain, hostdep=false))

    # check state variables
    stateexplicit_vars, stateexplicit_sms_vars =
        PB.get_host_variables(global_domain, PB.VF_StateExplicit, match_deriv_suffix="_sms")
    @test length(stateexplicit_vars) == 2  # A and O

    stateexplicit_vars, stateexplicit_sms_vars =
        PB.get_host_variables(ocean_domain, PB.VF_StateExplicit, match_deriv_suffix="_sms")
    @test length(stateexplicit_vars) == 1

    # get host-dependent variables
    global_hostdep_vars_vec =  PB.get_variables(global_domain, hostdep=true)
    @test length(global_hostdep_vars_vec) == 4
    global_hostdep_data = Dict(var.name=>PB.get_data(var, modeldata) for var in global_hostdep_vars_vec)

    ocean_hostdep_vars_vec =  PB.get_variables(ocean_domain, hostdep=true)
    @test length(ocean_hostdep_vars_vec) == 2
    ocean_hostdep_data = Dict(var.name=>PB.get_data(var, modeldata) for var in ocean_hostdep_vars_vec)

    PB.initialize_reactiondata!(model, modeldata)
    
    @info "dispatch_setup"
    PB.dispatch_setup(model, :setup, modeldata)
    PB.dispatch_setup(model, :norm_value, modeldata)
    PB.dispatch_setup(model, :initial_value, modeldata)

    @info "global host-dependent variables:\n"  global_hostdep_data
    @test global_hostdep_data["A"].v[]  == 3.193e18   # A
    @test global_hostdep_data["A"].v_moldelta[]  == 2.0*3.193e18    # A_moldelta
    
    @info "global const variables:"
    @test PB.get_data(global_modelcreated_vars_dict["ConstS"], modeldata)[] == 1.0

    @info "ocean host-dependent variables:\n"  ocean_hostdep_data
    @test ocean_hostdep_data["T"][1]           == 1.0*10.0    

    # take a time step
    dispatchlists = modeldata.dispatchlists_all
    PB.do_deriv(dispatchlists)

    @info "global model-created variables:\n" global_modelcreated_vars_dict
    @test  PB.get_data(global_modelcreated_vars_dict["A_norm"], modeldata)[] == 10.0
    @test  PB.get_data(global_modelcreated_vars_dict["A_delta"], modeldata)[]== 2.0

    @info "ocean model-created variables:\n" ocean_modelcreated_vars_dict
    @test  PB.get_data(ocean_modelcreated_vars_dict["T_conc"], modeldata) == fill(1.0, ocean_length)
    PB.do_deriv(dispatchlists) # repeat do_deriv to check total is reinitialized each time step
    @test  PB.get_data(ocean_modelcreated_vars_dict["T_total"], modeldata)[] == 1.0*10.0*ocean_length

    @info "sum variables:"
    @test  PB.get_data(ocean_modelcreated_vars_dict["vectorsum"], modeldata) == fill(30.0, ocean_length)
    @test  PB.get_data(global_modelcreated_vars_dict["scalarsum"], modeldata)[] == 4.0

    @info "ocean constant concentrations"
    const_conc_data = PB.get_data(ocean_modelcreated_vars_dict["const_conc"], modeldata)
    @test const_conc_data.v             == fill(0.1, ocean_length)
    @test const_conc_data.v_moldelta    == fill(-0.2, ocean_length)

    @info "weighted mean"
    @test PB.get_data(ocean_modelcreated_vars_dict["const_conc_mean"], modeldata)[] == PB.IsotopeLinear(0.1, -0.2)

    @info "volume in range"
    @test PB.get_data(ocean_modelcreated_vars_dict["frac"], modeldata)[] == 1.0


    @info "test complete"
    # close(logfile)

end
