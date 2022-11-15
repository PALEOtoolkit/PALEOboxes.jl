
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

    # test variable attribute access
    @test PB.get_variable_attribute(model, "global.A", :norm_value) == 3.193e17

    # allocate internal variables
    modeldata =  PB.create_modeldata(model)
    PB.allocate_variables!(model, modeldata, 1; hostdep=false)

    @test length( PB.get_unallocated_variables(global_domain, modeldata, 1)) == 4
    @test  PB.check_ready(global_domain, modeldata, throw_on_error=false) == false    

    # allocate arrays for host dependencies and set data pointers
    PB.allocate_variables!(model, modeldata, 1; hostdep=true)
    @test  PB.check_ready(global_domain, modeldata) == true

    # get all variable data arrays
    all_vars = PB.VariableAggregatorNamed(modeldata)
    @info "all_vars: $all_vars\n"
    all_data = all_vars.values

    # check state variables
    stateexplicit_vars, stateexplicit_sms_vars =
        PB.get_host_variables(global_domain, PB.VF_StateExplicit, match_deriv_suffix="_sms")
    @test length(stateexplicit_vars) == 2  # A and O
    # get global state variable aggregator
    global_stateexplicit_va = PB.VariableAggregator(stateexplicit_vars, fill(nothing, length(stateexplicit_vars)), modeldata, 1)
    
    stateexplicit_vars, stateexplicit_sms_vars =
        PB.get_host_variables(ocean_domain, PB.VF_StateExplicit, match_deriv_suffix="_sms")
    @test length(stateexplicit_vars) == 1
    # get ocean state variable aggregator
    ocean_stateexplicit_va = PB.VariableAggregator(stateexplicit_vars, fill(nothing, length(stateexplicit_vars)), modeldata, 1)

    # get host-dependent variables
    global_hostdep_vars_vec =  PB.get_variables(global_domain, hostdep=true)
    @test length(global_hostdep_vars_vec) == 4
 
    ocean_hostdep_vars_vec =  PB.get_variables(ocean_domain, hostdep=true)
    @test length(ocean_hostdep_vars_vec) == 2
 
    PB.initialize_reactiondata!(model, modeldata; create_dispatchlists_all=true)
    
    @info "dispatch_setup"
    PB.dispatch_setup(model, :setup, modeldata)
    PB.dispatch_setup(model, :norm_value, modeldata)
    PB.dispatch_setup(model, :initial_value, modeldata)

    @info "global stateexplicit variables:\n"
    @test all_data.global.A.v[]  == 3.193e18   # A
    @test all_data.global.A.v_moldelta[]  == 2.0*3.193e18    # A_moldelta
    # test access via flattened vector (as used eg by an ODE solver)
    global_A_indices = PB.get_indices(global_stateexplicit_va, "global.A")
    @test global_A_indices == 2:3
    stateexplicit_vector = PB.get_data(global_stateexplicit_va)
    @test stateexplicit_vector[global_A_indices] == [all_data.global.A.v[], all_data.global.A.v_moldelta[]]
    # test VariableAggregatorNamed created from VariableAggregator
    global_stateexplicit_data = PB.VariableAggregatorNamed(global_stateexplicit_va).values
    @test global_stateexplicit_data.global.A == all_data.global.A

    @info "global const variables:"
    @test all_data.global.ConstS[] == 1.0

    @info "ocean host-dependent variables:\n"
    @test all_data.ocean.T[1]           == 1.0*10.0    

    # take a time step
    dispatchlists = modeldata.dispatchlists_all
    PB.do_deriv(dispatchlists)

    @info "global model-created variables:\n"
    @test  all_data.global.A_norm[] == 10.0
    @test  all_data.global.A_delta[] == 2.0

    @info "ocean model-created variables:\n"
    @test  all_data.ocean.T_conc == fill(1.0, ocean_length)
    PB.do_deriv(dispatchlists) # repeat do_deriv to check total is reinitialized each time step
    @test  all_data.ocean.T_total[] == 1.0*10.0*ocean_length

    @info "sum variables:"
    @test  all_data.ocean.vectorsum == fill(30.0, ocean_length)
    @test  all_data.global.scalarsum[] == 4.0

    @info "ocean constant concentrations"
    const_conc_data = all_data.ocean.const_conc
    @test const_conc_data.v             == fill(0.1, ocean_length)
    @test const_conc_data.v_moldelta    == fill(-0.2, ocean_length)

    @info "weighted mean"
    @test all_data.ocean.const_conc_mean[] == PB.IsotopeLinear(0.1, -0.2)

    @info "volume in range"
    @test all_data.ocean.frac[] == 1.0


    @info "test complete"
    # close(logfile)

end
