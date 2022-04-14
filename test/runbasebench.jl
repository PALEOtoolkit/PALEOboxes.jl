import PALEOboxes as PB

using BenchmarkTools

include("ReactionPaleoMockModule.jl")

model = PB.create_model_from_config("configbase.yaml", "model1")

ocean_domain = PB.get_domain(model, "ocean")
ocean_domain_size = 1000
# domain_size = 10
ocean_domain.grid = PB.Grids.UnstructuredVectorGrid(ncells=ocean_domain_size)

modeldata = PB.create_modeldata(model)

PB.allocate_variables!(model, modeldata)

PB.set_default_solver_view!(model, modeldata) 

PB.set_data!(PB.get_variable(ocean_domain, "scalar_dep"), modeldata, [42.0])

PB.check_ready(model, modeldata)

PB.initialize_reactiondata!(model, modeldata)

PB.dispatch_setup(model, :initial_value, modeldata)

dispatchlists = modeldata.dispatchlists_all

println("@btime  PB.dispatch_list initialize")
@btime PB.dispatch_methodlist($dispatchlists.list_initialize)

println("@btime  PB.dispatch_list do")
@btime PB.dispatch_methodlist($dispatchlists.list_do)

println("@btime  PB.do_deriv")
@btime PB.do_deriv($dispatchlists)

