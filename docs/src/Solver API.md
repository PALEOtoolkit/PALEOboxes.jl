# Solver API

```@meta
CurrentModule = PALEOboxes
```

A [`Model`](@ref) contains [`Domain`](@ref)s, each of which contain [Variables](@ref) defining [`Field`](@ref)s, and [Reactions](@ref) with [`ReactionMethod`](@ref)s that operate on the [`Field`](@ref)s to calculate model time evolution.


## Initialization
```@meta
CurrentModule = PALEOboxes
```
```@docs
create_model_from_config(config_file::AbstractString, configmodel::AbstractString)
create_modeldata(model::Model, vectype::DataType=Array{Float64,1})
ModelData
allocate_variables!

check_ready(model::Model, modeldata::AbstractModelData; throw_on_error=true)
check_configuration(model::Model)
initialize_reactiondata!(model::Model, modeldata::AbstractModelData)
dispatch_setup
create_dispatch_methodlists
ReactionMethodDispatchList
```

## Attaching numerical solvers
High-level access to aggregated collections of state Variables and derivatives (see [Accessing model objects](@ref) for low-level access).
```@docs
SolverView
VariableAggregator
VariableAggregator(vars, cellranges, modeldata)
create_solver_view
set_default_solver_view!
copy_norm!
set_statevar!
get_statevar_sms!
```

## Defining CellRanges
```@meta
CurrentModule = PALEOboxes
```
```@docs
AbstractCellRange
CellRange
CellRangeColumns
create_default_cellrange
```
```@meta
CurrentModule = PALEOboxes.Grids
```
```@docs
create_default_cellrange
cellrange_cartesiantile
```

## Main loop
```@meta
CurrentModule = PALEOboxes
```
```@docs
set_tforce!
do_deriv
dispatch_methodlist
```


## Diagnostics
```@meta
CurrentModule = PALEOboxes
```
```@docs
show_methods_setup
show_methods_initialize
show_methods_do
show_variables
show_links
show_parameters
```

## Accessing model objects

### [`Model`](@ref)
```@docs
get_domain
get_reaction(model::Model, domainname, reactionname)
```

### [`Domain`](@ref)s
```@docs
get_variable
get_variables
get_host_variables
get_unallocated_variables
get_reaction(domain::Domain, reactname::String)
```

### VariableDomain
Low-level access to individual [Variables](@ref).
```@meta
CurrentModule = PALEOboxes
```
```@docs
set_data!
get_data
get_data_output
```