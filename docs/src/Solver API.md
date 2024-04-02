# Solver API

```@meta
CurrentModule = PALEOboxes
```

A [`Model`](@ref) contains [`Domain`](@ref)s, each of which contain [Variables](@ref) defining [`Field`](@ref)s, and `Reactions` with [`ReactionMethod`](@ref)s that operate on the [`Field`](@ref)s to calculate model time evolution.


## Create and initialise
```@meta
CurrentModule = PALEOboxes
```
```@docs
create_model_from_config(config_file::AbstractString, configmodel::AbstractString)
create_modeldata(model::Model, vectype::DataType=Array{Float64,1})
ModelData
add_arrays_data!
push_arrays_data!
allocate_variables!

check_ready(model::Model, modeldata::AbstractModelData; throw_on_error=true)
check_configuration(model::Model)
initialize_reactiondata!(model::Model, modeldata::AbstractModelData)
dispatch_setup
create_dispatch_methodlists
ReactionMethodDispatchList
```

## Attaching numerical solvers
High-level access to aggregated collections of Variables is provided by [`VariableAggregator`](@ref) and [`VariableAggregatorNamed`](@ref) (see [Accessing model objects](@ref) for low-level access).
```@docs
VariableAggregator
get_indices
copyto!(dest::VariableAggregator, src::AbstractVector; sof::Int=1)
copyto!(dest::AbstractVector, src::VariableAggregator; dof::Int=1)

VariableAggregatorNamed
```

Aggregated collections of a subset of Parameters as a flattened Vector (eg for sensitivity studies) is provided by [`ParameterAggregator`](@ref):
```@docs
ParameterAggregator
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
set_parameter_value!
get_parameter_value
set_variable_attribute!
get_variable_attribute

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
