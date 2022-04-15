# PALEOboxes

```@meta
CurrentModule = PALEOboxes
```

PALEOboxes provides a coupler for model components.

A [`Model`](@ref) contains [`Domain`](@ref)s, each of which contain [Variables](@ref) defining [`Field`](@ref)s, and [Reactions](@ref) with [`ReactionMethod`](@ref)s that operate on the [`Field`](@ref)s to calculate model time evolution.

## Model
```@docs
Model
get_domain
get_reaction(model::Model, domainname, reactionname)
```
### Initialization
```@docs
create_model_from_config(config_file::AbstractString, configmodel::AbstractString)
create_modeldata(model::Model, vectype::DataType=Array{Float64,1})
ModelData
allocate_variables!
set_default_solver_view!
check_ready(model::Model, modeldata::AbstractModelData; throw_on_error=true)
check_configuration(model::Model)
initialize_reactiondata!(model::Model, modeldata::AbstractModelData)
dispatch_setup
create_dispatch_methodlists
ReactionMethodDispatchList
```

### Attaching numerical solvers
```@docs
SolverView
VariableAggregator
VariableAggregator(vars, cellranges, modeldata)
create_solver_view
set_statevar!
get_statevar_sms!
```

### Main loop
```@docs
set_tforce!
do_deriv
dispatch_methodlist
```

### Diagnostics
```@meta
CurrentModule = PALEOboxes
```
```@docs
show_methods_setup
show_methods_initialize
show_methods_do
show_variables
show_parameters
```



## Domains
```@docs
Domain
set_data_dimension!
get_variable
get_variables
get_host_variables
get_unallocated_variables
get_reaction(domain::Domain, reactname::String)
```

### Grids
```@meta
CurrentModule = PALEOboxes
```
```@docs
AbstractMesh
```
```@meta
CurrentModule = PALEOboxes.Grids
```
```@docs
UnstructuredVectorGrid
UnstructuredColumnGrid
CartesianLinearGrid
CartesianArrayGrid
```

### CellRange
```@meta
CurrentModule = PALEOboxes
```
```@docs
AbstractCellRange
CellRange
CellRangeColumns
```
```@meta
CurrentModule = PALEOboxes.Grids
```
```@docs
create_default_cellrange
cellrange_cartesiantile
```

### Subdomains
```@meta
CurrentModule = PALEOboxes
```
```@docs
AbstractSubdomain
```
```@meta
CurrentModule = PALEOboxes.Grids
```
```@docs
BoundarySubdomain
InteriorSubdomain
set_subdomain!
subdomain_view
subdomain_indices
```


## Variables

```@meta
CurrentModule = PALEOboxes
```
PALEO `Variables` exist in [`Domain`](@ref)s and represent biogeochemical or other quantities.
They are defined by `Reactions` as [`VariableReaction`](@ref)s which are then linked to create [`VariableDomain`](@ref)s,
and contain [`Field`](@ref)s which represent data (a subtype of [`AbstractData`](@ref)) in a function space (a subtype of [`AbstractSpace`](@ref)) with dimensions defined by the Domain grid (a subtype of [`AbstractMesh`](@ref))

Dataflow dependency between `Variables` is represented by two pairings:
- [`VariableReaction`](@ref)s of [`VariableType`](@ref) `VT_ReactProperty` and `VT_ReactDependency` which are linked to create [`VariableDomPropDep`](@ref)s
- [`VariableReaction`](@ref)s of [`VariableType`](@ref) `VT_ReactTarget` and `VT_ReactContributor` which are linked to create  [`VariableDomContribTarget`](@ref)s.

Field, data and function space are defined by Variable [`Attribute`](@ref)s in combination with an [`AbstractMesh`](@ref)
- `:field_data` and optionally `:data_dims` define the subtype of [`AbstractData`](@ref)
- `:space` defines the subtype of [`AbstractSpace`](@ref)

Examples:
- `:field_data` = [`ScalarData`](@ref), `:space` = [`ScalarSpace`](@ref) defines a Domain scalar (0D) quantity.
- `:field_data` = [`IsotopeLinear`](@ref), `:space` = [`ScalarSpace`](@ref) defines a Domain scalar (0D) isotope quantity.
- `:field_data` = [`ScalarData`](@ref), `:space` = [`CellSpace`](@ref) defines a per-cell quantity in a spatial Domain.
- `:field_data` = [`ArrayScalarData`](@ref), `:data_dims =("wgrid",)` `:space` = [`ScalarSpace`](@ref) defines a Domain scalar (0D) quantity on a wavelength grid, a Vector of length given by the value of the Domain "wgrid" data dimension (see [`set_data_dimension!`](@ref))

[`FieldArray`](@ref) provides a generic array type with named dimensions [`NamedDimension`](@ref)s and optional coordinates [`FixedCoord`](@ref) for processing of model output.

```@docs
VariableBase
VariableType
```

TODO standard Variable names.

Use NetCDF CF standard names <https://cfconventions.org/Data/cf-standard-names/76/build/cf-standard-name-table.html> ?

### Data types
```@docs
AbstractData
UndefinedData
ScalarData
ArrayScalarData
```

#### Isotopes
```@docs
AbstractIsotopeScalar

IsotopeLinear

isotope_totaldelta
get_total
get_delta
```

### Spaces
```@docs
AbstractSpace
ScalarSpace
CellSpace
```

### Fields
```@docs
Field
FieldRecord
get_field
```

### FieldArray
```@docs
FieldArray
get_array
FixedCoord
NamedDimension
```

### Variable Attributes
```@docs
Attribute
StandardAttributes
get_attribute
num_components
VariableFunction
VariablePhase
```

### VariableDomain
```@meta
CurrentModule = PALEOboxes
```
```@docs
VariableDomain
VariableDomPropDep
VariableDomContribTarget

set_data!
get_data
get_data_output
```

## Reactions

```@meta
CurrentModule = PALEOboxes
```
```@docs
AbstractReaction
ReactionBase
```
### Registering with PALEOboxes framework
```@docs
add_reaction_factory
reaction_factory
```

### Parameters
```@docs
AbstractParameter
Parameter
VecParameter
VecVecParameter
ParametersTuple
add_par
setvalue!
```

### Initialisation callbacks
```@docs
set_model_geometry
register_methods!
register_dynamic_methods!
```

### VariableReactions
```@docs
VariableReaction
CreateVariableReaction
parse_variablereaction_namestr
set_attribute!
Reservoirs.ReservoirLinksVector
```

### ReactionMethods
```@docs
ReactionMethod
add_method_setup!
add_method_initialize!
add_method_do!
```
#### Defining collections of Variables
```@docs
AbstractVarList
VarList_single
VarList_namedtuple
VarList_tuple
VarList_vector
VarList_vvector
VarList_nothing
VarList_tuple_nothing
create_accessors
create_accessor
```

### Optimising loops over cells using explicit SIMD instructions
Reactions with simple loops over `cellindices` that implement time-consuming per-cell calculations 
may be optimised by using explicit SIMD instructions.
```@docs
SIMDutils.SIMDIter
```

### Predefined ReactionMethods

#### Setup and initialization of Variables
```@docs
add_method_setup_initialvalue_vars_default!
add_method_initialize_zero_vars_default!
```

#### Adding totals for Variables
```@docs
add_method_do_totals_default!
```

#### Chemical reactions
```@docs
RateStoich
create_ratestoich_method
```

## Reservoirs
```@meta
CurrentModule = PALEOboxes.Reservoirs
```
```@docs
ReactionReservoirScalar
ReactionReservoir
ReactionReservoirWellMixed
ReactionReservoirConst
ReactionReservoirForced
ReactionConst
```

## Variable Statistics
```@meta
CurrentModule = PALEOboxes.VariableStats
```
```@docs
ReactionSum
ReactionWeightedMean
ReactionAreaVolumeValInRange
```

## Fluxes
```@meta
CurrentModule = PALEOboxes
```
```@docs
Fluxes
```
```@meta
CurrentModule = PALEOboxes.Fluxes
```
```@docs
ReactionFluxTarget
ReactionFluxTransfer
```

### Contributing to flux couplers from a Reaction
```@docs
FluxContrib
```