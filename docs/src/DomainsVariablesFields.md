# Domains, Variables, Fields, and Data arrays.

```@meta
CurrentModule = PALEOboxes
```
A [`Model`](@ref) contains [`Domain`](@ref)s, each of which contain [Variables](@ref) defining [`Field`](@ref)s which contain `Data` arrays, and `Reactions` with [`ReactionMethod`](@ref)s that operate on the [`Field`](@ref)s to calculate model time evolution.

## Model
```@meta
CurrentModule = PALEOboxes
```
```@docs
Model
```

## Domains
```@meta
CurrentModule = PALEOboxes
```
```@docs
Domain
```

### Grids
```@meta
CurrentModule = PALEOboxes
```
```@docs
AbstractMesh

internal_size
cartesian_size
```
```@meta
CurrentModule = PALEOboxes.Grids
```
```@docs
set_subdomain!
get_subdomain

UnstructuredVectorGrid
UnstructuredColumnGrid
CartesianLinearGrid
CartesianArrayGrid
```

#### Subdomains
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
subdomain_view
subdomain_indices
```

#### Regions, Dimensions and Coordinates
```@meta
CurrentModule = PALEOboxes.Grids
```
```@docs
get_region
```

Grids may define name dimensions and attach coordinates for the convenience of output visualisation. Any coordinate information required by Reactions should be supplied as Variables.
```@meta
CurrentModule = PALEOboxes
```
```@docs
NamedDimension
FixedCoord
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


```@docs
VariableBase
VariableType
VariableDomain
VariableDomPropDep
VariableDomContribTarget
```

### Variable Attributes
```@meta
CurrentModule = PALEOboxes
```
```@docs
Attribute
StandardAttributes
```

```@example
import PALEOboxes as PB # hide
import DataFrames # hide
show( # hide
    DataFrames.sort!( # hide
        DataFrames.DataFrame((s=>getproperty.(PB.StandardAttributes, s) for s in (:name, :default_value, :required, :units, :description))...), # hide
        [DataFrames.order(:required, rev=true), :name] # hide
    ); # hide
    allcols=true, allrows=true, truncate=96, show_row_number=false #hide
) # hide
```

```@docs
get_attribute
VariableFunction
VariablePhase
```

TODO standard Variable names.

Use NetCDF CF standard names <https://cfconventions.org/Data/cf-standard-names/76/build/cf-standard-name-table.html> ?

## Fields
```@meta
CurrentModule = PALEOboxes
```
Field, data and function space are defined by Variable [`Attribute`](@ref)s in combination with an [`AbstractMesh`](@ref)
- `:field_data` and optionally `:data_dims` define the subtype of [`AbstractData`](@ref)
- `:space` defines the subtype of [`AbstractSpace`](@ref)

Examples:
- `:field_data` = [`ScalarData`](@ref), `:space` = [`ScalarSpace`](@ref) defines a Domain scalar (0D) quantity.
- `:field_data` = [`IsotopeLinear`](@ref), `:space` = [`ScalarSpace`](@ref) defines a Domain scalar (0D) isotope quantity.
- `:field_data` = [`ScalarData`](@ref), `:space` = [`CellSpace`](@ref) defines a per-cell quantity in a spatial Domain.
- `:field_data` = [`ArrayScalarData`](@ref), `:data_dims =("wgrid",)` `:space` = [`ScalarSpace`](@ref) defines a Domain scalar (0D) quantity on a wavelength grid, a Vector of length given by the value of the Domain "wgrid" data dimension (see [`set_data_dimension!`](@ref))

```@docs
Field
get_field
wrap_field
```

## Spaces
```@meta
CurrentModule = PALEOboxes
```
```@docs
AbstractSpace
ScalarSpace
CellSpace
ColumnSpace
```

## Data types
```@meta
CurrentModule = PALEOboxes
```
```@docs
AbstractData
UndefinedData
ScalarData
ArrayScalarData


allocate_values
check_values
zero_values!
dof_values
get_values_output


init_values!
copyfieldto!
copytofield!
add_field!
add_field_vec!

num_components(field_data::Type{<:AbstractData})
get_components
```

### Isotopes
```@meta
CurrentModule = PALEOboxes
```
```@docs
AbstractIsotopeScalar

IsotopeLinear

isotope_totaldelta
get_total
get_delta
```






