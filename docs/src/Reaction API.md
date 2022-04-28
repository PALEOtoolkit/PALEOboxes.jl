# Reaction API
```@meta
CurrentModule = PALEOboxes
```
API calls used by [`AbstractReaction`](@ref) implementations.

Implementations should create a subclass of [`AbstractReaction`](@ref) containing a [`ReactionBase`](@ref) and
[`Parameter`](@ref)s, register a reaction factory, and implement callbacks to define [`VariableReaction`](@ref)s and [`ReactionMethod`](@ref)s.

## Reaction struct
```@meta
CurrentModule = PALEOboxes
```
```@docs
AbstractReaction
ReactionBase
```

## Parameters
```@docs
AbstractParameter
Parameter
VecParameter
VecVecParameter
ParametersTuple
add_par
setvalue!
```

## Registering with PALEOboxes framework
```@docs
create_reaction(ReactionType::Type{<:AbstractReaction}, base::ReactionBase)

find_reaction
find_all_reactions
```

## Optional initialisation callbacks to define Domain Grids and array sizes
One Reaction per domain may implement [`set_model_geometry`](@ref) to create and attach
a grid, and [`set_data_dimension!`](@ref).
```@docs
set_model_geometry
set_data_dimension!
```

## Initialisation callbacks
```@docs
register_methods!
register_dynamic_methods!
```

### Defining [`VariableReaction`](@ref)s
```@docs
CreateVariableReaction
parse_variablereaction_namestr
set_attribute!
```

### Defining collections of VariableReactions
```@docs
AbstractVarList
VarList_single
VarList_namedtuple
VarList_tuple
VarList_vector
VarList_vvector
VarList_nothing
VarList_tuple_nothing
```

### Registering ReactionMethods
```@docs
ReactionMethod
add_method_setup!
add_method_initialize!
add_method_do!
```

## Predefined ReactionMethods

### Setup and initialization of Variables
```@docs
add_method_setup_initialvalue_vars_default!
add_method_initialize_zero_vars_default!
```

### Adding totals for Variables
```@docs
add_method_do_totals_default!
```

### Chemical reactions
```@docs
RateStoich
create_ratestoich_method
```

## Implementing Reaction method functions
Method functions should iterate over cells using the supplied `cellrange.indices`.

### Optimising loops over cells using explicit SIMD instructions
Reactions with simple loops over `cellindices` that implement time-consuming per-cell calculations 
may be optimised by using explicit SIMD instructions.
```@docs
SIMDutils.SIMDIter
```

## Internal details of Variable arrays accessor generation
[`VariableReaction`](@ref)s in the [`AbstractVarList`](@ref)s for a [`ReactionMethod`](@ref) are processed by [`create_accessor`](@ref) to supply `view`s on arrays as corresponding arguments to the [`ReactionMethod`](@ref) function.
```@docs
create_accessor
create_accessors
```
