# Reaction API
```@meta
CurrentModule = PALEOboxes
```
API calls used by [`AbstractReaction`](@ref) implementations.

Implementations should create a subclass of [`AbstractReaction`](@ref) containing a [`ReactionBase`](@ref) and
[`Parameter`](@ref)s, optionally provide a [`create_reaction`](@ref), and implement callbacks to define [`VariableReaction`](@ref)s and register [`ReactionMethod`](@ref)s.

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

## Registering Reactions with the PALEOboxes framework

All subtypes of [`AbstractReaction`](@ref) that are present in loaded modules are
available to the PALEO framework.  Available Reactions can be listed with [`find_reaction`](@ref)
and [`find_all_reactions`](@ref).  The default [`create_reaction`](@ref) is called to create Reactions
when the model is created (this method can be overridden if needed).

```@docs
create_reaction(ReactionType::Type{<:AbstractReaction}, base::ReactionBase)

find_reaction
find_all_reactions
```

## Defining Domain Grids and array sizes
```@docs
set_model_geometry
set_data_dimension!
```

## Creating and registering Reaction methods
All Reactions should implement [`register_methods!`](@ref), and may optionally implement [`register_dynamic_methods!`](@ref).
```@docs
register_methods!
register_dynamic_methods!
```

These methods should then define one or more [`ReactionMethod`](@ref)s, which requires:
- Defining collections of [`VariableReaction`](@ref)s (see [Defining VariableReactions](@ref),  [Defining collections of VariableReactions](@ref)).
- Implementing a function to iterate over model cells and calculate the required fluxes etc (see [Implementing method functions](@ref))
- Adding methods (see [Adding ReactionMethods](@ref)).

In addition it is possible to add [Predefined ReactionMethods](@ref) for some common operations (Variable initialisation, calculating totals, etc).

### Defining VariableReactions
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

### Implementing method functions

Reaction method functions should iterate over the cells in the Domain supplied by `cellrange` argument and calculate appropriate biogeochemical fluxes etc (which may include the model time derivative and any intermediate or diagnostic output).

#### Iterating over cells

The simplest case is a method function that iterates over individual cells, with skeleton form:

    function do_something_cellwise(m::PB.AbstractReactionMethod, pars, (vars, ), cellrange::PB.AbstractCellRange, deltat)

        @inbounds for i in cellrange.indices
            vars.A[i]  = something*vars.B[i]*vars.C[i]  # in general A is some function of B, C, etc
            # etc
        end

        return nothing
    end

#### Iterating over cells in columns

If necessary (eg to calculate vertical transport), provided the model grid and cellrange allow,
it is possible to iterate over columns and then cells within columns (in order from top to bottom):

    function do_something_columnwise(m::PB.AbstractReactionMethod, pars, (vars, ), cellrange::PB.AbstractCellRange, deltat)

        @inbounds for (icol, colindices) in cellrange.columns
            accum = zero(vars.A[first(colindices)]) # accumulator of appropriate type
            for i in colindices  # order is top to bottom
                accum += vars.A[i]
                vars.C[i] = accum  # C = sum of A in cells above                 
                # etc
            end

            vars.floor_C[icol] = vars.C[last(colindices)] # assumes model has a floor domain with one floor cell per column in the interior domain
        end

        return nothing
    end

Iteration from bottom to top within a column can be implemented using `Iterators.reverse`, eg

    function do_something_columnwise(m::PB.AbstractReactionMethod, pars, (vars, ), cellrange::PB.AbstractCellRange, deltat)
        @inbounds for (icol, colindices) in cellrange.columns
            colreverse = Iterators.reverse(colindices)
            for i in colreverse  # order is bottom to top
                # etc
            end
        end

        return nothing
    end

!!! note
    The method function shouldn't make any assumptions about `colindices` other than that it is a list of indices ordered from top to bottom in a column.  Depending on the grid in use, the indices may not be contiguous, and may not be integers.

!!! note
    The example above made the additional assumption that a `floor` domain had been defined (containing Variable `floor_C`) with one floor cell per column. This is determined by the model configuration, and is not true in general.

In rare cases where it is necessary to operate on a Vector representing a quantity for the whole column (rather than just iterate through it), this can be implemented using `view`, eg

    function do_something_columnwise(m::PB.AbstractReactionMethod, pars, (vars, ), cellrange::PB.AbstractCellRange, deltat)
        @inbounds for (icol, colindices) in cellrange.columns
            A_col = view(vars.A, colindices)  # A_col is an AbstractVector with contiguous indices 1:length(colindices)
            B_col = view(vars.B, colindices)  # B_col is an AbstractVector with contiguous indices 1:length(colindices)
            
            # do something that needs a vector of cells for a whole column
        end

        return nothing
    end


#### Optimising loops over cells using explicit SIMD instructions
Reactions with simple loops over `cellindices` that implement time-consuming per-cell calculations 
may be optimised by using explicit SIMD instructions.
```@docs
SIMDutils.SIMDIter
```

### Adding ReactionMethods
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

## Internal details of Variable arrays accessor generation
[`VariableReaction`](@ref)s in the [`AbstractVarList`](@ref)s for a [`ReactionMethod`](@ref) are processed by [`create_accessor`](@ref) to supply `view`s on arrays as corresponding arguments to the [`ReactionMethod`](@ref) function.
```@docs
create_accessor
create_accessors
```
