

#####################################################
# Attribute accessors
#####################################################

"""
    set_attribute!(var::VariableBase, name::Symbol, value, allow_create=false) -> var

Set Variable attribute.
"""
function set_attribute!(var::VariableBase, name::Symbol, value; allow_create=false)
    _set_attribute!(var.attributes, name, value, allow_create)
    return var
end

"""
    has_attribute(var::VariableBase, name::Symbol)

Test Variable attribute present.
"""
has_attribute(var::VariableBase, name::Symbol) = haskey(var.attributes, name)

"""
    get_attribute(var::VariableBase, name::Symbol, missing_value=missing) -> value

Get Variable attribute.
"""
function get_attribute(var::VariableBase, name::Symbol, missing_value=missing)
    return get(var.attributes, name, missing_value)
end

"""
    num_components(var::VariableBase) -> Int

Number of components in `var`.
"""
num_components(var::VariableBase) = num_components(get_attribute(var, :field_data))


###################################
# enum-valued VariableAttributes
####################################

"""
    @enum VariableFunction

Allowed values of `:vfunction` Variable [`Attribute`](@ref), 
defining the Variable function to the host ODE or DAE solver.

Explicit ODE problems with dS/dt = F(S) consist of pairs of S::VF_StateExplicit, F::VF_Deriv Variables.

An implicit ODE problem with dU/dt = F(U) where Total variables U are functions U(S) of State variables S
will consist of pairs of U::VF_Total and F::VF_Deriv Variables, and also the same number of S::VF_State (in 
no particular order).

Algebraic constraints  C(S) = 0 include variables C::VF_Constraint and the same number of S::VF_State,
with no corresponding VF_Deriv.

Not all solvers support all combinations.
"""
@enum VariableFunction::Cint begin
    VF_Undefined        = 0
    VF_StateExplicit    = 1
    VF_Total            = 2
    VF_Constraint       = 3
    VF_State            = 4
    VF_Deriv            = 5
end

"parse eg \"VF_Deriv\" as Enum VF_Deriv"
function Base.parse(::Type{VariableFunction}, str::AbstractString)   
    return getproperty(@__MODULE__, Symbol(str))::VariableFunction
end


"""
    @enum VariablePhase

Allowed values of `:vphase` Variable [`Attribute`](@ref), defining the component phase this Variable belongs to 
for multiphase cells.
"""
@enum VariablePhase::Cint begin
    VP_Undefined    = 0
    VP_Solute       = 1
    VP_Solid        = 2
end

"parse eg \"VP_Solute\" as Enum VP_Solute"
function Base.parse(::Type{VariablePhase}, str::AbstractString)   
    return getproperty(@__MODULE__, Symbol(str))::VariablePhase
end

###################################################
# Standard Attributes
#####################################################

"""
    Attribute{T}

Definition for Variable attribute `name`. Defines a data type `T`, a `default_value`, `required` (`true` if always present ie
added with `default_value` when Variable is created), `units`, and an optional `description`.

Note that Variable attributes are stored as a per-Variable `Dict(name => value)`, these definitions are only used to provide defaults,
check types, and provide descriptive metadata.

`ParseFromString` should usually be `false`: a value of `Type T` is then required when calling [`set_attribute!`](@ref).
If `ParseFromString` is `true`, then [`set_attribute!`](@ref) will accept an `AbstractString` and call `Base.parse(T, strvalue)`
to convert to `T`. This allows eg an enum-valued Attribute to be defined by Attribute{EnumType, true} and implementing
parse(EnumType, rawvalue::AbstractString)
"""
struct Attribute{T, ParseFromString}
    name::Symbol
    default_value::T
    required::Bool
    units::String
    description::String
end

attributeType(::Attribute{T}) where{T} = T

"""
    StandardAttributes

List of standard Variable attributes.

Some of these follow netCDF COARDS/CF conventions:

COARDS:<https://ferret.pmel.noaa.gov/Ferret/documentation/coards-netcdf-conventions>
- units   (where possible should follow the Unidata udunits package <https://docs.unidata.ucar.edu/udunits/current/>)
- long_name

CF conventions: <https://cfconventions.org/cf-conventions/cf-conventions.html#_description_of_the_data>
- standard_name   <http://cfconventions.org/Data/cf-standard-names/current/src/cf-standard-name-table.xml>

""" 
const StandardAttributes = [
    #                                   name                    default_value    required    units       description    

    Attribute{Type, AbstractData}(        :field_data,           UndefinedData,   true,       "",         "AbstractData type Variable contains")
    Attribute{Tuple{Vararg{String}}, Tuple{Vararg{String}}}(
                                          :data_dims,            (),              true,       "",         "Tuple of variable data dimension names, or empty for a scalar")
    Attribute{Type, AbstractSpace}(       :space,                CellSpace,       true,       "",         "function space Variable is defined on")
    Attribute{String, Nothing}(           :mesh,                 "default",       true,       "",         "grid mesh on which Variable is defined (empty for Domain spatial scalar)")
    Attribute{Bool, Nothing}(             :check_length,         true,            false,      "",         "true to check length matches length of linked VariableDomain")
    Attribute{Bool, Nothing}(             :is_constant,          false,           true,      "",          "true if variable is not changed after initialisation")
    Attribute{VariableFunction, VariableFunction}(
                                          :vfunction,             VF_Undefined,   true,       "",         "host function (to label state variables etc)")
    Attribute{Vector{Int}, Nothing}(      :operatorID,            Int[],          false,      "",         "Reaction operatorIDs that modify this Variable")
    Attribute{VariablePhase, VariablePhase}(
                                          :vphase,                VP_Undefined,   false,      "",         "phase for concentrations in multiphase cells")
    Attribute{String, Nothing}(           :totalname,             "",             false,      "",         "total Variable name for this species")
    Attribute{String, Nothing}(           :safe_name,             "",             false,      "",         "optional short or escaped name for compatibility with other software")
    Attribute{String, Nothing}(           :long_name,             "",             false,      "",         "netcdf long descriptive name")
    Attribute{String, Nothing}(           :units,                 "",             true,       "",         "where possible should follow netcdf conventions")
    Attribute{String, Nothing}(           :description,           "",             true,       "",         "human-readable description")
    Attribute{String, Nothing}(           :standard_name,         "",             false,      "",         "netcdf CF conventions standard name")
    Attribute{Bool, Nothing}(             :advect,                false,          false,      "",         "true to apply advective transport to this tracer")
    Attribute{Float64, Nothing}(          :advect_zmin,           0.0,            false,      "m",        "minimum height for transport")
    # Attribute{Bool, Nothing}(            :optional,             false,          true,       "",         "")
    Attribute{Bool, Nothing}(             :initialize_to_zero,    false,          true,       "",         "request initialize to zero at start of each timestep.")
    Attribute{Float64, Nothing}(          :vertical_movement,     0.0,            false,      "m d-1",    "vertical advective transport (+ve upwards) in column")
    Attribute{String, Nothing}(           :diffusivity_speciesname, "",           false,      "",         "species name to define diffusivity")
    Attribute{Union{Float64,Vector{Float64}}, Nothing}(
                                          :initial_value,         0.0,            true,       "",         "initial value to be applied eg to state or constant variable")
    Attribute{Union{Float64,Vector{Float64}}, Nothing}(
                                          :norm_value,            1.0,            true,       "",         "normalisation value to be passed through to solver and optionally provided as model variable")
    Attribute{Float64, Nothing}(          :initial_delta,         0.0,            true,       "per mil",  "initial value for isotope variable to be applied eg to state or constant variable")
    Attribute{Float64, Nothing}(          :specific_light_extinction, 0.0,        false,      "m^2 mol-1","wavelength-independent specific extinction in water column")
    Attribute{Union{Float64,Missing}, Nothing}(
                                          :deposition_velocity,   missing,        false,      "cm s-1",   "surface deposition velocity for atmospheric tracer")
    Attribute{Union{Float64,Missing}, Nothing}(
                                          :rainout,               missing,        false,      "",         "normalized rainout rate for atmospheric tracer")
]


"""
    default_variable_attributes()

Required attributes that should be present in all `VariableReaction`
"""
function default_variable_attributes()
    return Dict(a.name => a.default_value for a in StandardAttributes if a.required)
end


"""
    is_standard_attribute(attribute::Symbol) -> Bool

```jldoctest; setup = :(import PALEOboxes)
julia> PALEOboxes.is_standard_attribute(:data_dims)
true
julia> PALEOboxes.is_standard_attribute(:foo)
false
```
"""
function is_standard_attribute(attribute::Symbol) 
    return !isnothing(findfirst(a->a.name==attribute, StandardAttributes))
end

###################################
# Internal implementation
##################################

"Set attribute, with type conversion and checking"
function _set_attribute!(attributes::Dict{Symbol, Any}, name::Symbol, value, allow_create)
    if haskey(attributes, name) || allow_create
        attributes[name] = _convert_standard_attribute_type(name, value)        
    else
        error("invalid attempt to create new attribute ",name,"=",value)
    end
    return nothing
end

"""
NB: conversion of :initial_value or :norm_value from an integer fails, eg
julia> PALEOboxes._convert_standard_attribute_type(:initial_value, 1) 
as
julia> convert(Union{Float64,Vector{Float64}}, 1)
fails.

```jldoctest; setup = :(import PALEOboxes)
julia> PALEOboxes._convert_standard_attribute_type(:initial_delta, 1)
1.0

julia> PALEOboxes._convert_standard_attribute_type(:my_new_attribute, 1)
1

julia> PALEOboxes._convert_standard_attribute_type(:vphase, "VP_Solute")
VP_Solute::VariablePhase = 1
```
"""
function _convert_standard_attribute_type(attribute::Symbol, value)
    a_idx = findfirst(a->a.name==attribute, StandardAttributes)
    if isnothing(a_idx)
        return value # no type conversion
    else
        stdatt = StandardAttributes[a_idx]
        return convert(attributeType(stdatt), _parsevalue(stdatt, value))
    end
end

# Default - no attempt to parse value from String
function _parsevalue(attribute::Attribute{T, Nothing}, value) where {T}
    return value
end

function _parsevalue(attribute::Attribute{T, Nothing}, value::AbstractString) where {T}
    return value
end

# any (Data)Type - any subtype of PT is acceptable
function _parsevalue(attribute::Attribute{T, PT}, value::Type{V}) where {T, PT, V <: PT}
    return value
end

# enum is a value of type T
function _parsevalue(attribute::Attribute{T, T}, value::T) where {T}
    return value
end

# attempt to parse value from AbstractString into type PT
function _parsevalue(attribute::Attribute{T, PT}, value::AbstractString) where {T, PT}
    return parse(PT, value)
end



