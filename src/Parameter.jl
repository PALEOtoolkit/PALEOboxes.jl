
import Preferences

"""
    Parameter{T, ParseFromString}

A reaction parameter of type `T`.
Create using short names `ParDouble`, `ParInt`, `ParBool`, `ParString`.
Read value as `<par>.v::T`, set with [`setvalue!`](@ref).

Parameters with `external=true` may be set from the Model-level Parameters list, if `name` is present in that list.

`ParseFromString` should usually be `Nothing`: a value of `Type T` is then required when calling [`setvalue!`](@ref).
If `ParseFromString` is not `Nothing`, then [`setvalue!`](@ref) will accept an `AbstractString` and call `Base.parse(ParseFromString, strvalue)`.
This allows eg an enum-valued Parameter to be defined by Parameter{EnumType, EnumType} and implementing
parse(EnumType, rawvalue::AbstractString)
"""
mutable struct Parameter{T, ParseFromString} <: AbstractParameter
    name::String
    description::String
    units::String
    "value"
    v::T    
    default_value::T
    allowed_values::Vector{T}
    frozen::Bool
    external::Bool
end

"""
    ParDouble(name,  defaultvalue::Float64; units="", description="", external=false)
"""
function ParDouble(
    name,  defaultvalue::Float64; 
    units="", description="", external=false
)
    return Parameter{Float64, Nothing}(
        name, description, units, defaultvalue, defaultvalue, Vector{Float64}(), false, external
    )
end

"""
    ParInt(name,  defaultvalue::Integer; units="", description="", allowed_values=Vector{Int}(), external=false)
"""
function ParInt(
    name,  defaultvalue::Integer;
    units="", description="", allowed_values=Vector{Int}(), external=false
)
    # splat allowed_values to allow Tuple etc
    return Parameter{Int, Nothing}(
        name, description, units, defaultvalue, defaultvalue, [allowed_values...], false, external
    )
end

"""
    ParEnum(name,  defaultvalue::T; units="", description="", allowed_values=Vector{T}(), external=false)

T can be an Enum (or any Type) that implements Base.parse.
"""
function ParEnum(
    name,  defaultvalue::T;
    units="", description="", allowed_values=Vector{T}(), external=false
) where {T}
    # splat allowed_values to allow Tuple etc
    return Parameter{T, T}(
        name, description, units, defaultvalue, defaultvalue, [allowed_values...], false, external
    )
end

"""
    ParType(::Type{T}, name,  defaultvalue; units="", description="", allowed_values=Vector{T}(), external=false)

T can be an abstract Type that implements Base.parse.
"""
function ParType(
    ::Type{T}, name,  defaultvalue::Type{D};
    units="", description="", allowed_values=Vector{Type}(), external=false
) where {T, D <: T}
    # splat allowed_values to allow Tuple etc
    return Parameter{Type, T}(
        name, description, units, defaultvalue, defaultvalue, [allowed_values...], false, external
    )
end

"""
    ParBool(name,  defaultvalue::Bool;  description="", external=false)
"""
function ParBool(
    name,  defaultvalue::Bool;
    description="", external=false
)
    return Parameter{Bool, Nothing}(
        name, description, "", defaultvalue, defaultvalue, Vector{Bool}(), false, external
    )
end

"""
    ParString(name,  defaultvalue::String;  description="", allowed_values=Vector{String}(), external=false)
"""
function ParString(
    name,  defaultvalue::String;
    description="", allowed_values=Vector{String}(), external=false
)
    # splat allowed_values to allow Tuple etc
    par = Parameter{String, Nothing}(
        name, description, "", defaultvalue, defaultvalue, [allowed_values...], false, external
    )
    setvalue!(par, defaultvalue) # catch attempt to set a defaultvalue not in allowed_values
    return par
end


"""
    VecParameter{T, ParseFromString}

A reaction parameter of type Vector{T}.
Create using short names `ParDoubleVec`, `ParStringVec`.
Values are `<par>.v::Vector{T}`.
Set using standard yaml syntax for a vector eg [1, 2, 3]

`ParseFromString` should usually be `false` (see [`Parameter`](@ref)).
"""
mutable struct VecParameter{T, ParseFromString} <: AbstractParameter
    name::String
    description::String
    units::String
    v::Vector{T}
    default_value::Vector{T}
    allowed_values::Vector{T}
    frozen::Bool
    external::Bool
end

"""
    ParDoubleVec(name, defaultvalue::Vector{Float64}; units="", description="", external=false)
"""
function ParDoubleVec(
    name, defaultvalue::Vector{Float64};
    units="", description="", external=false
)
    return VecParameter{Float64, Nothing}(
        name, description, units, defaultvalue, defaultvalue, Vector{Float64}(), false, external
    )
end

"""
    ParIntVec(name,  defaultvalue::Vector{Integer}; units="", description="", allowed_values=Vector{Int}(), external=false)
"""
function ParIntVec(
    name, defaultvalue::Vector{Int};
    units="", description="", allowed_values=Vector{Int}(), external=false
)
    # splat allowed_values to allow Tuple etc
    par = VecParameter{Int, Nothing}(
        name, description, units, defaultvalue, defaultvalue, [allowed_values...], false, external
    )
    setvalue!(par, defaultvalue) # catch attempt to set defaultvalue not in allowed_values
    return par
end

"""
    ParBoolVec(name,  defaultvalue::Vector{Bool};  description="", external=false)
"""
function ParBoolVec(
    name, defaultvalue::Vector{Bool};
    description="", external=false
)
    return VecParameter{Bool, Nothing}(
        name, description, "", defaultvalue, defaultvalue, Vector{Bool}(), false, external
    )
end

"""
    ParStringVec(name, defaultvalue::Vector{String};  description="", external=false)
"""
function ParStringVec(
    name, defaultvalue::Vector{String}=String[];
    description="", allowed_values=Vector{String}(), external=false
)
    # splat allowed_values to allow Tuple etc
    par = VecParameter{String, Nothing}(
        name, description, "", defaultvalue, defaultvalue, [allowed_values...], false, external
    )
    setvalue!(par, defaultvalue) # catch attempt to set defaultvalue not in allowed_values
    return par
end

"""
    VecVecParameter{T, ParseFromString}

A reaction parameter of type Vector{Vector{T}}.
Create using short names `ParDoubleVecVec`.
Values are `<par>.v::Vector{Vector{T}}`.
Set using standard yaml syntax for a vector of vectors eg [[1, 2, 3], [4, 5, 6]]

`ParseFromString` should usually be `false` (see [`Parameter`](@ref)).
"""
mutable struct VecVecParameter{T, ParseFromString} <: AbstractParameter
    name::String
    description::String
    units::String
    v::Vector{Vector{T}}
    default_value::Vector{Vector{T}}
    allowed_values::Vector{T} # Int, String only
    frozen::Bool
    external::Bool
end

"""
    ParDoubleVecVec(name, defaultvalue::Vector{Vector{Float64}}; units="", description="", external=false)
"""
function ParDoubleVecVec(
    name, defaultvalue::Vector{Vector{Float64}} = Vector{Vector{Float64}}();
    units="", description="", external=false
)
    par = VecVecParameter{Float64, Nothing}(
        name, description, units, defaultvalue, defaultvalue, Vector{Float64}(), false, external
    )
    return par
end

"""
    check_parameter_sum(parameter::VecParameter, ncells) -> sumok::Bool

Check whether `parameter` is of length `ncells` and correctly normalized to sum to 1.0
"""
function check_parameter_sum(parameter::VecParameter, ncells; tol=1e-3)
    sumok = true
    
    isnothing(ncells) || length(parameter.v) == ncells || 
        (sumok = false; @error "config error: length($(parameter.name)) != ncells")
    sumpar = sum(parameter.v)
    abs(sumpar - 1.0) < tol || 
        (sumok = false; @error "config error: sum($(parameter.name)) = $sumpar != 1 +/- $(tol)")

    return sumok
end

"""
    vecvecpar_matrix(par::VecVecParameter) -> Matrix

Convert `par.v` to a Matrix

# Examples:
 ```jldoctest; setup = :(import PALEOboxes)
julia> p = PB.ParDoubleVecVec("pdvv");

julia> PB.setvalue!(p, [[1.0, 2.0], [3.0, 4.0]]);

julia> PB.vecvecpar_matrix(p)
2×2 Array{Float64,2}:
 1.0  2.0
 3.0  4.0
 ```
"""
function vecvecpar_matrix(par::VecVecParameter{T}) where T
    nrows = length(par.v)
    if nrows == 0
        return Matrix{T}(undef, 0, 0)
    end
    ncols = length(par.v[1])
    m = Matrix{T}(undef, nrows, ncols)
    for (i, rowvals) in enumerate(par.v)
        length(rowvals) == ncols || error("par is not convertable to a matrix: ", par)
        m[i, :] .= rowvals
    end
    return m
end

"Prevent modification via properties"
function Base.setproperty!(par::AbstractParameter, s::Symbol, v)
    error("setproperty! attempt to set Parameter $(par.name).$s=$v (use setvalue! to set :v)")
end

# Default - no attempt to parse value from String
function _parsevalue(par::Union{Parameter{T, Nothing}, VecParameter{T, Nothing}, VecVecParameter{T, Nothing}}, value) where {T}
    return value
end

function _parsevalue(par::Union{Parameter{T, Nothing}, VecParameter{T, Nothing}, VecVecParameter{T, Nothing}}, value::AbstractString) where {T}
    return value
end

# any (Data)Type - any subtype of PT is acceptable
function _parsevalue(par::Union{Parameter{T, PT}, VecParameter{T, PT}, VecVecParameter{T, PT}}, value::Type{V}) where {T, PT, V <: PT}
    return value
end

# enum is a value of type T
function _parsevalue(par::Union{Parameter{T, T}, VecParameter{T, T}, VecVecParameter{T, T}}, value::T) where {T}
    return parse(T, value)
end

# attempt to parse value::PT from AbstractString
function _parsevalue(par::Union{Parameter{T, PT}, VecParameter{T, PT}, VecVecParameter{T, PT}}, value::AbstractString) where {T, PT}
    return parse(PT, value)
end

"""
    setvalue!(par::Parameter, value)

Set Parameter to `value`.

Optionally (if [`Parameter`](@ref) has Type parameter `ParseFromString != Nothing`) parse `value` from a String.
"""
function setvalue!(par::Parameter, rawvalue)
    value = _parsevalue(par, rawvalue)
    if !isempty(par.allowed_values) && !(value in par.allowed_values)
        error("setvalue! attempt to set Parameter $(par.name) to invalid value=$value (allowed_values=$(par.allowed_values))")
    end
    if par.frozen
        error("setvalue! Parameter $(par.name) can no longer be modified")
    end
    setfield!(par, :v, value)
    return nothing
end



function setvalue!(par::VecParameter, values)  
    if !isempty(par.allowed_values)
        for  rawv in values
            v = _parsevalue(par, rawv)
            if !(v in par.allowed_values)
                error("setvalue! attempt to set VecParameter $(par.name) to $values with invalid value(s) (allowed_values=$(par.allowed_values))")
            end
        end
    end
    if par.frozen
        error("setvalue! Parameter $(par.name) can no longer be modified")
    end
    setfield!(par, :v, values)
    return nothing
end

# allow an empty Vector{Any} (as that is what yaml provides for [])
function setvalue!(par::VecParameter{T, ParseFromString}, values::Vector{Any}) where{T, ParseFromString}
    if isempty(values)
        return setvalue!(par, T[])
    else
        return setvalue!(par, values) # will error
    end
end

function setvalue!(par::VecVecParameter, values)
    if !isempty(par.allowed_values)
        for  rawv in Iterators.Flatten(values)
            v = _parsevalue(par, rawv)
            if !(v in par.allowed_values)
                error("setvalue! attempt to set VecVecParameter $(par.name) to $values with invalid value(s) (allowed_values=$(par.allowed_values))")
            end
        end
    end
    if par.frozen
        error("setvalue! Parameter $(par.name) can no longer be modified")
    end
    setfield!(par, :v, values)
    return nothing
end

"""
    setvalueanddefault!(par::Parameter, value; freeze=false)

Set Parameter value and default to `value`.

Optionally (if [`Parameter`](@ref) has Type parameter `ParseFromString != Nothing`) parse `value` from a String.
"""
function setvalueanddefault!(par::Union{Parameter, VecParameter, VecVecParameter}, rawvalue; freeze=false)
    value = _parsevalue(par, rawvalue)
    setvalue!(par, value)
    setfield!(par, :default_value, value)
    if freeze
        setfrozen!(par)
    end
    return nothing
end

"""
    externalvalue(rawvalue::AbstractString, external_parameters) -> value

Replace `"external%somename"`` with `external_parameters["somename"]`
"""
function externalvalue(rawvalue::AbstractString, external_parameters)
    # substitute 'external%parname' with value of external parameter 'parname'
    if length(rawvalue) > 9 && findfirst("external%", rawvalue) == 1:9 
        externalkey = rawvalue[10:end]
        if haskey(external_parameters, externalkey)
            value = external_parameters[externalkey]
            @info "      expandvalue: $rawvalue -> $value"
        else
            error("      $(rawvalue) external parameter $(externalkey) not found")
        end
    else
        value = rawvalue        
    end
    
    return value
end

"Non-string values returned unmodified"
externalvalue(rawvalue, external_parameters) = rawvalue

"""
    substitutevalue(module, rawvalue::AbstractString) -> value


Substitute "\$MyLocalSetting\$" with the value of the MyLocalSetting key in the 
Julia LocalPreferences.toml file, from the section for `module` 
(for a Reaction Parameter value, this should be the Julia module of the Reaction containing the Parameter).
"""
function substitutevalue(mod::Module, rawvalue::AbstractString; dontsub=("\$fluxname\$",)) 
    value = rawvalue
    for m in eachmatch(r"\$\w+\$", rawvalue)
        if m.match in dontsub
            # dont substitute
            continue
        end
        
        # look up in LocalPreferences.toml
        subval = m.match[2:end-1] # omit $ $

        # if a Reaction is not part of a package (eg development code), use Preferences package instead to look for a key
        if isnothing(Base.PkgId(mod).uuid)
            @warn "substitutevalue $subval - Module $(mod) does not correspond to a loaded package, using a key from [PALEOboxes] instead"
            mod = PALEOboxes
        end

        Preferences.has_preference(mod, String(subval)) || 
            error("substitutevalue: key [$(Base.PkgId(mod).name)] $subval not found in LocalPreferences.toml")
        replacestr = Preferences.load_preference(mod, String(subval)) 
       
            
        value = replace(value, m.match => replacestr)
    end

    return value
end

"Return non-string values unmodified"
substitutevalue(mod::Module, rawvalue) = rawvalue


function setfrozen!(parameters::AbstractParameter...)
    for par in parameters
        setfield!(par, :frozen, true)
    end
end

"""
    ParametersTuple(parameters::AbstractParameter...) -> NamedTuple
    ParametersTuple(parameters) -> NamedTuple

Create a NamedTuple of Parameters.
"""
function ParametersTuple(parameters::AbstractParameter...)
    return NamedTuple{Tuple(Symbol(par.name) for par in parameters)}(parameters)
end
ParametersTuple(parameters) = ParametersTuple(parameters...)
