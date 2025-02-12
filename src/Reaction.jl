
"""
    AbstractReaction

Abstract base Type for Reactions.

# Implementation
Derived types should include a field `base::`[`ReactionBase`](@ref), and usually a [`ParametersTuple`](@ref), eg

    Base.@kwdef mutable struct ReactionHello{P} <: PB.AbstractReaction
        base::PB.ReactionBase

        pars::P = PB.ParametersTuple(
            PB.ParDouble("my_par", 42.0, units="yr", 
                description="an example of a Float64-valued scalar parameter called 'my_par'"),
        )

        some_additional_field::Float64   # additional fields to eg cache data read from disk etc
    end

Derived types should implement [`register_methods!`](@ref), and may optionally implement [`create_reaction`](@ref), 
[`set_model_geometry`](@ref), [`check_configuration`](@ref), [`register_dynamic_methods!`](@ref).

Methods should be registered using [`add_method_setup!`](@ref), [`add_method_initialize!`](@ref), [`add_method_do!`](@ref).
"""
AbstractReaction

"""
    ReactionBase

Base Type for a biogeochemical Reaction.

# Implementation
Include as field `base` in types derived from [`AbstractReaction`](@ref)
"""
Base.@kwdef mutable struct ReactionBase
    "domain containing this Reaction"
    domain::Union{Nothing, Domain}                  = nothing
    "operator ID (to allow operator splitting)"
    operatorID::Vector{Int}                         = [1]
    "name of this instance (from .yaml config)"         
    name::String
    "classname   (as defined for .yaml config)"                                                 
    classname::String
    "external parameters and values supplied from Model or Domain"
    external_parameters::Dict{String, Any}          = Dict{String, Any}()

    methods_setup::Vector{AbstractReactionMethod}   = Vector{AbstractReactionMethod}()
    methods_initialize::Vector{AbstractReactionMethod}=Vector{AbstractReactionMethod}()
    methods_do::Vector{AbstractReactionMethod}      = Vector{AbstractReactionMethod}()

    # temporary storage for Variable configuration (read from .yaml config and applied after Variables linked)
    _conf_variable_links = nothing
    _conf_variable_attributes = nothing
end

struct NoReaction <: AbstractReaction
    base::ReactionBase
    NoReaction() = new(ReactionBase(name="", classname="NoReaction"))
end

# helps type stability: even though we are using @nospecialize(AbstractReaction) everywhere, 
# we still know the type of base
base(react::AbstractReaction) = getfield(react, :base)::ReactionBase

"""
    Base.getproperty(react::AbstractReaction, s::Symbol)

Forward to `react.base::ReactionBase` to define additional properties.
"""
function Base.getproperty(react::AbstractReaction, s::Symbol)
    if s == :base
        return base(react)
    elseif s == :name
        return base(react).name
    elseif s == :classname
        return base(react).classname
    elseif s == :domain
        return base(react).domain
    elseif s == :operatorID
        return base(react).operatorID
    elseif s == :external_parameters
        return base(react).external_parameters
    elseif s == :methods_setup
        return base(react).methods_setup
    elseif s == :methods_initialize
        return base(react).methods_initialize
    elseif s == :methods_do
        return base(react).methods_do
    else
        return getfield(react, s)
    end
end

##########################################################
# Query methods for Parameters, ReactionMethods, Variables
##########################################################

"Get all parameters"
function get_parameters(@nospecialize(reaction::AbstractReaction))
    if hasfield(typeof(reaction), :pars)
        all_parameters = AbstractParameter[v for v in getfield(reaction, :pars)::NamedTuple]
    else
        all_parameters = AbstractParameter[]
    end

    return all_parameters
end

"Get parameter by name"
function get_parameter(reaction::AbstractReaction, parname::AbstractString; allow_not_found=false)

    if hasfield(typeof(reaction), :pars) && hasfield(typeof(reaction.pars), Symbol(parname))
        par = getfield(reaction.pars, Symbol(parname))
    else
        allow_not_found ||
            error("configuration error, Reaction $(fullname(reaction)) $(reaction)\n",
                "has no parameter name='$(parname)' (available parameters $([p.name for p in get_parameters(reaction)]))")
        par = nothing
    end

    return par
end

set_parameter_value!(@nospecialize(reaction::AbstractReaction), parname::AbstractString, value) = 
    setvalue!(get_parameter(reaction, parname, allow_not_found=false), value)

get_parameter_value(@nospecialize(reaction::AbstractReaction), parname::AbstractString) = 
    get_parameter(reaction, parname, allow_not_found=false).v


"Get method by name"
get_method_setup(@nospecialize(reaction::AbstractReaction), methodname::AbstractString; allow_not_found=false) =
    _get_method(reaction.base.methods_setup, "setup", reaction, methodname, allow_not_found)

get_method_initialize(@nospecialize(reaction::AbstractReaction), methodname::AbstractString; allow_not_found=false) =
    _get_method(reaction.base.methods_initialize, "initialize", reaction, methodname, allow_not_found)

get_method_do(@nospecialize(reaction::AbstractReaction), methodname::AbstractString; allow_not_found=false) =
    _get_method(reaction.base.methods_do, "do", reaction, methodname, allow_not_found)

function _get_method(methodlist, methodtype, reaction, methodname, allow_not_found)

    matchmethods = filter(m -> m.name==methodname, methodlist)
    
    length(matchmethods) <= 1 ||
        error("coding error: duplicate method $(methodtype) name '$(name)' for Reaction: ", reaction)
    
    !isempty(matchmethods) || allow_not_found ||
        error("configuration error in $(fullname(reaction)) parameters: $(reaction) has no method $(methodtype) name='$(methodname)'")

    return isempty(matchmethods) ? nothing : matchmethods[1]
end

"""
    get_variables(reaction, localname) -> Vector{VariableReaction}
    get_variables(reaction; [filterfn=v->true]) -> Vector{VariableReaction}

Get matching Variables from all ReactionMethods.
"""
function get_variables(
    @nospecialize(reaction::AbstractReaction), localname::AbstractString
)
    return get_variables(reaction, filterfn = v -> v.localname==localname)
end

function get_variables(
    @nospecialize(reaction::AbstractReaction); 
    filterfn = v -> true
)
    matchvars = VariableReaction[]
    for methodlist in (reaction.methods_setup, reaction.methods_initialize, reaction.methods_do)
        for m in methodlist
            append!(matchvars, get_variables(m;  filterfn=filterfn))
        end
    end

    return matchvars
end

"""
    get_variable(reaction, methodname, localname) -> VariableReaction or nothing

Get a single VariableReaction or nothing if match not found.
"""
function get_variable(
    @nospecialize(reaction::AbstractReaction), methodname::AbstractString, localname::AbstractString;
    allow_not_found=false
)
    matchvars = get_variables(reaction, filterfn = v -> (v.method.name == methodname && v.localname==localname))

    length(matchvars) <= 1 ||
        error("duplicate variable localname", localname)

    !isempty(matchvars) || allow_not_found || 
        error("method ", fullname(method), " no variable localname=", localname)

    return isempty(matchvars) ? nothing : matchvars[1]
end


##################################################
# Initialization callbacks
#################################################

"""
    set_model_geometry(reaction, model)

Optional: define [`Domain`](@ref) `grid`, `data_dims`.

One Reaction per Domain may create a grid (an [`AbstractMesh`](@ref) subtype) and set the `domain.grid` field.

Multiple Reactions per domain may call [`set_data_dimension!`](@ref) to define (different) named data dimensions (eg wavelength grids).
"""
function set_model_geometry(reaction::AbstractReaction, model::Model)
end

""" 
    check_configuration(reaction, model) -> Bool

Optional: check configuration is valid, log errors and return `false` if not
"""
function check_configuration(reaction::AbstractReaction, model::Model)
    return true
end

"""
    register_methods!(reaction)
    register_methods!(reaction, model)

Add [`ReactionMethod`](@ref)s, using [`add_method_setup!`](@ref), [`add_method_initialize!`](@ref)
[`add_method_do!`](@ref). See also [`register_dynamic_methods!`](@ref).
"""
register_methods!(reaction::AbstractReaction, model::Model) = register_methods!(reaction)
register_methods!(reaction::AbstractReaction) =
    error("$(typename(reaction)) register_methods! not implemented")


"""
    register_dynamic_methods!(reaction)
    register_dynamic_methods!(reaction, model)

Optional: called after first variable link pass, to
add [`ReactionMethod`](@ref)s that depend on Variables generated by other Reactions
(see [`register_methods!`](@ref)).
"""
register_dynamic_methods!(reaction::AbstractReaction, model::Model) = 
    register_dynamic_methods!(reaction)
function register_dynamic_methods!(reaction::AbstractReaction)
end

########################################################################################
# Helper functions to allow a Reaction implementation to populate Parameter and ReactionMethod lists
#########################################################################################


"""
    add_method_setup!(reaction::AbstractReaction, method::AbstractReactionMethod)
    add_method_setup!(reaction::AbstractReaction, methodfn::Function, vars::Tuple{Vararg{AbstractVarList}}; kwargs...) -> ReactionMethod

Add or create-and-add a setup method (called before main loop) eg to set persistent data or initialize state variables.
`methodfn`, `vars`, `kwargs` are passed to [`ReactionMethod`](@ref).
"""
function add_method_setup!(@nospecialize(reaction::AbstractReaction), method::AbstractReactionMethod)
    push!(reaction.base.methods_setup, method)
    return nothing
end

add_method_setup!(@nospecialize(reaction::AbstractReaction), methodfn::Function, vars::Tuple{Vararg{AbstractVarList}}; kwargs...) = 
    _add_method!(reaction, methodfn, vars, add_method_setup!; kwargs...)

"""
    add_method_initialize!(reaction::AbstractReaction, method::AbstractReactionMethod)
    add_method_initialize!(reaction::AbstractReaction, methodfn::Function, vars::Tuple{Vararg{AbstractVarList}}; kwargs...) -> ReactionMethod

Add or create-and-add an initialize method (called at start of each main loop iteration) 
eg to zero out accumulator Variables.
`methodfn`, `vars`, `kwargs` are passed to [`ReactionMethod`](@ref).
"""
function add_method_initialize!(@nospecialize(reaction::AbstractReaction), method::AbstractReactionMethod)
    push!(reaction.base.methods_initialize, method)
    return nothing
end

add_method_initialize!(@nospecialize(reaction::AbstractReaction), methodfn::Function, vars::Tuple{Vararg{AbstractVarList}}; kwargs...) = 
    _add_method!(reaction, methodfn, vars, add_method_initialize!; kwargs...)


"""
    add_method_do!(reaction::AbstractReaction, method::AbstractReactionMethod)
    add_method_do!(reaction::AbstractReaction, methodfn::Function, vars::Tuple{Vararg{AbstractVarList}}; kwargs...) -> ReactionMethod

Add or create and add a main loop method.
`methodfn`, `vars`, `kwargs` are passed to [`ReactionMethod`](@ref).
"""
function add_method_do!(@nospecialize(reaction::AbstractReaction), @nospecialize(method::AbstractReactionMethod))
    push!(reaction.base.methods_do, method)
    return nothing
end

add_method_do!(@nospecialize(reaction::AbstractReaction), @nospecialize(methodfn::Function), @nospecialize(vars::Tuple{Vararg{AbstractVarList}}); kwargs...) = 
    _add_method!(reaction, methodfn, vars, add_method_do!; kwargs...)

default_preparefn(m, vardata) = vardata

function _add_method!(
    @nospecialize(reaction::AbstractReaction), @nospecialize(methodfn::Function), @nospecialize(vars::Tuple{Vararg{AbstractVarList}}), add_method_fn;
    name=string(methodfn),
    p=nothing,
    preparefn=default_preparefn,
    operatorID=reaction.operatorID,
    domain=reaction.domain
)

    method = ReactionMethod(
        methodfn, 
        reaction,
        name,
        vars,
        p,
        copy(operatorID),
        domain,
        preparefn=preparefn
    )
    add_method_fn(reaction, method)
    return method
end



###################################
# creation from _cfg.yaml
##################################


"Create new ReactionXXXX of specified classname from config file, set parameters,
read variable mapping and initial values (to be set later by _configure_variables!)"
function create_reaction_from_config(
    classname::AbstractString,
    rdict::Dict{String, Type},
    domain::Domain,
    name,
    conf_reaction,
    external_parameters::Dict{String, Any}
)
 
    io = IOBuffer()

    local newreaction, conf_parameters  # make available outside try block

    try
        println(io, "create_reaction_from_config: $(domain.name).$name classname $classname")

        newreaction = create_reaction(rdict, classname, name, external_parameters)
        newreaction.base.domain = domain

        for k in keys(conf_reaction)
            if !(k in ("operatorID", "parameters", "variable_links", "variable_attributes"))
                error("reaction $(fullname(newreaction)) configuration error invalid key '$k'")
            end
        end

        operatorID = get(conf_reaction, "operatorID", nothing)
        if !isnothing(operatorID)
            println(io, "    operatorID=$(operatorID)")
            newreaction.base.operatorID = operatorID
        end

        # set parameters
        allpars = get_parameters(newreaction)
        conf_parameters_raw = get(conf_reaction, "parameters", Dict{Any,Any}())  # empty 'parameters:' will return nothing
        conf_parameters = isnothing(conf_parameters_raw) ? Dict{Any, Any}() : copy(conf_parameters_raw)
        for par in allpars
            rawvalue = par.v
            par_modified = false
            if par.external && haskey(newreaction.base.external_parameters, par.name)
                rawvalue = newreaction.base.external_parameters[par.name]
                par_modified = true
            end
            if haskey(conf_parameters, par.name)
                rawvalue = pop!(conf_parameters, par.name)
                par_modified = true
            end
            par_modstr = ("[Default]     ", "[config.yaml] ")[1+par_modified]
            println(io, "    set parameters: $par_modstr $(rpad(par.name,20))=$(rawvalue)")
            value = substitutevalue(parentmodule(typeof(newreaction)), externalvalue(io, rawvalue, newreaction.base.external_parameters))
            if isa(rawvalue, AbstractString) && value != rawvalue
                par_modified = true
                println(io, "      after substitution $(par.name)=$(value)")
            end
            if par_modified
                setvalue!(par, value)
            end
        end
    finally
        @info String(take!(io))
    end

    if !isempty(conf_parameters)
        io = IOBuffer()
        println(io, "reaction $(fullname(newreaction)) has no Parameter(s):")
        for (k, v) in conf_parameters
            println(io, "    $k:    $v")
        end
        error(String(take!(io)))
    end

    # Read Variable configuration
    # These are applied later by _configure_variables! after Variables are created
    newreaction.base._conf_variable_links = get(conf_reaction, "variable_links", nothing)
    newreaction.base._conf_variable_attributes = get(conf_reaction, "variable_attributes", nothing)

    return newreaction
end


function _register_methods!(@nospecialize(reaction::AbstractReaction), model::Model)
    empty!(reaction.methods_setup)
    empty!(reaction.methods_initialize)
    empty!(reaction.methods_do)
    register_methods!(reaction, model)

    return nothing
end

"apply Variable configuration. 
NB: runs twice:
 - after register_methods! (with allow_missing true)
 - after register_dynamic_methods! (with allow_missing false)"
function _configure_variables(@nospecialize(reaction::AbstractReaction); allow_missing::Bool, dolog::Bool)
    function sortstarfirst(x, y) 
        if occursin("*", x) && !occursin("*", y)
            lt = true
        elseif !occursin("*", x) && occursin("*", y)
            lt = false
        else 
            lt = x < y
        end
        return lt
    end

    io = devnull

    if !isnothing(reaction.base._conf_variable_links) # missing or empty 'variable_links:' will return nothing
        # sort Dict so wild cards (ending in *) are processed first, so they can be selectively overridden
        cvl = sort!(OrderedCollections.OrderedDict(reaction.base._conf_variable_links), lt=sortstarfirst)
        if dolog
            io = IOBuffer()
            println(io, "_configure_variables: $(nameof(typeof(reaction))) $(fullname(reaction)) variable_links:")
        end
        for (name, fullmapnameraw) in cvl
            try
                fullmapname = externalvalue(io, fullmapnameraw, reaction.base.external_parameters)
                linkreq_domain, linkreq_subdomain, mapname = split_link_name(fullmapname)
                
                match_vars = _gen_var_names(get_variables(reaction), name, mapname)
                uniquelocalnames = Set{String}()

                if isempty(match_vars)
                    allow_missing || 
                        error("    set variable_links:  $name -> $fullmapname no variables match $name\n",
                            "      available Variable local names: ", unique([v.localname for v in get_variables(reaction)]))
                else
                    for (var, newname) = match_vars
                        linkreq_fullname = combine_link_name(linkreq_domain, linkreq_subdomain, newname)

                        # Variables may appear in multiple ReactionMethods, so just print a log message for the first one
                        dolog && !(var.localname in uniquelocalnames) && 
                            println(io, "    set variable_links: $(rpad(var.localname,20)) -> $linkreq_fullname")
                        push!(uniquelocalnames, var.localname)

                        var.linkreq_name = newname
                        var.linkreq_domain = linkreq_domain
                        var.linkreq_subdomain = linkreq_subdomain
                    end
                end
            catch
                io === devnull || @info String(take!(io))
                @warn "_configure_variables: error setting Variable link for $(nameof(typeof(reaction))) $(fullname(reaction)) $name"
                rethrow()
            end
        end
    end

    # set variable attributes
    if !isnothing(reaction.base._conf_variable_attributes)
        if dolog 
            io = (io === devnull) ? IOBuffer() : io
            println(io, "_configure_variables: $(nameof(typeof(reaction))) $(fullname(reaction)) variable_attributes:")
        end
        cva = reaction.base._conf_variable_attributes 
        for (nameattrib, rawvalue) in cva
            try
                split_na = split(nameattrib, (':', '%'))
                length(split_na) == 2 || error("   invalid variable:attribute or variable%attribute $nameattrib")
                name, attrib = split_na
                match_vars = _gen_var_names(get_variables(reaction), name, "not used") # no wild cards (trailing *) allowed
                uniquelocalnames = Set{String}()
                if isempty(match_vars)
                    allow_missing || error("    $nameattrib = $rawvalue no variables match $name")
                else
                    for (var, dummy) = match_vars
                        # if isnothing(var.linkvar)
                        #    @warn "    unable to set $(Symbol(attrib)) attribute, VariableReaction $(fullname(var)) not linked"
                        # else
                        #    set_attribute!(var.linkvar, Symbol(attrib), value)
                        # end

                        # Variables may appear in multiple ReactionMethods, so just print a log message for the first one
                        iofirstvar = (dolog && !(var.localname in uniquelocalnames)) ? io : devnull
                        push!(uniquelocalnames, var.localname)

                        println(iofirstvar, "    set attribute: $(rpad(var.localname,20)) :$(rpad(attrib,20)) = $(rpad(rawvalue, 20)) ")

                        value = externalvalue(iofirstvar, rawvalue, reaction.base.external_parameters)
                        set_attribute!(var, Symbol(attrib), value)                   
                    end
                end
            catch
                io === devnull || @info String(take!(io))
                @warn "_configure_variables: error setting Variable attribute for $(nameof(typeof(reaction))) $(fullname(reaction)) $nameattrib"
                rethrow()
            end
        end
    end
    
    io === devnull ||  @info String(take!(io))

    return nothing
end

"""
    _gen_var_names(variables, matchroot::AbstractString, newroot::AbstractString) -> Vector{Pair{VariableReaction, String}}

Return a list of (variable => newname)  for `variables` where localname matches `matchroot`,
generating `newname` from `newroot`.

If `matchroot` and `newroot` contain a *, treat this as a wildcard.
If `matchroot` contains a * and `newroot` doesn't (legacy form), append a * to `newroot` first.
"""
function _gen_var_names(variables, matchroot::AbstractString, newroot::AbstractString)
    #                       variable => newname
    match_vars = Vector{Pair{VariableReaction, String}}()   

    count("*", matchroot) <= 1 ||
        error("matchroot $matchroot contains more than one *")

    if contains(newroot, "*")
        newrootstar = newroot
    else
        # legacy form without a * on the RHS - interpret as if had a trailing *
        newrootstar = newroot*"*"
    end

    for var in variables
        if contains(matchroot, "*") && (length(var.localname) >= length(matchroot) - 1)
            starpos = findfirst('*', matchroot)
            beforestar = matchroot[1:starpos-1]
            afterstar = matchroot[starpos+1:end]
            if startswith(var.localname, beforestar) && endswith(var.localname, afterstar)
                starmatched = var.localname[length(beforestar)+1:end-length(afterstar)]  # string that * matched
                push!(match_vars, var => replace(newrootstar, "*"=>starmatched))
            end
        else
            if var.localname == matchroot
                push!(match_vars, var => newroot)
            end
        end
    end

    return match_vars
end



###########################################
# Pretty printing
############################################

function Base.show(io::IO, react::AbstractReaction)
    print(
        io, 
        typename(react),
        "(name='", react.name, 
        "', classname='", react.classname, 
        "', domain='", domainname(react) ,
        ", operatorID=", react.operatorID,
        ")"
    )
end

function Base.show(io::IO, ::MIME"text/plain", react::AbstractReaction)
    dump_reaction(io, react)
end

function dump_reaction(io::IO, react::AbstractReaction; prefix="", show_parameters::Bool=true)
    println(io, prefix, typename(react))
    println(io, prefix, "  name='", react.name, "'")
    println(io, prefix, "  classname='", react.classname, "'")
    println(io, prefix, "  domain='", domainname(react), "'")
    println(io, prefix, "  operatorID=", react.operatorID)
    show_parameters && println(io, prefix, "  parameters=", get_parameters(react))
    println(io, prefix, "  methods_setup=", react.methods_setup)
    println(io, prefix, "  methods_initialize=", react.methods_initialize)
    println(io, prefix, "  methods_do=", react.methods_do)
end

"""
    show_parameters(react::AbstractReaction) -> DataFrame
    show_parameters(classname::AbstractString) -> DataFrame

list all parameters for a Reaction `react` instance or `classname`
"""
function show_parameters(react::AbstractReaction)
    pars = get_parameters(react)
    dfreact = DataFrames.DataFrame()
    dfreact.name = [p.name for p in pars]
    dfreact.v   = [p.v for p in pars]
    dfreact.type = [typeof(p) for p in pars]
    dfreact.units = [p.units for p in pars]
    dfreact.description = [p.description for p in pars]
            
    DataFrames.sort!(dfreact, [:name])
    return dfreact
end

show_parameters(classname::AbstractString) = 
    show_parameters(_create_reaction(classname, "test", Dict{String,Any}()))

"fully-qualified name"
fullname(react::AbstractReaction) = domainname(react)*"."*react.name
"safe Reaction Domain name"
domainname(react::AbstractReaction) = isnothing(react.domain) ? "<no domain>" : react.domain.name
"type name, excluding (verbose) template arguments"
typename(react::AbstractReaction) = join(Base.fullname(parentmodule(typeof(react))), ".")*"."*String(nameof(typeof(react)))
