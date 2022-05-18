

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
    "Reaction parameters"
    parameters::Vector{AbstractParameter}           = Vector{AbstractParameter}()

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



"""
    Base.getproperty(react::AbstractReaction, s::Symbol)

Forward to `react.base::ReactionBase` to define additional properties.
"""
function Base.getproperty(react::AbstractReaction, s::Symbol)
    if s == :name
        return react.base.name
    elseif s == :classname
        return react.base.classname
    elseif s == :domain
        return react.base.domain
    elseif s == :operatorID
        return react.base.operatorID
    elseif s == :external_parameters
        return react.base.external_parameters
    elseif s == :methods_setup
        return react.base.methods_setup
    elseif s == :methods_initialize
        return react.base.methods_initialize
    elseif s == :methods_do
        return react.base.methods_do
    else
        return getfield(react, s)
    end
end

##########################################################
# Query methods for Parameters, ReactionMethods, Variables
##########################################################

"Get all parameters"
get_parameters(reaction::AbstractReaction) = reaction.base.parameters

"Get parameter by name"
function get_parameter(reaction::AbstractReaction, parname::AbstractString; allow_not_found=false)
    matchpars = filter(p -> p.name==parname, reaction.base.parameters)
    
    length(matchpars) <= 1 ||
        error("coding error: duplicate parameter name '$(name)' for Reaction: ", reaction)
    
    !isempty(matchpars) || allow_not_found ||
        error("configuration error, Reaction $(fullname(reaction)) $(reaction)\n",
            "has no parameter name='$(parname)' (available parameters $([p.name for p in get_parameters(reaction)]))")

    return isempty(matchpars) ? nothing : matchpars[1]
end

"convenience function to set Parameter value"
function set_parameter_value!(reaction::AbstractReaction, parname::AbstractString, value)

    setvalue!(get_parameter(reaction, parname, allow_not_found=false), value)  
    
    return nothing
end

"Get method by name"
get_method_setup(reaction::AbstractReaction, methodname::AbstractString; allow_not_found=false) =
    _get_method(reaction.base.methods_setup, "setup", reaction, methodname, allow_not_found)

get_method_initialize(reaction::AbstractReaction, methodname::AbstractString; allow_not_found=false) =
    _get_method(reaction.base.methods_initialize, "initialize", reaction, methodname, allow_not_found)

get_method_do(reaction::AbstractReaction, methodname::AbstractString; allow_not_found=false) =
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
    reaction::AbstractReaction, localname::AbstractString
)
    return get_variables(reaction, filterfn = v -> v.localname==localname)
end

function get_variables(
    reaction::AbstractReaction; 
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
    reaction::AbstractReaction, methodname::AbstractString, localname::AbstractString;
    allow_not_found=false
)
    matchvars = get_variables(reaction, filterfn = v -> (v.method.name == methodname && v.localname==localname))

    length(matchvars) <= 1 ||
        error("duplicate variable localname", localname)

    !isempty(matchvars) || allow_not_found || 
        error("method ", fullname(method), " no variable localname=", localname)

    return isempty(matchvars) ? nothing : matchvars[1]
end

"""
    get_rate_stoichiometry(rj::PB.AbstractReaction) -> Vector[(ratevarname, process, Dict(speciesname=>stoich, ...)), ...]

Optionally provide rate Variable name(s), a process name ("photolysis", "reaction", ...) and stoichiometry of reactants and products, for postprocessing of results.
"""
function get_rate_stoichiometry(rj::AbstractReaction)
    rs = []
    for m in rj.methods_do
        append!(rs, get_rate_stoichiometry(m))
    end
    return rs
end

##################################################
# Initialization callbacks
#################################################

"""
    set_model_geometry(reaction, model)

Optional: define model `Domain` size, Subdomains etc (only implemented by
Reactions that define grids).
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
    add_par(reaction::AbstractReaction, par::AbstractParameter)
    add_par(reaction::AbstractReaction, objectwithpars)

Add a single parameter or parameters from fields of `objectwithpars` to a new Reaction.

Not usually needed: Parameters in `pars::ParametersTuple`` will be added automatically, only needed if there are additional
Parameters that are not members of `pars`.
"""
function add_par(reaction::AbstractReaction, par::AbstractParameter)
    
    if isnothing(get_parameter(reaction, par.name, allow_not_found=true))
        push!(reaction.base.parameters, par)
    else
        error("attempt to add duplicate parameter name=''", par.name, "'' to reaction", reaction)
    end

    return nothing
end

function add_par(reaction::AbstractReaction, objectwithpars)
    for f in fieldnames(typeof(objectwithpars))
        if getfield(objectwithpars, f) isa AbstractParameter
            add_par(reaction, getfield(objectwithpars, f))
        end
    end
end


"""
    add_method_setup!(reaction::AbstractReaction, method::AbstractReactionMethod)
    add_method_setup!(reaction::AbstractReaction, methodfn::Function, vars::Tuple{Vararg{AbstractVarList}}; kwargs...) -> ReactionMethod

Add or create-and-add a setup method (called before main loop) eg to set persistent data or initialize state variables.
`methodfn`, `vars`, `kwargs` are passed to [`ReactionMethod`](@ref).
"""
function add_method_setup!(reaction::AbstractReaction, method::AbstractReactionMethod)
    push!(reaction.base.methods_setup, method)
    return nothing
end

add_method_setup!(reaction::AbstractReaction, methodfn::Function, vars::Tuple{Vararg{AbstractVarList}}; kwargs...) = 
    _add_method!(reaction, methodfn, vars, add_method_setup!; kwargs...)

"""
    add_method_initialize!(reaction::AbstractReaction, method::AbstractReactionMethod)
    add_method_initialize!(reaction::AbstractReaction, methodfn::Function, vars::Tuple{Vararg{AbstractVarList}}; kwargs...) -> ReactionMethod

Add or create-and-add an initialize method (called at start of each main loop iteration) 
eg to zero out accumulator Variables.
`methodfn`, `vars`, `kwargs` are passed to [`ReactionMethod`](@ref).
"""
function add_method_initialize!(reaction::AbstractReaction, method::AbstractReactionMethod)
    push!(reaction.base.methods_initialize, method)
    return nothing
end

add_method_initialize!(reaction::AbstractReaction, methodfn::Function, vars::Tuple{Vararg{AbstractVarList}}; kwargs...) = 
    _add_method!(reaction, methodfn, vars, add_method_initialize!; kwargs...)


"""
    add_method_do!(reaction::AbstractReaction, method::AbstractReactionMethod)
    add_method_do!(reaction::AbstractReaction, methodfn::Function, vars::Tuple{Vararg{AbstractVarList}}; kwargs...) -> ReactionMethod

Add or create and add a main loop method.
`methodfn`, `vars`, `kwargs` are passed to [`ReactionMethod`](@ref).
"""
function add_method_do!(reaction::AbstractReaction, method::AbstractReactionMethod)
    push!(reaction.base.methods_do, method)
    return nothing
end

add_method_do!(reaction::AbstractReaction, methodfn::Function, vars::Tuple{Vararg{AbstractVarList}}; kwargs...) = 
    _add_method!(reaction, methodfn, vars, add_method_do!; kwargs...)

function _add_method!(
    reaction::AbstractReaction, methodfn::Function, vars::Tuple{Vararg{AbstractVarList}}, add_method_fn;
    name=string(methodfn),
    p=nothing,
    preparefn=nothing,
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
       
    newreaction = create_reaction(rdict, classname, name, external_parameters)
    newreaction.base.domain = domain

    for k in keys(conf_reaction)
        if !(k in ("operatorID", "parameters", "variable_links", "variable_attributes"))
            error("reaction $(fullname(newreaction)) configuration error invalid key '$k'")
        end
    end

    operatorID = get(conf_reaction, "operatorID", nothing)
    if !isnothing(operatorID)
        @info "reaction $(fullname(newreaction)) operatorID=$(operatorID)"
        newreaction.base.operatorID = operatorID
    end

    # set parameters
    allpars = get_parameters(newreaction)
    conf_parameters = copy(get(conf_reaction, "parameters", Dict{Any,Any}()))
    for par in allpars
        rawvalue = par.v
        par_modified = false
        if par.external && haskey(newreaction.base.external_parameters, par.name)
            rawvalue = newreaction.base.external_parameters[par.name]
            par_modified = true
        end
        if !isnothing(conf_parameters) && haskey(conf_parameters, par.name) # empty 'parameters:' will return nothing
            rawvalue = pop!(conf_parameters, par.name)
            par_modified = true
        end
        par_modstr = ("[Default]     ", "[config.yaml] ")[1+par_modified]
        @info "    set parameters: $par_modstr $(rpad(par.name,20))=$(rawvalue)" 
        value = substitutevalue(parentmodule(typeof(newreaction)), externalvalue(rawvalue, newreaction.base.external_parameters))
        if isa(rawvalue, AbstractString) && value != rawvalue
            par_modified = true
            @info "      after substitution $(par.name)=$(value)"
        end
        if par_modified
            setvalue!(par, value)
        end
    end    

    if !isempty(conf_parameters)
        io = IOBuffer()
        write(io, "reaction $(fullname(newreaction)) has no Parameter(s):\n")
        for (k, v) in conf_parameters
            write(io, "    $k:    $v\n")
        end
        error(String(take!(io)))
    end

    # Read Variable configuration
    # These are applied later by _configure_variables! after Variables are created
    newreaction.base._conf_variable_links = get(conf_reaction, "variable_links", nothing)
    newreaction.base._conf_variable_attributes = get(conf_reaction, "variable_attributes", nothing)

    return newreaction
end


function _register_methods!(reaction::AbstractReaction, model::Model)
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
function _configure_variables(reaction::AbstractReaction; allow_missing::Bool, dolog::Bool)
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

    if !isnothing(reaction.base._conf_variable_links) # missing or empty 'variable_links:' will return nothing
        # sort Dict so wild cards (ending in *) are processed first, so they can be selectively overriden
        cvl = sort(reaction.base._conf_variable_links, lt=sortstarfirst)
        dolog && @info "_configure_variables: $(nameof(typeof(reaction))) $(fullname(reaction))"
        for (name, fullmapnameraw) in cvl
            fullmapname = externalvalue(fullmapnameraw, reaction.base.external_parameters)
            linkreq_domain, linkreq_subdomain, mapname = split_link_name(fullmapname)
            
            match_vars = _gen_var_names(get_variables(reaction), name, mapname)
            uniquelocalnames = Set{String}()

            if isempty(match_vars)
                allow_missing || 
                    error("    $(nameof(typeof(reaction))) $(fullname(reaction)) set variable_links:  $name -> $fullmapname no variables match $name\n",
                          "      available Variable local names: ", unique([v.localname for v in get_variables(reaction)]))
            else
                for (var, newname) = match_vars
                    linkreq_fullname = combine_link_name(linkreq_domain, linkreq_subdomain, newname)

                    # Variables may appear in multiple ReactionMethods, so just print a log message for the first one
                    dolog && !(var.localname in uniquelocalnames) && 
                        @info "    set variable_links: $(rpad(var.localname,20)) -> $linkreq_fullname"
                    push!(uniquelocalnames, var.localname)

                    var.linkreq_name = newname
                    var.linkreq_domain = linkreq_domain
                    var.linkreq_subdomain = linkreq_subdomain
                end
            end
        end
    end

    # set variable attributes
    if !isnothing(reaction.base._conf_variable_attributes)
        cva = reaction.base._conf_variable_attributes 
        for (nameattrib, value) in cva
            split_na = split(nameattrib, (':', '%'))
            length(split_na) == 2 || error("   $(nameof(typeof(reaction))) $(fullname(reaction)) set variable_attributes: invalid attribute $nameattrib")
            name, attrib = split_na
            match_vars = _gen_var_names(get_variables(reaction), name, "not used") # no wild cards (trailing *) allowed
            uniquelocalnames = Set{String}()
            if isempty(match_vars)
                allow_missing || error("    $(nameof(typeof(reaction))) $(fullname(reaction)) $nameattrib = $value no variables match $name")
            else
                for (var, dummy) = match_vars
                    # if isnothing(var.linkvar)
                    #    @warn "    unable to set $(Symbol(attrib)) attribute, VariableReaction $(fullname(var)) not linked"
                    # else
                    #    set_attribute!(var.linkvar, Symbol(attrib), value)
                    # end

                    # Variables may appear in multiple ReactionMethods, so just print a log message for the first one
                    dolog && !(var.localname in uniquelocalnames) && 
                        @info "    set attribute: $(rpad(var.localname,20)) :$(rpad(attrib,20)) = $(rpad(value, 20)) "
                    push!(uniquelocalnames, var.localname)

                    set_attribute!(var, Symbol(attrib), value)                   
                end
            end
        end
    end
    
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

"compact form"
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

"multiline form"
function Base.show(io::IO, ::MIME"text/plain", react::AbstractReaction)
    println(io, typename(react))
    println(io, "  name='", react.name, "'")
    println(io, "  classname='", react.classname, "'")
    println(io, "  domain='", domainname(react), "'")
    println(io, "  operatorID=", react.operatorID)
    println(io, "  parameters=", get_parameters(react))
    println(io, "  methods_setup=", react.methods_setup)
    println(io, "  methods_initialize=", react.methods_initialize)
    println(io, "  methods_do=", react.methods_do)
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
