import Infiltrator

"""
    Model

A biogeochemical model consisting of [`Domain`](@ref)s, created from a [YAML](https://en.wikipedia.org/wiki/YAML)
configuration file using [`create_model_from_config`](@ref).
"""
Base.@kwdef mutable struct Model
    name::String
    config_files::Vector{String}
    parameters::Dict{String, Any}                    
    domains::Vector{AbstractDomain}                 = Vector{AbstractDomain}()
    sorted_methods_setup                            = nothing
    sorted_methods_initialize                       = nothing
    sorted_methods_do                               = nothing
end



"Get number of model domains"
function get_num_domains(model::Model)
    return length(model.domains)
end

"""
    get_domain(model::Model, name::AbstractString) -> Domain or nothing
    get_domain(model::Model, domainid) -> Domain

Get Domain by `name` (may be nothing if `name` not matched)
or `domainid` (range 1:num_domains).
"""
function get_domain(model::Model, name::AbstractString)

    domainidx = findfirst(d -> d.name==name, model.domains)
    
    if isnothing(domainidx)
        return nothing
    else
        return model.domains[domainidx]
    end
end

function get_domain(model::Model, domainid)
    if domainid > length(model.domains)
        error("get_domain invalid domainid=", domainid)
    end

    return model.domains[domainid]    
end

function get_mesh(model::Model, domainname::AbstractString)
    dom = get_domain(model, domainname)
    !isnothing(dom) ||
        throw(ArgumentError("no Domain $domainname"))
    return dom.grid
end

"""
    get_reaction(model::Model, domainname, reactionname) -> Reaction or nothing

Get Reaction by domainname and reaction name
"""
function get_reaction(model::Model, domainname, reactionname)
    domain = get_domain(model, domainname)
    if !isnothing(domain)
        return get_reaction(domain, reactionname)
    end
    return nothing
end

"""
    get_variable(model::Model, varnamefull) -> VariableDomain or nothing

Get Variable by name of form `<domainname>.<variablename>`
"""
function get_variable(model::Model, varnamefull::AbstractString)

    varsplit = split(varnamefull, ".")
    length(varsplit) == 2 || 
        throw(ArgumentError("varnamefull $varnamefull is not of form <domainname>.<variablename>"))
    domainname, variablename = varsplit

    domain = get_domain(model, domainname)
    if !isnothing(domain)
        return get_variable(domain, variablename)
    end
    return nothing
end

"""
    get_field(model::Model, modeldata, varnamefull) -> Field

Get [`Field`](@ref) by Variable name
"""
function get_field(model::Model, modeldata::AbstractModelData, varnamefull::AbstractString)
    var = get_variable(model, varnamefull)
    !isnothing(var) ||
        throw(ArgumentError("Variable $varnamefull not found"))
    return get_field(var, modeldata)
end


"""
    set_variable_attribute!(model::Model, varnamefull, attributename::Symbol, value)

Set `varnamefull` (of form <domain name>.<var name>) `attributename` to `value`.
"""
function set_variable_attribute!(
    model::Model, 
    varnamefull::AbstractString,
    attributename::Symbol,
    value
)
    domvar = get_variable(model, varnamefull)
    !isnothing(domvar) || 
        throw(ArgumentError("Variable $varnamefull not found"))

    set_attribute!(domvar, attributename, value)  
    
    return nothing
end

set_variable_attribute!(
    model::Model, 
    domainname::AbstractString,
    variablename::AbstractString,
    attributename::Symbol,
    value
) = set_variable_attribute!(model, domainname*"."*variablename, attributename, value)

"convenience function to get VariableReaction with 'localname' from all ReactionMethods of 'reactionname'"
function get_reaction_variables(
    model::Model, domainname::AbstractString, reactionname::AbstractString, localname::AbstractString
)

    react = get_reaction(model, domainname, reactionname)
    !isnothing(react) || error("Reaction $domainname.$reactionname not found")

    return get_variables(react, localname)
end

"convenience function to get VariableReaction matching 'filterfn' from all ReactionMethods of 'reactionname'"
function get_reaction_variables(
    model::Model, domainname::AbstractString, reactionname::AbstractString;
    filterfn = v -> true,
)

    react = get_reaction(model, domainname, reactionname)
    !isnothing(react) || error("Reaction $domainname.$reactionname not found")

    return get_variables(react, filterfn=filterfn)
end


"""
    set_parameter_value!(model::Model, domainname, reactionname, parname, value)

Convenience function to set Parameter value.
"""
function set_parameter_value!(
    model::Model, 
    domainname::AbstractString,
    reactionname::AbstractString,
    parname::AbstractString,
    value
)
    react = get_reaction(model, domainname, reactionname)
    !isnothing(react) || error("Reaction $domainname.$reactionname not found")

    set_parameter_value!(react, parname, value)  
    
    return nothing
end

###################################
# creation from _cfg.yaml
##################################

"""
    create_model_from_config(
        config_file::AbstractString, configmodel::AbstractString;
        modelpars::Dict = Dict()
    ) -> model::Model
    
    create_model_from_config(
        config_files,
        configmodel::AbstractString;
        modelpars::Dict = Dict()
    ) -> model::Model

Construct model from a single [YAML](https://en.wikipedia.org/wiki/YAML) `config_file`, or from a collection of `config_files`, 
which are read in order and concatenated before being parsed as yaml.

Optional argument `modelpars` provides parameters that override those in `<configmodel>:` `parameters:` section.
"""
function create_model_from_config(
    config_file::AbstractString,
    configmodel::AbstractString; 
    kwargs...
)
    return create_model_from_config([config_file], configmodel; kwargs...)
end

function create_model_from_config(
    config_files, 
    configmodel::AbstractString;
    modelpars::Dict=Dict(),
    sort_methods_algorithm=group_methods,
)

    @info lpad("", 80, "=")
    @info "create_model_from_config file(s) $(abspath.(config_files)) configmodel $configmodel"
    @info lpad("", 80, "=")
    yamltext = ""
    for fn in config_files
        isfile(fn) || error("config file $(abspath(fn)) not found")
        yamltext *= read(fn, String)
    end

    data = nothing
    try
        data = YAML.load(yamltext)
    catch e
        error("$(typeof(e)) while reading .yaml config file(s) $(abspath.(config_files)).\n"*
            "If the error isn't obvious by looking at the file(s) (often this is a whitespace issue), "*
            "try an online YAML validator eg http://www.yamllint.com")
    end

    conf_model = data[configmodel]
    for k in keys(conf_model)
        if !(k in ("parameters", "domains", "geometry_precedence"))
            error("Model configuration error invalid key '$k'")
        end
    end
    
    parameters = get(conf_model, "parameters", Dict{String, Any}())  
    if isnothing(parameters)
        @warn "empty Model.parameters"
        parameters = Dict{String,Any}()
    end
    # extrapars can override parameter from configuration file, but not create new parameters (to catch typos)
    for (parname, value) in modelpars
        if haskey(parameters, parname)
            @info "    Resetting Model parameter $(rpad(parname, 20)) = $value (configuration file had value=$value)"
            parameters[parname] = value
        else
            error("configuration error: modelpars parameter $parname not present in configuration file for Model parameters:")
        end
    end
    for (parname, value) in parameters
        @info "    Model parameter: $(rpad(parname,20)) = $value"
    end

    @timeit "Model" model = Model(
        config_files=[abspath(fn) for fn in config_files],
        name=configmodel,
        parameters=parameters
    )

    conf_domains = conf_model["domains"]

    # create domains
    @timeit "create Domains" begin
    rdict = find_all_reactions()
    @info "generated Reaction catalog with $(length(rdict)) Reactions"
    for (name, conf_domain) in conf_domains
        nextDomainID = length(model.domains) + 1
        @info "creating domain '$(name)' ID=$nextDomainID" 
        if isnothing(conf_domain) # empty domain will return nothing
            conf_domain = Dict()
        end       
        push!(
            model.domains, 
            create_domain_from_config(name, nextDomainID, conf_domain, model.parameters, rdict)
        )     
    end
    end # timeit
    # request configuration of Domain sizes, Subdomains
    for dom in model.domains
        for react in dom.reactions
            set_model_geometry(react, model)
        end
    end

    # register methods
    @timeit "_register_reaction_methods" _register_reaction_methods!(model)

    # Link variables
    @timeit "_link_variables" _link_variables!(model)

    # sort methods
    @timeit "_sort_method_dispatch" _sort_method_dispatch!(model, sort_methods_algorithm=sort_methods_algorithm)

    return model
end


###################################
# setup
##################################

"""
    create_modeldata(model::Model [, eltype] [; threadsafe]) -> modeldata::ModelData

Create a new [`ModelData`](@ref) struct for model variables of element type `eltype`.
"""
function create_modeldata(
    model::Model, eltype::DataType=Float64;
    threadsafe=false,
    allocatenans=true, # fill Arrays with NaN when first allocated
)
    modeldata = ModelData{eltype}(model=model, threadsafe=threadsafe, allocatenans=allocatenans)

    for dom in model.domains
        push!(modeldata.domain_data, DomainData(dom))      
    end

    modeldata.cellranges_all = create_default_cellrange(model)

    return modeldata
end

"""
    create_default_cellrange(model::Model [; operatorID=0]) -> Vector{AbstractCellRange}

Create a `Vector` of `CellRange` instances covering the entire model.
"""
function create_default_cellrange(model::Model ; operatorID=0)
    cellranges_all = Vector{AbstractCellRange}()

    for dom in model.domains       
        push!(cellranges_all, Grids.create_default_cellrange(dom, dom.grid, operatorID=operatorID))
    end

    return cellranges_all
end

"""
    allocate_variables!(model, modeldata; [hostdep] [,eltypemap])
    
Allocate memory for Domain variables for every Domain in `model`.

See [`allocate_variables!(domain::Domain, modeldata::AbstractModelData)`](@ref).
"""
function allocate_variables!(
    model::Model, modeldata::AbstractModelData; 
    hostdep::Union{Bool,Nothing}=nothing,
    eltypemap=Dict{String, DataType}()
)
    @info lpad("", 80, "=")
    @info "allocate_variables!"
    @info lpad("", 80, "=")

    check_modeldata(model, modeldata)
    for dom in model.domains
        allocate_variables!(dom, modeldata, hostdep=hostdep, eltypemap=eltypemap)
    end
   
    return nothing
end



"""
    check_ready(model, modeldata; [throw_on_error=true] [, check_hostdep_varnames=true] [, expect_hostdep_varnames=["global.tforce"]]) -> ready

Check all variable pointers set, and no unexpected host-dependent non-state Variables (ie unlinked Variables)
"""
function check_ready(
    model::Model, modeldata::AbstractModelData;
    throw_on_error=true,
    check_hostdep_varnames=true,
    expect_hostdep_varnames=["global.tforce"],
)
    check_modeldata(model, modeldata)
    ready = true
    for dom in model.domains
        # don't exit on first error so display report from all Domains
        ready = ready && check_ready(dom, modeldata, throw_on_error=false)
    
        if check_hostdep_varnames
            dom_hdv, _ = get_host_variables(dom, VF_Undefined)
            for hv in dom_hdv
                fullname_hv = fullname(hv)
                if !(fullname_hv in expect_hostdep_varnames)
                    io = IOBuffer()
                    write(io, "check_ready: unexpected host-dependent Variable $fullname_hv (usually an unlinked Variable due to eg a "*
                    "missing renaming in the :variable_links sections in the .yaml file, a spelling mistake either "*
                    "in a Variable default name in the code or renaming in the .yaml file, or a missing Reaction)\n")
                    show_links(io, hv)
                    @warn String(take!(io))
                    ready = false
                end
            end
        end
    end    

    if !ready 
        @error "check_ready failed"
        throw_on_error && error("check_ready failed")
    end

    return ready
end

"Check configuration"
function check_configuration(
    model::Model;
    throw_on_error=true
)
    configok = true
    for dom in model.domains
        configok = configok && check_configuration(dom, model)
    end

    if !configok
        @error "check_configuration failed"
        throw_on_error && error("check_configuration failed")
    end

    return configok
end


"""
    initialize_reactiondata!(model::Model, modeldata::AbstractModelData)
    
Processes `VarList_...`s from ReactionMethods and populates `modeldata.sorted_methodsdata_...`
with sorted lists of ReactionMethods and corresponding Variable accessors.

Calls `create_dispatch_methodlists(model, modeldata, modeldata.cellranges_all)`
to set `modeldata.dispatchlists_all` to (optional) default  ReactionMethodDispatchLists for entire model.
"""
function initialize_reactiondata!(
    model::Model, modeldata::AbstractModelData;
    method_barrier=nothing,
    create_dispatchlists_all=true,
)
    @info lpad("", 80, "=")
    @info "initialize_reactiondata!"
    @info lpad("", 80, "=")

    check_modeldata(model, modeldata)

    # TODO using Ref here seems to trade off time to create ReactionMethodDispatchList
    # and time for first call to do_deriv ??
    # (Ref gives fast ReactionMethodDispatchList creation, but slow first do_deriv)
    # NB: passing Ref to call_method seems to speed up first do_deriv 

    modeldata.sorted_methodsdata_setup = 
        [
            Any[Ref(m), Ref(m.preparefn(m, create_accessors(m, modeldata)))]
            for m in get_methods(model.sorted_methods_setup)
        ]

    modeldata.sorted_methodsdata_initialize = 
        [
            Any[Ref(m), Ref(m.preparefn(m, create_accessors(m, modeldata)))]
            for m in get_methods(model.sorted_methods_initialize, method_barrier=method_barrier)
        ]

    modeldata.sorted_methodsdata_do = 
        [
            Any[Ref(m), Ref(m.preparefn(m, create_accessors(m, modeldata)))]
            for m in get_methods(model.sorted_methods_do, method_barrier=method_barrier)
        ]

    if create_dispatchlists_all
        # create default dispatchlists_all corresponding to cellranges_all
        modeldata.dispatchlists_all =
            create_dispatch_methodlists(model, modeldata, modeldata.cellranges_all)
    end

    return nothing
end


"""
    dispatch_setup(model, attribute_name, modeldata, cellranges=modeldata.cellranges_all)

Call setup methods, eg to initialize data arrays (including state variables).

`attribute_name` defines the setup operation performed. `dispatch_setup` should be called in sequence with `attribute_name` = :
- `:setup`: initialise Reactions and set up any non-state Variables (eg model grid Variables)
- `:norm_value`: set state Variable values from `:norm_value` attribute in .yaml file, and initialise any Reaction state that requires this value.
- `:initial_value` (optional): set state Variable values from `:initial_value` attribute in .yaml file.
"""
function dispatch_setup(
    model::Model, attribute_name, modeldata::AbstractModelData, cellranges=modeldata.cellranges_all
)
    @info lpad("", 80, "=")
    @info "dispatch_setup :$attribute_name"
    @info lpad("", 80, "=")

    check_modeldata(model, modeldata) 

    for (method, vardata) in modeldata.sorted_methodsdata_setup
        for cr in _dispatch_cellranges(method[], cellranges)
           call_method(method[], vardata[], cr, attribute_name)
        end
    end

    return nothing  
end

"""
    create_dispatch_methodlists(model::Model, modeldata::AbstractModelData, cellranges) 
        -> (;list_initialize::ReactionMethodDispatchList, list_do::ReactionMethodDispatchList)
 
Compile lists of `initialize` and `do` methods + corresponding `cellrange` for main loop [`do_deriv`](@ref).

Subset of methods and cellrange to operate on are generated from supplied `cellranges`.
"""
function create_dispatch_methodlists(
    model::Model, modeldata::AbstractModelData, cellranges;
    verbose=false,
    generated_dispatch=true,
)
    check_modeldata(model, modeldata)

    verbose && @info "list_initialize:\n"
    @timeit "list_initialize" list_initialize = _create_dispatch_methodlist(
        modeldata.sorted_methodsdata_initialize,
        cellranges,
        generated_dispatch,
    )
    verbose && @info "list_do:\n"
    @timeit "list_do" list_do = _create_dispatch_methodlist(
        modeldata.sorted_methodsdata_do,
        cellranges,
        generated_dispatch,
    )
    
    return (;list_initialize, list_do)
end


function _create_dispatch_methodlist(methodsdata, cellranges, generated_dispatch::Bool)
    methods, vardatas, crs = [], [], []
   
    @timeit "create arrays" begin
    for (method, vardata) in methodsdata
        for cr in _dispatch_cellranges(method[], cellranges)
            push!(methods, method)
            push!(vardatas, vardata)
            push!(crs, cr)
        end
    end
    end # timeit

    if generated_dispatch
        @timeit "ReactionMethodDispatchList" rmdl = ReactionMethodDispatchList(methods, vardatas, crs)
    else
        @timeit "ReactionMethodDispatchListNoGen" rmdl = ReactionMethodDispatchListNoGen(methods, vardatas, crs)
    end
    return rmdl
end

function _dispatch_cellranges(@nospecialize(method::AbstractReactionMethod), cellranges)
    if is_do_barrier(method)
        return (nothing, )
    else
        domain = method.domain
        operatorID = method.operatorID
        return Iterators.filter(
            cr -> (cr.domain === domain) && 
                (cr.operatorID == 0 || cr.operatorID in operatorID),
            cellranges
        )
    end
end

####################################
# Main loop methods
#####################################

"""
    do_deriv(dispatchlists, deltat::Float64=0.0)

Wrapper function to calculate entire derivative (initialize and do methods) in one call.
`dispatchlists` is from [`create_dispatch_methodlists`](@ref).
"""
function do_deriv(dispatchlists, deltat::Float64=0.0)

    dispatch_methodlist(dispatchlists.list_initialize)
    
    dispatch_methodlist(dispatchlists.list_do, deltat)

    return nothing
end


"""
    dispatch_methodlist(dl::ReactionMethodDispatchList, deltat::Float64=0.0)

Dispatch to a list of methods.

# Implementation

As an optimisation, uses @generated for Type stability
and to avoid dynamic dispatch, instead of iterating over lists.

[`ReactionMethodDispatchList`](@ref) fields are Tuples hence are fully Typed, the @generated
function emits unrolled code with a function call for each Tuple element. 
"""
@generated function dispatch_methodlist(
    dl::ReactionMethodDispatchList{M, V, C}, 
    deltat::Float64=0.0
) where {M, V, C}

    # See https://discourse.julialang.org/t/manually-unroll-operations-with-objects-of-tuple/11604
     
    ex = quote ; end  # empty expression
    for j=1:fieldcount(M)
        push!(ex.args,
            quote
                # let
                # call_method(dl.methods[$j][], dl.vardatas[$j][], dl.cellranges[$j], deltat)
                # pass Ref to function to reduce compile time
                call_method(dl.methods[$j], dl.vardatas[$j], dl.cellranges[$j], deltat)
                # end
            end
            )
    end
    push!(ex.args, quote; return nothing; end)
    
    return ex
end

function dispatch_methodlist(
    dl::ReactionMethodDispatchListNoGen, 
    deltat::Float64=0.0
)

   for (m, v, c) in zip(dl.methods, dl.vardatas, dl.cellranges)
        call_method(m[], v[], c, deltat)
   end

   return nothing
end

#################################
# Pretty printing
################################

"compact form"
function Base.show(io::IO, model::Model)
    print(io, "Model(config_files='", model.config_files,"', name='", model.name,"')")
end
"multiline form"
function Base.show(io::IO, ::MIME"text/plain", model::Model)
    println(io, "Model")
    println(io, "\tname='", model.name,"'")
    println(io, "\tconfig_files='", model.config_files,"'")  
    println(io, "\tdomains=", model.domains)
end

"""
    show_methods_setup(model::Model)

Display ordered list of Reaction setup methods (registered by [`add_method_setup!`](@ref),
called by [`dispatch_setup`](@ref))
"""
function show_methods_setup(model::Model)
    println("All methods_setup:")
    println(model.sorted_methods_setup)
    return nothing    
end

"""
    show_methods_initialize(model::Model)

Display ordered list of Reaction initialize methods (registered by [`add_method_initialize!`](@ref),
called by [`do_deriv`](@ref) at start of each model timestep).
"""
function show_methods_initialize(model::Model)
    println("All methods_initialize:")
    println(model.sorted_methods_initialize)
    return nothing    
end

"""
    show_methods_do(model::Model)

Display ordered list of Reaction do methods (registered by [`add_method_do!`](@ref),
called by [`do_deriv`](@ref) for each model timestep).
"""
function show_methods_do(model::Model)
    println("All methods_do:")
    println(model.sorted_methods_do)
    return nothing    
end

"""
    show_variables(model::Model; [attributes], [filter], showlinks=false, modeldata=nothing) -> DataFrame
    show_variables(model::Model, domainname; [attributes], [filter], showlinks=false, modeldata=nothing) -> DataFrame

Show table of Domain Variables. Optionally get variable links, data.

# Keywords:
See [`show_variables(domain::Domain, kwargs...)`](@ref) 

# Examples:
Display all model Variables using VS Code table viewer:

    julia> vscodedisplay(PB.show_variables(run.model))

Write out all model Variables as csv for import to a spreadsheet:

    julia> CSV.write("vars.csv", PB.show_variables(run.model, modeldata=modeldata, showlinks=true))
"""
function show_variables(model::Model; kwargs...)
    df = DataFrames.DataFrame()
    for dom in model.domains
        dfdom = show_variables(dom; kwargs...)
        # prepend domain name
        DataFrames.insertcols!(dfdom, 1, :domain=>fill(dom.name, size(dfdom,1)))
        # append to df
        df = vcat(df, dfdom)
    end
    DataFrames.sort!(df, [:domain, :name])
    return df
end

function show_variables(model::Model, domainname::AbstractString; kwargs...)
    dom = get_domain(model, domainname)
    !isnothing(dom) ||
        throw(ArgumentError("no Domain $domainname"))
    return show_variables(dom; kwargs...)
end

"""
    show_parameters(model) -> DataFrame

Get parameters for all reactions in model.

# Examples:
Show all model parameters using VS Code table viewer:

    julia> vscodedisplay(PB.show_parameters(run.model))

Write out all model parameters as csv for import to a spreadsheet:

    julia> CSV.write("pars.csv", PB.show_parameters(run.model))
"""
function show_parameters(model::Model)
    df = DataFrames.DataFrame()
    for dom in model.domains
        for react in dom.reactions
            dfreact = show_parameters(react)
            DataFrames.insertcols!(dfreact, 1, :reaction=>fill(react.name, size(dfreact,1)))
            DataFrames.insertcols!(dfreact, 2, :classname=>fill(react.classname, size(dfreact,1)))
            DataFrames.insertcols!(dfreact, 1, :domain=>fill(dom.name, size(dfreact,1)))
            
            # append to df
            df = vcat(df, dfreact)
        end
    end
    DataFrames.sort!(df, [:domain, :reaction, :name])
    return df
end
    
######################################
# Variable linking
#######################################

"""
    _link_variables!(model::Model) -> nothing

Populate `Domain` variable lists, renaming and linking Reaction variables
"""
function _link_variables!(model::Model)

    # First pass - variables defined in config file    
    foreach(_link_clear!, model.domains)

    @info "link_variables first pass"
    # _link_variables!(model, _link_print)
    _link_variables!(model, _link_create, false)
    _link_variables!(model, _link_create_contrib, false)
    _link_variables!(model, _link_create_dep, false)
    _link_variables!(model, _link_link, false)

    # Allow reactions to define additional variables, based on 
    # the model structure
    @info "link_variables initialize dynamic variables"
    _register_reaction_dynamic_methods!(model)

    # Second pass including additional variables
    foreach(_link_clear!, model.domains)

    @info "link_variables second pass:"
    _link_variables!(model, _link_print, true)
    _link_variables!(model, _link_create, true)
    _link_variables!(model, _link_create_contrib, true)
    _link_variables!(model, _link_create_dep, true)
    _link_variables!(model, _link_link, true)
    @info "link_variables unlinked variables:"
    _link_variables!(model, _link_print_not_linked, true)
    # 
   
    return nothing
end

"""
    relink_variables!(model::Model) -> nothing

Optional relink `Reaction` variables and repopulate `Domain` Variable lists.
Only needed if modifying `variable_links` after calling [`create_model_from_config`](@ref)

NB: dynamic Variables are not recreated, so will not reflect any changes
"""
function relink_variables!(model::Model)
    foreach(_link_clear!, model.domains)

    @info "relink_variables!:"
    _link_variables!(model, _link_print, true)
    _link_variables!(model, _link_create, true)
    _link_variables!(model, _link_create_contrib, true)
    _link_variables!(model, _link_create_dep, true)
    _link_variables!(model, _link_link, true)
    @info "relink_variables! unlinked variables:"
    _link_variables!(model, _link_print_not_linked, true)
    # 
   
    return nothing
end

function _link_variables!(model::Model, oper, dolog)
    for dom in model.domains
        _link_variables!(dom, model, oper, dolog)
    end
end

############################
# Method registration and ordering
##############################

function _register_reaction_methods!(model::Model)
    
    for dom in model.domains
        for r in dom.reactions
            _register_methods!(r, model)
            # not all methods registered and ReactionVariables available, so no error if missing
            _configure_variables(r; allow_missing=true, dolog=false)
        end
    end

    return nothing
end

function _register_reaction_dynamic_methods!(model::Model)
    
    for dom in model.domains
        for r in dom.reactions            
            register_dynamic_methods!(r, model)
            _configure_variables(r; allow_missing=false, dolog=true)
        end
    end

    return nothing
end

get_methods_setup(model; kwargs...) = _get_methods(model, :methods_setup)
get_methods_initialize(model; kwargs...) = _get_methods(model, :methods_initialize)
get_methods_do(model; kwargs...) = _get_methods(model, :methods_do)

function _get_methods(model::Model, mfield::Symbol; filterfn=m->!is_do_nothing(m))
    methods=ReactionMethod[]
    for dom in model.domains
        for r in dom.reactions
            append!(methods, filter(filterfn, getproperty(r, mfield)))
        end
    end
    return methods
end

"""
    _sort_method_dispatch!(model::Model) -> nothing

Sort methods into dispatch order, based on variable dependencies
"""
function _sort_method_dispatch!(model::Model; sort_methods_algorithm=group_methods)
    # create unsorted list of all Reactions in model
    all_reacts = Vector{AbstractReaction}()
    for dom in model.domains
        append!(all_reacts, dom.reactions)
    end
      
    # get all domain Variables
    all_domvars = Vector{VariableDomain}()
    for dom in model.domains
        append!(all_domvars, get_variables(dom))
    end

    methods_setup = get_methods_setup(model)
    model.sorted_methods_setup = sort_methods_algorithm(methods_setup, all_domvars)
    
    methods_initialize = get_methods_initialize(model)
    # TODO no sort needed as no dependencies allowed
    model.sorted_methods_initialize = sort_methods_algorithm(methods_initialize, all_domvars)

    methods_do = get_methods_do(model)
    model.sorted_methods_do = sort_methods_algorithm(methods_do, all_domvars)

    return nothing
end

