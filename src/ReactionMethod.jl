




"""
    ReactionMethod

Defines a callback function `methodfn` with Variables `vars`, to be called from the Model framework
either during setup or as part of the main loop.

# Fields
- `methodfn::Function`: callback from Model framework
- `reaction::AbstractReaction`: the Reaction that created this ReactionMethod.
- `name::String`: a descriptive name, eg generated from the name of methodfn
- `vars::Tuple{Vararg{AbstractVarList}}`: [`AbstractVarList`](@ref)s of [`VariableReaction`](@ref)s.
  Corresponding Variable accessors (views on Arrays) will be provided to the `methodfn` callback.
- `p::P`: optional context field (of arbitrary type) to store data needed by methodfn.
- `operatorID::Vector{Int64}`: operator ID (to allow operator splitting)
- `domain::Domain`: the Domain that `methodfn` operates on (ie will be called with `cellrange` for this `domain`).
- `preparefn::Function`: optional function to eg define buffers before main loop starts.

# methodfn
The `methodfn` callback is:

    methodfn(m::ReactionMethod, vardata::Tuple, cellrange::AbstractCellRange, modelctxt)

With arguments:
- `m::ReactionMethod`: context is available as `m.reaction::AbstractReaction` (the Reaction that defined the `ReactionMethod`),
  and `m.p` (an arbitrary extra context field supplied when `ReactionMethod` created).
- `vardata`: A Tuple of collections of views on Domain data arrays, defined by `vars`
- `cellrange::AbstractCellRange`: range of cells to calculate.
- `modelctxt`:
    - for a setup method, `:initial_value` or `:norm_value` defining the type of setup requested
    - for a main loop method `deltat` providing timestep information eg for rate throttling.

# preparefn
An optional `preparefn` callback can be supplied eg to allocate buffers that require knowledge of
the data types of `vardata` or to cache expensive calculations:

    preparefn(m::ReactionMethod, vardata::Tuple) -> vardata

This is called after model arrays are allocated, and prior to setup.
"""
mutable struct ReactionMethod{M, R, P, V} <: AbstractReactionMethod
    methodfn::M 
    reaction::R
    name::String
    vars::V
    p::P
    operatorID::Vector{Int64}
    domain::Domain
    preparefn::Union{Nothing, Function} # NB: not concretely typed as not performance-critical

    function ReactionMethod(
        methodfn::M,
        reaction::R,
        name,
        vars::V,
        p::P, 
        operatorID::Vector{Int64}, 
        domain; 
        preparefn = nothing
    ) where {M <: Function, R <: AbstractReaction, P, V <: Tuple{Vararg{AbstractVarList}}}        
        
        copy_vars = deepcopy(vars)

        newmethod = new{M, R, P, V}(
            methodfn, 
            reaction,
            name,
            copy_vars,
            p,
            operatorID, 
            domain, 
            preparefn
        )

        for va in copy_vars
            for v in get_variables(va)
                v.method = newmethod
                set_attribute!(v, :operatorID, newmethod.operatorID, allow_create=true)
            end
        end

        return newmethod
    end
end

"Get VariableReactions from `vars` as a Tuple"
function get_variables_tuple(method::AbstractReactionMethod)
    vars_tuple = []
    for va in method.vars
        push!(vars_tuple, get_variables(va))
    end
    return Tuple(vars_tuple)
end

"""
    get_variables(method::AbstractReactionMethod; filterfn = v -> true) -> Vector

Get VariableReactions from `method.vars` as a flat Vector, optionally restricting to those that match `filterfn`
"""
function get_variables(
    method::AbstractReactionMethod;
    filterfn = v -> true
)
    vars = VariableReaction[]
    for va in method.vars
        append!(vars, filter(filterfn, get_variables(va)))
    end
    return vars
end

"Get a single VariableReaction by localname"
function get_variable(
    method::AbstractReactionMethod, localname::AbstractString; 
    allow_not_found=false
)
    matchvars = get_variables(method; filterfn = v -> v.localname==localname)
    length(matchvars) <= 1 || error("duplicate variable localname", localname)

    !isempty(matchvars) || allow_not_found || error("method ", fullname(method), " no variable localname=", localname)

    return isempty(matchvars) ? nothing : matchvars[1]
end

fullname(method::AbstractReactionMethod) = fullname(method.reaction)*"."*method.name

is_method_setup(method::AbstractReactionMethod) = (method in method.reaction.methods_setup)
is_method_initialize(method::AbstractReactionMethod) = (method in method.reaction.methods_initialize)
is_method_do(method::AbstractReactionMethod) = (method in method.reaction.methods_do)

get_rate_stoichiometry(m::ReactionMethod) = []

###########################################
# Pretty printing
############################################

"compact form"
function Base.show(io::IO, method::ReactionMethod)
    print(
        io, 
        "ReactionMethod(fullname='", fullname(method), 
            "', methodfn=", string(method.methodfn), 
            ", domain='", method.domain.name ,
            "', operatorID=", method.operatorID,
            ")",
    )
end




##############################
# Default methods (optional)
##############################

"""
    add_method_setup_initialvalue_vars_default!(react::AbstractReaction, variables [; kwargs...])

Create and add a default method to initialize Variables matching `filterfn` (defaults to state Variables)
to Attribute :initial_value at beginning of integration.

# Keywords
- `filterfn`: set to f(var)::Bool to override the default selection on :vfunction for state variables.
- `transfer_attribute_vars=[]`: Set to a list of the same length as `variables` to initialise `variables`
  from attributes of `transfer_attribute_vars`.
- `setup_callback=(method, attribute_value, var, vardata) -> nothing`: Set to a function that is called 
  after each Variable initialisation eg to store norm values.
- `convertvars=[]`
- `convertfn = (convertvars_tuple, i) -> 1.0`
- `convertinfo = ""`

# Including volume etc conversion
Set `convertvars` to a Vector of Variables (eg for cell volume) and supply `convertfn` and `convertinfo` 
to initialize to `:initial_value*convertfn(convertvars_tuple, i)` where the argument of `convertfn`
is the Tuple generated by `VarList_tuple(convertvars)`.

Example: To interpret :initial_value as a concentration-like quantity:

    convertvars = [volume], 
    convertfn = ((volume, ), i) -> volume[i], 
    convertinfo = " * volume"

"""
function add_method_setup_initialvalue_vars_default!(
    react::AbstractReaction, variables;
    filterfn=var -> get_attribute(var, :vfunction) in (VF_StateExplicit, VF_State, VF_Total, VF_Constraint),
    transfer_attribute_vars = [],
    convertvars = [],
    convertfn = (convertvars_tuple, i) -> 1.0, 
    convertinfo = "",
    setup_callback = (method, attribute_value, var, vardata) -> nothing,
)
    isempty(transfer_attribute_vars) || length(variables) == length(transfer_attribute_vars) ||
        error("transfer_attribute_vars list is not the same length as variables list")

    setup_vars = []
    setup_transfer_vars = []
   
    for i in eachindex(variables)
        var = variables[i]
        if filterfn(var)
            push!(setup_vars, var)
            isempty(transfer_attribute_vars) || 
                push!(setup_transfer_vars, VarDep(transfer_attribute_vars[i]))
        end
    end
  
    add_method_setup!(
        react,
        setup_initialvalue_vars_default,
        (VarList_fields(setup_vars), VarList_tuple(setup_transfer_vars), VarList_tuple(convertvars)),
        p = (convertfn, convertinfo, setup_callback),
    )
    
    return nothing
end

"""
    setup_initialvalue_vars_default

Initialize Variables to (`:initial_value` or `:norm_value`) [* convertfn] at beginning of integration.
"""
function setup_initialvalue_vars_default(
    m::ReactionMethod, 
    (varfields, transfervardata, convertvarsdata),
    cellrange::AbstractCellRange,
    attribute_value
)    
    (convertfn, convertinfo, setup_callback) = m.p

    # VariableReactions corresponding to (vardata, transfervardata, convertvarsdata)
    vars, transfer_attribute_vars, convertvars = get_variables_tuple(m)

    domvars = [v.linkvar for v in vars] # get_attributes from Domain var, as only 'master' Reaction var will be updated
    if isempty(transfer_attribute_vars)
        attrbvars = domvars
    else
        attrbvars = [v.linkvar for v in transfer_attribute_vars]
    end

    for (rv, v, attrbv, vfield) in zip(vars, domvars, attrbvars, varfields)
      
        if v === attrbv
            trsfrinfo = ""
        else
            trsfrinfo = " [from $(fullname(attrbv))]"
        end

        init_field!(
            vfield, attribute_value, attrbv, convertfn, convertvarsdata, cellrange, (fullname(v), convertinfo, trsfrinfo)
        )
       
        setup_callback(m, attribute_value, rv, vfield.values)
    end

    return nothing
end

"""
    add_method_initialize_zero_vars_default!(react::AbstractReaction, variables=PB.get_variables(react))

Create and add a default method to initialize Variables to zero at beginning of each timestep. Defaults to
adding all Variables from `react` with `:initialize_to_zero` attribute `true`.

NB: TODO variables are converted to VarDep (no dependency checking or sorting needed, could define a VarInit or similar?)
"""
function add_method_initialize_zero_vars_default!(
    react::AbstractReaction, variables=get_variables(react); 
    filterfn=v->get_attribute(v, :initialize_to_zero)
)

    init_vars = []
    
    for var in variables        
        if filterfn(var)
            push!(init_vars, VarInit(var))
        end
    end

    add_method_initialize!(
        react,
        initialize_zero_vars_default,
        (
            VarList_fields(init_vars),
        ),
    )

    return nothing
end

"""
    initialize_zero_vars_default

Initialize variables at beginning of timestep
"""
function initialize_zero_vars_default(
    m::ReactionMethod, 
    (init_fields,), 
    cellrange::AbstractCellRange,
    deltat
)
    
    IteratorUtils.foreach_longtuple_p(zero_field!, init_fields, cellrange)

    return nothing
end



"""
    add_method_do_nothing!(react::AbstractReaction, variables)

Create and add a dummy method to define Variables only.
"""
function add_method_do_nothing!(react::AbstractReaction, variables)

    add_method_do!(
        react,
        methodfn_do_nothing,
        (VarList_tuple(variables),)
    )
    
    return nothing
end

methodfn_do_nothing(m::ReactionMethod, (vars, ), cellrange::AbstractCellRange, deltat) = nothing

is_do_nothing(m::ReactionMethod) = false
is_do_nothing(m::ReactionMethod{M, R, P, V}) where {M <: typeof(methodfn_do_nothing), R, P, V} = true


function do_method_barrier(m::ReactionMethod, _, _, _)
    (barrier, barrierfn) = m.p

    barrierfn(barrier)

    return nothing
end

"""
    reaction_method_thread_barrier(barrier)

Create a ReactionMethod holding a thread barrier `barrier`.
"""
function reaction_method_thread_barrier(barrier, barrierfn)
    
    return ReactionMethod(
        do_method_barrier,
        NoReaction(),
        "thread_barrier",
        (VarList_nothing(),),
        (barrier, barrierfn), # p 
        [-1], 
        Domain(name="thread_barrier", ID=-1, parameters=Dict{String, Any}()),
    )
end



is_do_barrier(m::ReactionMethod) = false
is_do_barrier(m::ReactionMethod{M, R, P, V}) where {M <: typeof(do_method_barrier), R, P, V} = true
