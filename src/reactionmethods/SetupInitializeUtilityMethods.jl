##############################
# Default setup methods (optional)
##############################

"""
    add_method_setup_initialvalue_vars_default!(react::AbstractReaction, variables [; kwargs...])

Create and add a default method to initialize Variables matching `filterfn` (defaults to state Variables)
at beginning of integration.

# Setup callbacks used
- State Variables and similar (`:vfunction != VF_Undefined`) are initialized in a setup callback with 
  `attribute_name in (:initial_value, :norm_value)`, with values from those Variable attributes.
- If `force_state_norm_value=false`, other Variables (with `:vfunction == VF_Undefined`) are initialized in a setup callback with 
  `attribute_name=:setup`, with values from the `:initial_value` Variable attribute. NB: `filterfn` must be set to include these Variables.
- If `force_initial_norm_value=true`, all Variables (including those with `:vfunction == VF_Undefined`) are initialised as state Variables
  
# Keywords
- `filterfn`: set to f(var)::Bool to override the default selection for state variables only
  (Variables with `:vfunction in (VF_StateExplicit, VF_State, VF_Total, VF_StateTotal, VF_Constraint)`)
- `force_initial_norm_value=false`: `true` to always use `:initial_value`, `:norm_value`, even for variables with `:vfunction=VF_Undefined`
- `transfer_attribute_vars=[]`: Set to a list of the same length as `variables` to initialise `variables`
  from attributes of `transfer_attribute_vars`.
- `setup_callback=(method, attribute_name, var, vardata) -> nothing`: Set to a function that is called 
  after each Variable initialisation eg to store norm values.
- `convertvars=[]`
- `convertfn = (convertvars_tuple, i) -> 1.0`
- `convertinfo = ""`

# Including volume etc conversion
Set `convertvars` to a Vector of Variables (eg for cell volume) and supply `convertfn` and `convertinfo` 
to initialize to `:initial_value*convertfn(convertvars_tuple, i)` where the argument of `convertfn`
is the Tuple generated by `VarList_tuple(convertvars)`.

Example: To interpret `:initial_value` as a concentration-like quantity:

    convertvars = [volume], 
    convertfn = ((volume, ), i) -> volume[i], 
    convertinfo = " * volume"

"""
function add_method_setup_initialvalue_vars_default!(
    react::AbstractReaction, variables;
    filterfn=var -> get_attribute(var, :vfunction) in (VF_StateExplicit, VF_State, VF_Total,VF_StateTotal, VF_Constraint),
    force_initial_norm_value=false,
    transfer_attribute_vars = [],
    convertvars = [],
    convertfn = (convertvars_tuple, i) -> 1.0, 
    convertinfo = "",
    setup_callback = (method, attribute_name, var, vardata) -> nothing,
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
        p = (force_initial_norm_value, convertfn, convertinfo, setup_callback),
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
    attribute_name
)    
    (force_initial_norm_value, convertfn, convertinfo, setup_callback) = m.p
   
    # VariableReactions corresponding to (vardata, transfervardata, convertvarsdata)
    vars, transfer_attribute_vars, convertvars = get_variables_tuple(m)

    domvars = [v.linkvar for v in vars] # get_attributes from Domain var, as only 'master' Reaction var will be updated
    if isempty(transfer_attribute_vars)
        attrbvars = domvars
    else
        attrbvars = [v.linkvar for v in transfer_attribute_vars]
    end

    first_var = true
    for (rv, v, attrbv, vfield) in IteratorUtils.zipstrict(vars, domvars, attrbvars, varfields)        
        vfunction = get_attribute(v, :vfunction, VF_Undefined)
        if vfunction == VF_Undefined && !force_initial_norm_value
            # non-state Variables are initialized in :setup, from :initial_value attribute
            attribute_name == :setup || continue
            attribute_to_read = :initial_value
        else
            # state Variables etc are initialized in :initial_value, :norm_value, from that attribute
            attribute_name in (:initial_value, :norm_value) || continue
            attribute_to_read = attribute_name
        end

        if first_var
            @info "$(fullname(m)):"
            first_var = false
        end

        if v === attrbv
            trsfrinfo = ""
        else
            trsfrinfo = " [from $(fullname(attrbv))]"
        end

        init_field!(
            vfield, attribute_to_read, attrbv, convertfn, convertvarsdata, cellrange, (fullname(v), convertinfo, trsfrinfo)
        )
       
        setup_callback(m, attribute_name, rv, vfield.values)
    end

    return nothing
end

##############################
# Default initialize methods (optional)
##############################

"""
    add_method_initialize_zero_vars_default!(react::AbstractReaction, variables=PB.get_variables(react))

Create and add a default method to initialize Variables to zero at beginning of each timestep. Defaults to
adding all Variables from `react` with `:initialize_to_zero` attribute `true`.

NB: TODO variables are converted to VarDep (no dependency checking or sorting needed, could define a VarInit or similar?)
"""
function add_method_initialize_zero_vars_default!(
    react::AbstractReaction, variables=get_variables(react); 
    filterfn=v->get_attribute(v, :initialize_to_zero, false)
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

################################
# do nothing method
#################################

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
is_do_nothing(m::ReactionMethod{typeof(methodfn_do_nothing)}) = true

###############################################
# Thread barrier method (PALEOboxes internal use only)
###############################################

function do_method_barrier(m::ReactionMethod, _, _, _)
    (barrier, barrierfn) = m.p

    barrierfn(barrier)

    return nothing
end

"""
    reaction_method_thread_barrier(barrier, barrierfn; [operatorID=[1]])

Create a ReactionMethod holding a thread barrier `barrier`.
"""
function reaction_method_thread_barrier(barrier, barrierfn; operatorID=[1])
    
    return ReactionMethod(
        do_method_barrier,
        NoReaction(),
        "thread_barrier",
        (VarList_nothing(),),
        (barrier, barrierfn), # p 
        operatorID, 
        Domain(name="thread_barrier", ID=-1, parameters=Dict{String, Any}()),
    )
end



is_do_barrier(m::ReactionMethod) = false
is_do_barrier(m::ReactionMethod{typeof(do_method_barrier)}) = true
