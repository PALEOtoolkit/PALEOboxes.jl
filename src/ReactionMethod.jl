"""
    ReactionMethod(
        methodfn::Function,
        reaction::AbstractReaction,
        name::String,
        varlists::Tuple{Vararg{AbstractVarList}},
        p, 
        operatorID::Vector{Int64}, 
        domain::AbstractDomain; 
        preparefn = (m, vardata) -> vardata
    ) -> m::ReactionMethod

Defines a callback function `methodfn` with Variables `varlists`, to be called from the Model framework
either during setup or as part of the main loop.

# Fields
$(FIELDS)

# methodfn
The `methodfn` callback is:

    methodfn(m::ReactionMethod, pars, vardata::Tuple, cellrange::AbstractCellRange, modelctxt)

or (if Parameters are not required):

    methodfn(m::ReactionMethod, vardata::Tuple, cellrange::AbstractCellRange, modelctxt)

With arguments:
- `m::ReactionMethod`: context is available as `m.reaction::AbstractReaction` (the Reaction that defined the `ReactionMethod`),
  and `m.p` (an arbitrary extra context field supplied when `ReactionMethod` created).
- `pars`: a struct with Parameters as fields (current just the ParametersTuple defined as `reaction.pars`)
- `vardata`: A Tuple of collections of views on Domain data arrays corresponding to [`VariableReaction`](@ref)s defined by `varlists`
- `cellrange::AbstractCellRange`: range of cells to calculate.
- `modelctxt`:
    - for a setup method, `:setup`, `:initial_value` or `:norm_value` defining the type of setup requested
    - for a main loop method `deltat` providing timestep information eg for rate throttling.

# preparefn
An optional `preparefn` callback can be supplied eg to allocate buffers that require knowledge of
the data types of `vardata` or to cache expensive calculations:

    preparefn(m::ReactionMethod, vardata::Tuple) -> modified_vardata::Tuple

This is called after model arrays are allocated, and prior to setup.
"""
struct ReactionMethod{M, R, P, Nargs} <: AbstractReactionMethod    
    "callback from Model framework"
    methodfn::M   # NB: store function object, even though we have M (see https://discourse.julialang.org/t/is-there-a-way-to-get-f-from-typeof-f/18818/12)

    "the Reaction that created this ReactionMethod"
    reaction::R

    "a descriptive name, eg generated from the name of methodfn"
    name::String

    "Tuple of VarLists, each representing a list of [`VariableReaction`](@ref)s.
    Corresponding Variable accessors `vardata` (views on Arrays) will be provided to the `methodfn` callback.
    NB: not concretely typed to reduce compile time, as not performance-critical"
    varlists::Tuple{Vararg{AbstractVarList}}

    "optional context field (of arbitrary type) to store data needed by methodfn."
    p::P

    operatorID::Vector{Int64}
    domain::Domain

    "preparefn(m::ReactionMethod, vardata::Tuple) -> modified_vardata::Tuple
     optionally modify `vardata` to eg add buffers. NB: not concretely typed as not performance-critical"
    preparefn::Function

    function ReactionMethod(
        methodfn::M,
        reaction::R,
        name,
        varlists::Tuple{Vararg{AbstractVarList}},
        p::P, 
        operatorID::Vector{Int64}, 
        domain; 
        preparefn = (m, vardata) -> vardata,
    ) where {M <: Function, R <: AbstractReaction, P}        
        
        # Find number of arguments that methodfn takes
        # (in order to support two forms of 'methodfn', with and without Parameters)
        nargs = fieldcount(methods(methodfn).ms[1].sig) - 1

        newmethod = new{M, R, P, nargs}(
            methodfn, 
            reaction,
            name,
            varlists,
            p,
            operatorID, 
            domain, 
            preparefn,
        )

        for v in get_variables(newmethod)
            v.method = newmethod
            set_attribute!(v, :operatorID, newmethod.operatorID, allow_create=true)
        end

        return newmethod
    end
end

get_nargs(method::ReactionMethod{M, R, P, Nargs}) where {M, R, P, Nargs} = Nargs
get_nargs(methodref::Ref{ReactionMethod{M, R, P, Nargs}}) where {M, R, P, Nargs} = Nargs
# deprecated form without pars
@inline call_method(method::ReactionMethod{M, R, P, 4}, vardata, cr, modelctxt) where {M, R, P} = 
    method.methodfn(method, vardata, cr, modelctxt)
# updated form with pars
@inline call_method(method::ReactionMethod{M, R, P, 5}, vardata, cr, modelctxt) where {M, R, P} = 
    method.methodfn(method, method.reaction.pars, vardata, cr, modelctxt)    

@noinline call_method(methodref::Ref{<: ReactionMethod}, vardataref::Ref, cr, modelctxt) = 
    call_method(methodref[], vardataref[], cr, modelctxt)

# for benchmarking etc: apply codefn to the ReactionMethod methodfn (without this, will just apply to the call_method wrapper)
# codefn=code_warntype, code_llvm, code_native
call_method_codefn(io::IO, codefn, method::ReactionMethod{M, R, P, 4}, vardata, cr, modelctxt; kwargs...) where {M, R, P} = 
    codefn(io, method.methodfn, (typeof(method), typeof(vardata), typeof(cr), typeof(modelctxt)); kwargs...)
call_method_codefn(io::IO, codefn, method::ReactionMethod{M, R, P, 5}, vardata, cr, modelctxt; kwargs...) where {M, R, P} = 
    codefn(io, method.methodfn, (typeof(method), typeof(method.reaction.pars), typeof(vardata), typeof(cr), typeof(modelctxt)); kwargs...)


"""
    get_variables_tuple(method::AbstractReactionMethod) -> (Vector{VariableReaction}, ...)
    
Get all [`VariableReaction`](@ref)s from `method` as a Tuple of `Vector{VariableReaction}`
"""
get_variables_tuple(@nospecialize(method::ReactionMethod)) = Tuple(get_variables(vl) for vl in method.varlists)

"""
    get_variables(method::AbstractReactionMethod; filterfn = v -> true) -> Vector{VariableReaction}

Get VariableReactions from `method.varlists` as a flat Vector, optionally restricting to those that match `filterfn`
"""
function get_variables(
    @nospecialize(method::ReactionMethod);
    filterfn = v -> true
)
    vars = VariableReaction[]
    for vl in method.varlists
        append!(vars, filter(filterfn, get_variables(vl)))
    end
    return vars
end

"""
    get_variable(
        method::AbstractReactionMethod, localname::AbstractString; 
        allow_not_found=false
    ) -> v

Get a single VariableReaction `v` by `localname`.

If `localname` not present, returns `nothing` if `allow_not_found==true` otherwise errors.
"""
function get_variable(
    @nospecialize(method::ReactionMethod), localname::AbstractString; 
    allow_not_found=false
)
    matchvars = get_variables(method; filterfn = v -> v.localname==localname)
    length(matchvars) <= 1 || error("duplicate variable localname", localname)

    !isempty(matchvars) || allow_not_found || error("method ", fullname(method), " no variable localname=", localname)

    return isempty(matchvars) ? nothing : matchvars[1]
end

fullname(@nospecialize(method::ReactionMethod)) = fullname(method.reaction)*"."*method.name

is_method_setup(@nospecialize(method::ReactionMethod)) = (method in base(method.reaction).methods_setup)
is_method_initialize(@nospecialize(method::ReactionMethod)) = (method in base(method.reaction).methods_initialize)
is_method_do(@nospecialize(method::ReactionMethod)) = (method in base(method.reaction).methods_do)

"""
    get_rate_stoichiometry(m::PB.ReactionMethod) -> Vector[(ratevarname, process, Dict(speciesname=>stoich, ...)), ...]

DEPRECATED: Optionally provide rate Variable name(s), a process name ("photolysis", "reaction", ...) and stoichiometry of reactants and products, for postprocessing of results.

Use attributes on rate variable instead:
- `rate_processname::String`: a process name (eg "photolysis", "reaction", ...)
- `rate_species::Vector{String}` names of reactant and product species
- `rate_stoichiometry::Vector{Float64}` stoichiometry of reactant and product species

"""
get_rate_stoichiometry(@nospecialize(m::ReactionMethod)) = []

###########################################
# Pretty printing
############################################

"compact form"
function Base.show(io::IO, @nospecialize(method::ReactionMethod))
    print(
        io, 
        "ReactionMethod(fullname='", fullname(method), 
            "', methodfn=", string(method.methodfn), 
            ", domain='", method.domain.name ,
            "', operatorID=", method.operatorID,
            ")",
    )
end
