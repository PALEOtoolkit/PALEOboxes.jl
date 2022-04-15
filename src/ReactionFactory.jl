import Printf

"""
    const reaction_factories = Dict{classname::String, Function}()

Dictionary of reaction factories.

A module containing Reactions should define a module `__init__()` function that includes a call
to [`add_reaction_factory`] for each Reaction. This registers functions that then return instances of the Reaction structs.

    "Install create_reactionXXX factories when module imported"
    function __init__()
        PB.add_reaction_factory(ReactionReservoirScalar)
    end

"""
const reaction_factories = Dict{String, Function}()

"""
    reaction_factory(ReactionType::Type{<:AbstractReaction}; [, operatorID=[1]] [,extra_par_fields], [, set_pars::Tuple]) -> Function

Return a function that will create Reaction of type `rt` when called.
    
Optionally, after creation, set and freeze Parameters in `set_pars` (a Tuple of Pairs of name=>value eg ("total"=>false,) )
"""
function reaction_factory(
    ReactionType::Type{<:AbstractReaction}; 
    operatorID=[1],
    extra_par_fields::Vector{Symbol}=Symbol[], # eg [:pars_stoich]
    set_pars::Tuple{Vararg{Pair}}=()
)

    function create_init_reaction(base::ReactionBase)
        rj = ReactionType(base=base)
        rj.base.operatorID = operatorID
       
        # Add parameters from pars field
        if hasfield(ReactionType, :pars)
            add_par(rj, rj.pars)
        end

        # Add parameters from any additional fields
        for ef in extra_par_fields
            add_par(rj, getfield(rj, ef))
        end

        # optionally set and freeze supplied parameters
        for (parname, value) in set_pars
            par = get_parameter(rj, parname)
            setvalueanddefault!(par, value)
            setfrozen!(par)
        end

        return rj
    end

    return create_init_reaction
end

"""
    add_reaction_factory(ReactionType::Type{<:AbstractReaction}; [, operatorID=[1]] [, set_pars::Tuple])
    add_reaction_factory(classname, ReactionType::Type{<:AbstractReaction}; [, operatorID=[1]] [, set_pars::Tuple])
    add_reaction_factory(classname, factory_fn)

Register a ReactionType to be created either by using the default [`reaction_factory`](@ref) and typename as classname,
or by the explicitly supplied `factory_fn` and `classname`.

Optionally, after creation, set and freeze Parameters in `set_pars` (a Tuple of Pairs of name=>value eg ("total"=>false,) )
"""
function add_reaction_factory(ReactionType::Type{<:AbstractReaction}; kwargs...)
    add_reaction_factory(string(nameof(ReactionType)), ReactionType; kwargs...)
    return nothing
end

function add_reaction_factory(classname::AbstractString, ReactionType::Type{<:AbstractReaction}; kwargs...)
    add_reaction_factory(classname, reaction_factory(ReactionType; kwargs...))
    return nothing
end

function add_reaction_factory(classname::AbstractString, factory_fn::Function)
    reaction_factories[classname] = factory_fn
    return nothing
end

"""
    _create_reaction(classname::String, name::String) -> AbstractReaction

Create a new instance of a reaction
"""
function _create_reaction(classname::String, name::String, external_parameters::Dict{String, Any})

    # Create a new ReactionXXX instance if a Reaction factory has been registered for this classname
    if haskey(reaction_factories, classname)
        base=ReactionBase(name=name, classname=classname, external_parameters=external_parameters)
        # use factory to create ReactionXXX instance
        rj = reaction_factories[classname](base)
        if ! (rj isa AbstractReaction)
            error("_create_reaction classname='$classname', name='$name' returned ", rj)
        end
        return rj
    else
        error("classname $classname not found")
        return nothing
    end
end

"""
    show_all_reactions(classfilter="", typenamefilter="")

List all registered Reactions with `classname` containing (`occursin`) `classfilter` and `typenamefilter`
A Reaction is loaded when the module that defines it is imported.

Examples:
- `PB.show_all_reactions(r"reservoir"i)` case-insensitive match for classname containing "reservoir".
- `PB.show_all_reactions("", "Reservoir")` all Reactions defined in a module name containing "Reservoir"
"""
function show_all_reactions(classnamefilter="", typenamefilter="")
    
    for (classname, ctorfn) in sort(reaction_factories)
        if occursin(classnamefilter, classname)
            println(classname)
            
            base=ReactionBase(name="test", classname=classname, external_parameters=Dict{String, Any}())
            r = ctorfn(base)       
            rt = typeof(r)
            rtstring = Printf.@sprintf("%s", rt)
            if  occursin(typenamefilter, rtstring )      
                println("    ", rt)
                
                # doc = Base.Docs.doc(rt)  # a Markdown object?
                # display(doc)
            end
        end
    end

end