import Printf
import InteractiveUtils

"""
    find_all_reactions() -> Dict{String, Type}

Use `InteractiveUtils.subtypes(AbstractReaction)` to find all currently loaded subtypes off
AbstractReaction, and create a `Dict` with last part of the name of the Type
as key (ie without the module prefix) and Type as value.

Any Types that generate non-unique keys (eg Module1.MyReactionType and Module2.MyReactionType) will generate
a warning, and no entry will be added to the Dict (so if this Reaction is present in a config file, it will
not be found and will error).
"""
function find_all_reactions()
    rtypes = InteractiveUtils.subtypes(AbstractReaction)

    rdict = Dict{String, Type}()
    duplicate_keys = []
    for ReactionType in rtypes
        rname = _classname(ReactionType)
        if haskey(rdict, rname)
            push!(duplicate_keys, (rname, ReactionType))
        end
        rdict[rname] = ReactionType
    end

    for (rname, ReactionType) in duplicate_keys
        @warn "Duplicate reaction name $rname for Type $ReactionType (removing from Dict)"
        if haskey(rdict, rname)
            @warn "Duplicate reaction name $rname for Type $(rdict[rname]) (removing from Dict)"
            delete!(rdict, rname)
        end
    end

    return rdict
end

"""
    _classname(ReactionType::Type{<:AbstractReaction}) -> String

Get Reaction classname from `ReactionType` (this is the Julia Type after stripping module prefixes and converting to String)
"""
_classname(ReactionType::Type{<:AbstractReaction}) = String(last(split(string(ReactionType), ".")))

"""
    find_reaction(class::AbstractString) -> ReactionType

Look up "class" in list of Reactions from [`find_all_reactions`](@ref), and return 
fully-qualified Reaction Type (including module prefixes).
"""
function find_reaction(class::AbstractString)
    rdict = find_all_reactions()
    if haskey(rdict, class)
        return rdict[class]
    else
        error("class \"$class\" not found")
    end
end

"""
    create_reaction(ReactionType::Type{<:AbstractReaction}, base::ReactionBase) -> reaction::AbstractReaction

Default method to create a `ReactionType` and set `base` field.

A reaction implementation may optionally implement a custom method eg to set additional fields
"""
function create_reaction(ReactionType::Type{<:AbstractReaction}, base::ReactionBase)
    return ReactionType(base=base)
end

"""
    create_reaction(rdict::Dict{String, Type}, classname::String, name::String, external_parameters::Dict{String, Any}) -> reaction::AbstractReaction
    create_reaction(ReactionType::Type{<:AbstractReaction}, name::String, external_parameters::Dict{String, Any}) -> reaction::AbstractReaction

Create and configure a reaction.

Sets `ReactionBase` with name, classname, external_parameters, and list of `Parameters` from `pars` field (if present)
""" 
function create_reaction(
    ReactionType::Type{<:AbstractReaction}, name::String, external_parameters::Dict{String, Any};
    classname=_classname(ReactionType),
)
    base=ReactionBase(;name, classname, external_parameters)
    rj = create_reaction(ReactionType, base)
    # Add parameters from pars field
    if hasproperty(rj, :pars)
        add_par(rj, rj.pars)
    end
    return rj
end

function create_reaction(
    rdict::Dict{String, Type}, classname::String, name::String, external_parameters::Dict{String, Any}
)   
    if haskey(rdict, classname)
        return create_reaction(rdict[classname], name, external_parameters; classname)
    else
        error("create_reaction: name $name classname $classname not found")
        return nothing
    end
end


function add_reaction_factory(ReactionType::Type{<:AbstractReaction})
    Base.depwarn("call to deprecated add_reaction_factory($ReactionType), this does nothing and can be removed", :add_reaction_factory, force=true)
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
    
    for (classname, ReactionType) in sort!(OrderedCollections.OrderedDict(find_all_reactions()))
        if occursin(classnamefilter, classname)
            println(classname)
            
            rtstring = Printf.@sprintf("%s", ReactionType)
            if  occursin(typenamefilter, rtstring )      
                println("    ", ReactionType)
                
                # doc = Base.Docs.doc(rt)  # a Markdown object?
                # display(doc)
            end
        end
    end

end