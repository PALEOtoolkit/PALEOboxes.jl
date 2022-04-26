import Printf
import InteractiveUtils


function find_all_reactions()
    rtypes = InteractiveUtils.subtypes(AbstractReaction)

    rdict = Dict{String, Type}()
    for ReactionType in rtypes
        rname = String(last(split(string(ReactionType), ".")))
        !haskey(rdict, rname) || 
            error("Duplicate reaction name $rname for Types $(rdict[rname])    $ReactionType")
        rdict[rname] = ReactionType
    end

    return rdict
end

"""
    create_reaction(ReactionType::Type{<:AbstractReaction}, base::ReactionBase) -> reaction::AbstractReaction

Create a `ReactionType` and set `base` field.

Default implementation may be overriden to eg set additional fields
"""
function create_reaction(ReactionType::Type{<:AbstractReaction}, base::ReactionBase)
    return ReactionType(base=base)
end

function create_reaction(rdict::Dict{String, Type}, classname::String, name::String, external_parameters::Dict{String, Any})
   
    if haskey(rdict, classname)
        base=ReactionBase(name=name, classname=classname, external_parameters=external_parameters)
        rj = create_reaction(rdict[classname], base)
        # Add parameters from pars field
        if hasproperty(rj, :pars)
            add_par(rj, rj.pars)
        end
        return rj
    else
        error("classname $classname not found")
        return nothing
    end
end


function add_reaction_factory(ReactionType::Type{<:AbstractReaction})
    @warn "call to deprecated add_reaction_factory($ReactionType)"
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
    
    for (classname, RType) in find_all_reactions()
        if occursin(classnamefilter, classname)
            println(classname)
            
            rtstring = Printf.@sprintf("%s", RType)
            if  occursin(typenamefilter, rtstring )      
                println("    ", RType)
                
                # doc = Base.Docs.doc(rt)  # a Markdown object?
                # display(doc)
            end
        end
    end

end