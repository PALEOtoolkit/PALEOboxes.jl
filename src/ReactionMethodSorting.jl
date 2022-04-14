
"""
    MethodSort

Sorted ReactionMethods.

Each group of methods in `groups` has no internal dependencies, but depends on methods
in previous groups in the Vector.

# Fields
- `methods::Vector{AbstractReactionMethod}`: (unsorted) ReactionMethods
- `groups::Vector{BitSet}`: indices in `methods` of sequentially dependent groups of ReactionMethods
  (ReactionMethods in each group depend on the previous group, no dependencies within group)
- `g::Graphs.SimpleDiGraph `: dependency graph (Vertices are numbered using indices in `methods`)
- `depvars::Dict{Tuple{Int64, Int64}, Vector{PALEOboxes.VariableDomain}}`:
  Dict((modifier_idx, dependent_idx)=>var):  Variables modified by method modifier_idx and
  required by method dependent_idx
"""
struct MethodSort
    methods::Vector{AbstractReactionMethod}
    groups::Vector{BitSet}
    g::Graphs.SimpleDiGraph
    depvars::Dict{Tuple{Int64, Int64}, Vector{PALEOboxes.VariableDomain}}
end

"""
    get_methods(ms::MethodSort; [method_barrier=nothing]) -> Vector{AbstractReactionMethod}

Flatten MethodSort to a Vector of ReactionMethods, optionally inserting
`method_barrier` between groups for a multithreaded timestepper.

NB: adds a `method_barrier` before first group, no `method_barrier` after last group.
"""
function get_methods(ms::MethodSort; method_barrier=nothing)
    methods = AbstractReactionMethod[]

    for grp in ms.groups
        isnothing(method_barrier) || push!(methods, method_barrier)
        append!(methods, (ms.methods[midx] for midx in grp))
    end
    return methods
end

function Base.show(io::IO, m::MethodSort)
    println(io, "MethodSort")
    for grp in m.groups
        println(io, "  group=", grp)
        for midx in grp
            println(io, "    ", fullname(m.methods[midx]))
        end
    end
end


"""
    dfs_methods(methods::Vector, all_domvars::Vector; verbose=false) -> MethodSort

Depth-first topological sort based on dependencies from Variables.

NB: No grouping of methods (returns each method in a group of 1).
"""
function dfs_methods(methods::Vector, all_domvars::Vector; verbose=false)

    g, depvars = _method_dependencies(methods, all_domvars)

    # topological sort
    groups = BitSet.(Graphs.topological_sort_by_dfs(g))

    MethodSort(methods, groups, g, depvars)
end

"""
    group_methods(methods::Vector, all_domvars::Vector; verbose=false) -> MethodSort

Group methods based on dependencies from Variables using a topological sort
(using Kahn's algorithm, similar to a breadth-first search).
Methods in each group have no dependencies (ie can run in any order). Each group must
run after the previous groups in the list.
"""
function group_methods(methods::Vector, all_domvars::Vector; verbose=false)

    g, depvars = _method_dependencies(methods, all_domvars)

    num_vertices = Graphs.nv(g)

    # number of incoming edges for each Vertex
    inn = [length(Graphs.inneighbors(g,v)) for v in Graphs.vertices(g)]

    groups = BitSet[]

    num_grouped = 0
    while true
        # find all vertices with no remaining incoming edges
        grp = BitSet(v for v in Graphs.vertices(g) if iszero(inn[v]))

        if isempty(grp)
            break
        else
            push!(groups, grp)
            num_grouped += length(grp)
            verbose && println(grp)
            for gv in grp
                verbose && println("    ", fullname(methods[gv]))
                inn[gv] -= 1 # now -1 (ie remove from further consideration)
                # remove an edge from all our dependents
                for ov in Graphs.outneighbors(g, gv)
                    inn[ov] -= 1
                end
            end
        end
    end

    num_grouped == num_vertices ||
        error("The input graph contains at least one loop.")

    return MethodSort(methods, groups, g, depvars)
end


"""
    _method_dependencies(methods::Vector, all_domvars::Vector) ->
        (g::Graphs.SimpleDiGraph, depvars::Dict((modifier_idx, dependent_idx)=>Vector{VariableDomain})

Create ReactionMethod dependency graph based on dependencies from `all_domvars`.

Vertices of directed graph `g` are numbered by index in `methods`.  'Out' edges from a method point to
(are `in` edges to) methods that are dependent on Variables modified by it.

Domain Variables generating dependencies are returned in `depvars::Dict((modifier_idx, dependent_idx)=>var)`
where `modifier_idx`, `dependent_idx` are indices in `methods`.
"""
function _method_dependencies(methods::Vector, all_domvars::Vector)

    # invert to get map from ReactionMethod  -> ReactionMethod index
    method_to_idx = Dict{AbstractReactionMethod, Int}(
        method=>idx for (idx, method) in enumerate(methods)
    )

    g = Graphs.SimpleDiGraph(length(method_to_idx))
    depvars = Dict{Tuple{Int64, Int64}, Vector{PALEOboxes.VariableDomain}}()

    # add edges for each variable Dependency
    for var in all_domvars
        modifying_methods = get_modifying_methods(var) # all methods (setup, initialize, do)
        dependent_methods = get_dependent_methods(var) # all methods (setup, initialize, do)
        for mod_m in modifying_methods
            for dep_m in dependent_methods
                # filter out methods to those we are interested in
                if haskey(method_to_idx, mod_m) && haskey(method_to_idx, dep_m)
                    mod_idx = method_to_idx[mod_m]
                    dep_idx = method_to_idx[dep_m]

                    # Add edge to graph - duplicate edges will not be added and will return false
                    Graphs.add_edge!(g, mod_idx, dep_idx)

                    # keep track of all Variables generating dependencies
                    push!(get!(depvars, (mod_idx, dep_idx), VariableDomain[]), var)
                end
            end
        end
    end

    return (g, depvars)
end
