


"""
DomainData

struct to hold VariableDomain data.
"""
Base.@kwdef mutable struct DomainData <: AbstractDomainData
    domain::Domain
    variable_data::Vector{Any}           = Vector{Any}()
end

function DomainData(domain)
    domaindata = DomainData(domain = domain, variable_data = fill(nothing, get_num_variables(domain)))
 
    return domaindata
end



"""
    ModelData(model::Model; arrays_eltype::DataType=Float64, allocatenans::Bool=true)

Create a ModelData struct containing `model` data arrays.

One set of data arrays is created with `eltype=arrays_eltype`, accessed with `arrays_idx=1`

Additional sets of data arrays may be added by [`push_arrays_data!`](@ref), eg in order to support automatic differentiation
which requires Dual numbers as the array element type.

# Fields
- `cellranges_all::Vector{AbstractCellRange}`: default cellranges covering all domains
- `dispatchlists_all`: default dispatchlists covering all domains
- `solver_view_all`: optional untyped context field for use by external solvers.
"""
Base.@kwdef mutable struct ModelData <: AbstractModelData
    model::Model
    domain_data::Vector{Tuple{DataType, String, Vector{DomainData}}} = Tuple{DataType, String, Vector{DomainData}}[]
    allocatenans::Bool # fill with NaN when allocating

    sorted_methodsdata_setup                    = nothing    # only for arrays_idx = 1
    sorted_methodsdata_initialize               = [] # one per arrays_idx
    sorted_methodsdata_do                       = [] # one per arrays_idx
    
    cellranges_all::Vector{AbstractCellRange}   = Vector{AbstractCellRange}()    
    dispatchlists_all                           = nothing    # only for arrays_idx = 1
    solver_view_all                             = nothing    # only for arrays_idx = 1
end

function ModelData(
    model::Model; 
    arrays_eltype::DataType=Float64, 
    allocatenans::Bool=true,
)
    modeldata = ModelData(;
        model,
        allocatenans,        
    )

    push_arrays_data!(modeldata, arrays_eltype, "base")

    return modeldata
end

"""
    push_arrays_data!(modeldata, arrays_eltype::DataType, arrays_tagname::AbstractString)

Add an (unallocated) additional array set with element type `arrays_eltype`.
"""
function push_arrays_data!(modeldata::ModelData, arrays_eltype::DataType, arrays_tagname::AbstractString)
    isnothing(find_arrays_idx(modeldata, arrays_tagname::AbstractString)) || throw(ArgumentError("arrays_tagname $arrays_tagname already exists"))

    push!(modeldata.domain_data, (arrays_eltype, arrays_tagname,  [DomainData(dom) for dom in modeldata.model.domains]))

    resize!(modeldata.sorted_methodsdata_initialize, num_arrays(modeldata))
    resize!(modeldata.sorted_methodsdata_do, num_arrays(modeldata))

    return nothing
end

"""
    pop_arrays_data!(modeldata)

Remove last array set (NB: base array set with `arrays_idx==1`) cannot be removed).
"""
function pop_arrays_data!(modeldata::ModelData)
    length(modeldata.domain_data) > 1 || error("cannot remove base domain_data")
    return pop!(modeldata.domain_data)
end

num_arrays(modeldata::ModelData) = length(modeldata.domain_data)
find_arrays_idx(modeldata::ModelData, arrays_eltype::DataType) = findfirst(x -> x[1] == arrays_eltype, modeldata.domain_data)
find_arrays_idx(modeldata::ModelData, arrays_tagname::AbstractString) = findfirst(x -> [2] == arrays_tagname, modeldata.domain_data)
   
"""
    copy_base_values!(modeldata::ModelData, array_indices=2:num_arrays(modeldata))

Copy array contents from base `arrays_idx=1` to other `array_indices`
"""
function copy_base_values!(modeldata::ModelData, array_indices=2:num_arrays(modeldata))

    base_domain_data = first(modeldata.domain_data)
    _, _, base_data = base_domain_data

    for ai in array_indices
        domain_data = modeldata.domain_data[ai]
        dd_eltype, dd_tagname, dd_data = domain_data
        @info "copy_base_values! copying to modeldata arrays_index $ai ($dd_eltype, $dd_tagname)"
        for (dbase, d) in IteratorUtils.zipstrict(base_data, dd_data)
            for (vardata_base, vardata) in IteratorUtils.zipstrict(dbase.variable_data, d.variable_data)
                if !isnothing(vardata) # non-base arrays may not be allocated
                    if vardata !== vardata_base # non-base array may just be a "link" to the same array
                        vardata.values .= vardata_base.values
                    end
                end
            end
        end
    end

    return nothing
end

"""
    eltype(modeldata::ModelData, arrays_idx::Int) -> DataType
    eltype(modeldata::ModelData, arrays_tagname::AbstractString)

Determine the type of Array elements for ModelData
"""
function Base.eltype(modeldata::ModelData, arrays_idx::Int)
    arrays_idx in 1:length(modeldata.domain_data) || throw(ArgumentError("arrays_idx $arrays_idx out of range"))
    return modeldata.domain_data[arrays_idx][1]
end

function Base.eltype(modeldata::ModelData, arrays_tagname::AbstractString)
    arrays_idx = find_arrays_idx(modeldata, arrays_tagname)
    !isnothing(arrays_idx) || throw(ArgumentError("arrays_tagname $arrays_tagname not found"))
    return eltype(modeldata, arrays_idx)
end

Base.eltype(modeldata::ModelData) = error("eltype(modeldata) not supported")

"Check ModelData consistent with Model"
function check_modeldata(model::Model, modeldata::ModelData)
    if model === modeldata.model
        return nothing
    else
        error("modeldata inconsistent with model - check for missing call to create_modeldata")
    end
end

# compact form
function Base.show(io::IO, md::ModelData)
    print(io, "ModelData(model=", md.model, ", domain_data=", [(x[1], x[2]) for x in md.domain_data], ")")
end
# multiline form
function Base.show(io::IO, ::MIME"text/plain", md::ModelData)
    println(io, "ModelData")
    println(io, "  model=", md.model)
    println(io, "  domain_data=", [(x[1], x[2]) for x in md.domain_data])
    println(io, "  solver_view_all=", md.solver_view_all)
end


function get_domaindata(modeldata::ModelData, domain::AbstractDomain, arrays_idx::Int)
    arrays_idx in 1:length(modeldata.domain_data) || throw(ArgumentError("arrays_idx $arrays_idx out of range"))
    return modeldata.domain_data[arrays_idx][3][domain.ID]
end





