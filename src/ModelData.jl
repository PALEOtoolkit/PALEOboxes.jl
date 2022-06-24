


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
    ModelData{T}

Contains `model` data arrays of element type T.

# Fields
- `cellranges_all::Vector{AbstractCellRange}`: default cellranges covering all domains
- `dispatchlists_all`: default dispatchlists covering all domains
- `solver_view_all`: optional untyped context field for use by external solvers.
"""
Base.@kwdef mutable struct ModelData{T<:Real} <: AbstractModelData
    model::Model
    threadsafe::Bool
    allocatenans::Bool                          = true # fill with NaN when allocating
    domain_data::Vector{DomainData}             = Vector{DomainData}()

    sorted_methodsdata_setup                    = nothing    
    sorted_methodsdata_initialize               = nothing
    sorted_methodsdata_do                       = nothing
    
    cellranges_all::Vector{AbstractCellRange}   = Vector{AbstractCellRange}()    
    dispatchlists_all                           = nothing    
    solver_view_all                             = nothing
end

"Determine the type of Array elements for ModelData"
Base.eltype(::Type{ModelData{T}}) where{T} = T

"Check ModelData consistent with Model"
function check_modeldata(model::Model, modeldata::ModelData)
    if model === modeldata.model
        return nothing
    else
        error("modeldata inconsistent with model - check for missing call to create_modeldata")
    end
end


"compact form"
function Base.show(io::IO, md::ModelData)
    print(io, "ModelData(model=", md.model, "eltype=", eltype(md), ")")
end
"multiline form"
function Base.show(io::IO, ::MIME"text/plain", md::ModelData)
    println(io, "ModelData")
    println(io, "  model=", md.model)
    println(io, "  eltype=", eltype(md))
    println(io, "  solver_view_all=", md.solver_view_all)
end


function get_domaindata(modeldata::ModelData, domain::AbstractDomain)
    return modeldata.domain_data[domain.ID]
end





