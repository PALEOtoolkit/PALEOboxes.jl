

"""
    CellRange

Defines a range of cells in a specified [`Domain`](@ref) as a linear list.

If `operatorID==0`, call all `Reaction`s, otherwise
only call those with matching `operatorID` (this enables operator splitting).
"""
Base.@kwdef mutable struct CellRange{T} <: AbstractCellRange
    domain::Domain
    operatorID::Int     = 0
    "may be any valid Julia indexing range thing eg 1:100, [1 2 3 4], etc"
    indices::T
end


"Add an array of indices to a CellRange instance"
function add_indices!(cellrange::CellRange{Vector{Int64}}, indicestoadd::Vector{Int64})   
    append!(cellrange.indices, indicestoadd)
     
    if length(unique(cellrange.indices)) != length(cellrange.indices)
        error("add_indices! duplicate indices")
    end   
end



"""
    CellRangeColumns

Defines a range of cells in a specified [`Domain`](@ref), organised by columns.
If `operatorID==0`, call all `Reaction`s, otherwise
only call those with matching `operatorID` (this enables operator splitting).
"""
Base.@kwdef mutable struct CellRangeColumns{T1, T2} <: AbstractCellRange
    domain::Domain
    operatorID::Int     = 0
    "iterator through cells in arbitrary order"
    indices::T1
    "iterator through columns: columns[n] returns a Pair icol=>cells where cells are ordered top to bottom"
    columns::T2
end

"replace a contiguous range of indices (as a Vector of indices) with a Range"
function replace_contiguous_range(indices)
    if indices == first(indices):last(indices)
        return first(indices):last(indices)
    else
        return indices
    end
end
