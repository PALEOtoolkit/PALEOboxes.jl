
"""
    AtomicScalar{T}

Provides an implementation of the Julia Array interface
for a length 1 Vector that holds an atomic scalar variable.
"""
mutable struct AtomicScalar{T} <: AbstractArray{T, 1} 
    value::Threads.Atomic{T}
    function AtomicScalar{T}() where T
        return new(Threads.Atomic{T}())
    end
end

Base.size(as::AtomicScalar) = (1, )
Base.IndexStyle(::Type{<:AtomicScalar}) = IndexLinear()

Base.getindex(as::AtomicScalar) = as.value[]
function Base.getindex(as::AtomicScalar, i::Integer)
    @boundscheck(checkbounds(as, i))
    return as.value[]
end

Base.setindex!(as::AtomicScalar, v) = (as.value[] = v) 
function Base.setindex!(as::AtomicScalar, v, i::Integer)
    @boundscheck(checkbounds(as, i))
    as.value[] = v
    return v
end

"""
    atomic_add!(as::AtomicScalar, v)
    atomic_add!(a::Array, v)

Add v to as[] using Thread-safe atomic operation, with fallback for
normal Julia Arrays.

Fallback adds v to a[] using standard addition.
"""
atomic_add!(as::AtomicScalar, v) = Threads.atomic_add!(as.value, v)

@Base.propagate_inbounds atomic_add!(a::Array, v) = (a[] += v)