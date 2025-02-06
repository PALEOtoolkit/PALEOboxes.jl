
"""
     SIMDutils

Helper functions to pack and unpack vector data to enable use of SIMD instructions.
"""
module SIMDutils

using ...PALEOboxes: @public

@public FP64P2, FP64P4, FP32P4, FP32P8
@public SIMDIter, vgatherind, vscatterind!, vaddind!

import SIMD
import Preferences


#####################################
# Short names for packed vector Types
#####################################

const FP64P4 = SIMD.Vec{4, Float64}
const FP64P2 = SIMD.Vec{2, Float64}
const FP32P8 = SIMD.Vec{8, Float32}
const FP32P4 = SIMD.Vec{4, Float32}

######################################
# Iteration 
######################################

"""
    SIMDIter(baseiter, Val{N})
    SIMDIter(baseiter, ::Type{SIMD.Vec{N, U}})
    SIMDIter(baseiter, ::Val{1}) # scalar fallback
    SIMDIter(baseiter, ::Type{U}) where {U <: Real} # scalar fallback

Iterator that takes up to `N` SIMD elements at a time from `baseiter` (which should represent
indices into a Vector). See Julia package [SIMD.jl](https://github.com/eschnett/SIMD.jl)

If `baseiter` contained 1 or more but less then `N` elements, then `indices` is filled with
repeats of the last available element.

Returns Tuple of indices (length `N`).

# Examples


    v_a = [1.0, 2.0, 3.0, 4.0, 5.0]
    v_b = similar(v_a)

    iter = eachindex(v_a) # iter should represent indices into a Vector

    # simplest version - Float64 x 4, ie type of v_a x 4
 
    for i in SIMDIter(iter, Val(4))
        x = v_a[i]  # x is a packed SIMD vector
        v_b[i] = x
    end


    # with type conversion - Float32 x 8, ie explicitly change Type of SIMD vector

    ST = SIMD.Vec{8, Float32}}    
    for i in SIMDIter(iter, ST)
        #   v = vec[i]  <--> vgatherind(ST, vec, i)
        #   vec[i] = v  <--> vscatterind!(v, vec, i)
        #   vec[i] += v <--> vaddind!(v, vec, i) 

        x = vgatherind(ST, v_a, i)  # x is a packed SIMD vector with type conversion to Float32
        vscatterind!(x, v_b, i)
    end

"""
struct SIMDIter{T, N}
    baseiter::T
end

SIMDIter(baseiter::T, ::Val{N}) where {T, N} = SIMDIter{T, N}(baseiter)
SIMDIter(baseiter::T, ::Type{SIMD.Vec{N, U}}) where {T, N, U} = SIMDIter{T, N}(baseiter)
# scalar fallback
SIMDIter(baseiter, ::Val{1}) = baseiter
SIMDIter(baseiter, ::Type{U}) where {U <: Real} = baseiter

import Base.iterate
iterate(iter::SIMDIter{T, N}) where {T, N} = _iterate(iter, iterate(iter.baseiter))
iterate(iter::SIMDIter{T, N}, state) where {T, N} = _iterate(iter, iterate(iter.baseiter, state))

_iterate(iter::SIMDIter{T, 2}, basenext::Nothing) where {T} = nothing

function _iterate(iter::SIMDIter{T, 2}, basenext) where {T}

    local i1::Int64, i2::Int64

    (i1, state) = basenext

    basenext = iterate(iter.baseiter, state)
    if isnothing(basenext); i2 = i1;  else;  (i2, state) = basenext; end   

    return (SIMD.Vec((i1, i2)), state)   
end

_iterate(iter::SIMDIter{T, 4}, basenext::Nothing) where {T} = nothing

function _iterate(iter::SIMDIter{T, 4}, basenext) where {T}

    local i1::Int64, i2::Int64, i3::Int64, i4::Int64

    (i1, state) = basenext
    
    basenext = iterate(iter.baseiter, state)
    if isnothing(basenext); i2 = i1;  else;  (i2, state) = basenext; end   
    basenext = iterate(iter.baseiter, state)
    if isnothing(basenext); i3 = i2;  else;  (i3, state) = basenext; end
    basenext = iterate(iter.baseiter, state)
    if isnothing(basenext); i4 = i3;  else;  (i4, state) = basenext; end

    return (SIMD.Vec((i1, i2, i3, i4)), state)
end

_iterate(iter::SIMDIter{T, 8}, basenext::Nothing) where {T} = nothing

function _iterate(iter::SIMDIter{T, 8}, basenext) where {T}

    local i1::Int64, i2::Int64, i3::Int64, i4::Int64, i5::Int64, i6::Int64, i7::Int64, i8::Int64

    (i1, state) = basenext
    
    basenext = iterate(iter.baseiter, state)
    if isnothing(basenext); i2 = i1; else;  (i2, state) = basenext; end   
    basenext = iterate(iter.baseiter, state)
    if isnothing(basenext); i3 = i2; else;  (i3, state) = basenext; end
    basenext = iterate(iter.baseiter, state)
    if isnothing(basenext); i4 = i3; else;  (i4, state) = basenext; end
    basenext = iterate(iter.baseiter, state)
    if isnothing(basenext); i5 = i4; else;  (i5, state) = basenext; end
    basenext = iterate(iter.baseiter, state)
    if isnothing(basenext); i6 = i5; else;  (i6, state) = basenext; end
    basenext = iterate(iter.baseiter, state)
    if isnothing(basenext); i7 = i6; else;  (i7, state) = basenext; end
    basenext = iterate(iter.baseiter, state)
    if isnothing(basenext); i8 = i7; else;  (i8, state) = basenext; end

    return (SIMD.Vec((i1, i2, i3, i4, i5, i6, i7, i8)), state)
end



#################################
# Fill SIMD type from scalar Vector
#####################################

"""
    vgatherind(vec, indices::Integer, mask) -> v # v[i], scalar fallback
    vgatherind(vec, indices::SIMD.Vec{N, <:Integer}, mask) -> v::SIMD.Vec{N, eltype(vec)}

Fill SIMD type from scalar Vector, without type conversion

See [`SIMDIter`](@ref) for example usage.
"""
@Base.propagate_inbounds function vgatherind(vec, indices::Integer, mask)
    return vec[indices]
end

@Base.propagate_inbounds function vgatherind(vec, indices::SIMD.Vec{N, <:Integer}, mask) where N
    return SIMD.vgather(vec, indices, mask)
end

"""
    vgatherind(::Type{T}, vec, indices::Integer, mask) -> v::T # v[i], scalar fallback
    vgatherind(::Type{SIMD.Vec{N, T}}, vec, indices::SIMD.Vec{N, <:Integer}) -> v::SIMD.Vec{N, T}

Fill SIMD type from scalar Vector, with explicit type conversion to scalar type `T`

NB: assumes all `indices` are valid

See [`SIMDIter`](@ref) for example usage.
"""
@Base.propagate_inbounds function vgatherind(::Type{T}, vec, indices::Integer) where {T}
    return convert(T, vec[indices])
end

@Base.propagate_inbounds function vgatherind(::Type{SIMD.Vec{2, T}}, vec, indices::SIMD.Vec{2, <:Integer}) where T
    return SIMD.Vec{2, T}((vec[indices[1]], vec[indices[2]]))
end

@Base.propagate_inbounds function vgatherind(::Type{SIMD.Vec{4, T}}, vec, indices::SIMD.Vec{4, <:Integer}) where T
    return SIMD.Vec{4, T}((vec[indices[1]], vec[indices[2]], vec[indices[3]], vec[indices[4]]))
end

@Base.propagate_inbounds function vgatherind(::Type{SIMD.Vec{8, T}}, vec, indices::SIMD.Vec{8, <:Integer}) where T
    return SIMD.Vec{8, T}(
        ( 
            vec[indices[1]], vec[indices[2]], vec[indices[3]], vec[indices[4]],
            vec[indices[5]], vec[indices[6]], vec[indices[7]], vec[indices[8]]
        )
    )
end

###############################################
# Unpack SIMD type and write to scalar Vector
#############################################

"""
    vscatterind!(v::T, vec::Array{T, N}, indices::Integer, mask) # scalar fallback, no type conversion
    vscatterind!(v::T, vec::Array{T, N}, indices::Integer) # scalar fallback, no type conversion
    vscatterind!(v::SIMD.Vec{N, T}, vec, indices::SIMD.Vec{N, <:Integer}, mask) # converts to eltype{vec}
    vscatterind!(v::SIMD.Vec{N, T}, vec, indices::SIMD.Vec{N, <:Integer}) # converts to eltype{vec} 

Unpack SIMD type and write to scalar Vector, with type conversion to `eltype(vec)`.

Versions without `mask` assume all indices are valid (but may be repeated).

See [`SIMDIter`](@ref) for example usage.
"""
@Base.propagate_inbounds function vscatterind!(v::T, vec::Array{T, N}, indices::Integer, mask) where {T,N}
    vec[indices] = v   
    return nothing
end

@Base.propagate_inbounds function vscatterind!(v::T, vec::Array{T,N}, indices::Integer) where {T,N}
    vec[indices] = v   
    return nothing
end

@Base.propagate_inbounds function vscatterind!(v::SIMD.Vec{N, T}, vec, indices::SIMD.Vec{N, <:Integer}, mask) where {N, T}
    SIMD.vscatter(convert(SIMD.Vec{N, eltype(vec)}, v), vec, indices, mask)
    return nothing
end

@Base.propagate_inbounds function vscatterind!(v::SIMD.Vec{N, T}, vec, indices::SIMD.Vec{N, <:Integer}) where {N, T}
    SIMD.vscatter(convert(SIMD.Vec{N, eltype(vec)}, v), vec, indices)
    return nothing
end

###############################################
# Unpack SIMD type and add to scalar Vector
#############################################

"""
    vaddind!(v::T,  vec::Array{T, N}, indices::Integer) # scalar fallback, no type conversion
    vaddind!(v::SIMD.Vec{N, T}, vec, indices::SIMD.Vec{N, <:Integer})

Unpack SIMD type and add to scalar Vector, with type conversion to `eltype(vec)`.

Assumes all indices are valid (but may be repeated, which will still have the effect of one add).

See [`SIMDIter`](@ref) for example usage.
"""
@Base.propagate_inbounds function vaddind!(v::T,  vec::Array{T, N}, indices::Integer) where {T, N}
    vec[indices] += v
    return nothing
end

@Base.propagate_inbounds function vaddind!(v::SIMD.Vec{N, T}, vec, indices::SIMD.Vec{N, <:Integer}) where {N, T}
    tmp = vgatherind(SIMD.Vec{N, T}, vec, indices) 
    vscatterind!(tmp + v, vec, indices)
    return nothing
end


######################################################
# functions to call SLEEF library Sleef <https://sleef.org>, using binary from SLEEF_jll.jl
# see <https://github.com/chriselrod/SLEEFPirates.jl/blob/master/src/svmlwrap32.jl> for ccall examples
# 
# NB: need to provide specific Base.exp etc for each type in order to override the fallbacks in SIMD.jl
#############################################################

# raw tuple data types held in Vec.data
const FP64P4_d = NTuple{4,Core.VecElement{Float64}}
const FP64P2_d = NTuple{2,Core.VecElement{Float64}}
const FP32P8_d = NTuple{8,Core.VecElement{Float32}}
const FP32P4_d = NTuple{4,Core.VecElement{Float32}}

if !@Preferences.has_preference("USE_SLEEF")
    @info "$(@__MODULE__) defining USE_SLEEF = false in LocalPreferences.toml"
    @Preferences.set_preferences!("USE_SLEEF"=>false)
end

const USE_SLEEF = @Preferences.load_preference("USE_SLEEF", false)
@static if USE_SLEEF
    @info "$(@__MODULE__) defining Vectorized SIMD functions log, exp, log10 functions from Sleef library $(SLEEF_jll.libsleef)"*
        " - to disable Sleef, set USE_SLEEF = false in LocalPreferences.toml and restart your Julia session"
    
    import SLEEF_jll

    # exp Sleef Vectorized double/single precision base-e exponential functions functions with 1.0 ULP error bound
    sleefexp(v::FP64P4_d)   = ccall((:Sleef_expd4_u10, SLEEF_jll.libsleef), FP64P4_d, (FP64P4_d,), v)
    Base.exp(v::FP64P4)     = SIMD.Vec(sleefexp(v.data))
    sleefexp(v::FP64P2_d)   = ccall((:Sleef_expd2_u10, SLEEF_jll.libsleef), FP64P2_d, (FP64P2_d,), v)
    Base.exp(v::FP64P2)     = SIMD.Vec(sleefexp(v.data))
    sleefexp(v::FP32P8_d)   = ccall((:Sleef_expf8_u10, SLEEF_jll.libsleef), FP32P8_d, (FP32P8_d,), v)
    Base.exp(v::FP32P8)     = SIMD.Vec(sleefexp(v.data))
    sleefexp(v::FP32P4_d)   = ccall((:Sleef_expf4_u10, SLEEF_jll.libsleef), FP32P4_d, (FP32P4_d,), v)
    Base.exp(v::FP32P4)     = SIMD.Vec(sleefexp(v.data))

    # log Sleef Vectorized double/single precision natural logarithmic functions with 1.0 ULP error bound
    sleeflog(v::FP64P4_d)   = ccall((:Sleef_logd4_u10, SLEEF_jll.libsleef), FP64P4_d, (FP64P4_d,), v)
    Base.log(v::FP64P4)     = SIMD.Vec(sleeflog(v.data))
    sleeflog(v::FP64P2_d)   = ccall((:Sleef_logd2_u10, SLEEF_jll.libsleef), FP64P2_d, (FP64P2_d,), v)
    Base.log(v::FP64P2)     = SIMD.Vec(sleeflog(v.data))
    sleeflog(v::FP32P8_d)   = ccall((:Sleef_logf8_u10, SLEEF_jll.libsleef), FP32P8_d, (FP32P8_d,), v)
    Base.log(v::FP32P8)     = SIMD.Vec(sleeflog(v.data))
    sleeflog(v::FP32P4_d)   = ccall((:Sleef_logf4_u10, SLEEF_jll.libsleef), FP32P4_d, (FP32P4_d,), v)
    Base.log(v::FP32P4)     = SIMD.Vec(sleeflog(v.data))
    # this should work, but isn't specific enough to override the fallbacks in SIMD.jl
    # Base.log(v::SIMD.Vec{N,T}) where {N,T} = SIMD.Vec(sleeflog(v.data))

    # Vectorized double/single precision base-10 logarithmic functions with 1.0 ULP error bound
    sleeflog10(v::FP64P4_d)   = ccall((:Sleef_log10d4_u10, SLEEF_jll.libsleef), FP64P4_d, (FP64P4_d,), v)
    Base.log10(v::FP64P4)     = SIMD.Vec(sleeflog10(v.data))
    sleeflog10(v::FP64P2_d)   = ccall((:Sleef_log10d2_u10, SLEEF_jll.libsleef), FP64P2_d, (FP64P2_d,), v)
    Base.log10(v::FP64P2)     = SIMD.Vec(sleeflog10(v.data))
    sleeflog10(v::FP32P8_d)   = ccall((:Sleef_log10f8_u10, SLEEF_jll.libsleef), FP32P8_d, (FP32P8_d,), v)
    Base.log10(v::FP32P8)     = SIMD.Vec(sleeflog10(v.data))
    sleeflog10(v::FP32P4_d)   = ccall((:Sleef_log10f4_u10, SLEEF_jll.libsleef), FP32P4_d, (FP32P4_d,), v)
    Base.log10(v::FP32P4)     = SIMD.Vec(sleeflog10(v.data))
else
    # @info "$(@__MODULE__) Using defaults (probably slow scalar fallbacks) for SIMD log, exp, log10 functions"*
    #    " - to enable Sleef library, set USE_SLEEF = true in LocalPreferences.toml and restart your Julia session"
end

end # module
