
using Test

import PALEOboxes as PB


function copy_simd_fillindex!(b, a, ::Type{T}) where T
  
    # N = Val{4}()
    # indvec = [0,0,0,0]
    iter = eachindex(a)

    @inbounds for indices in PB.SIMDutils.SIMDIter(iter, T)
        # println("  indices=", indices)

        xs = PB.SIMDutils.vgatherind(T, a, indices)
        PB.SIMDutils.vscatterind!(xs, b, indices)

    end
    return nothing
end

function copy_scalar_fillindex!(b, a)
  
    iter = eachindex(a)
    ST = eltype(a)

    for i in PB.SIMDutils.SIMDIter(iter, ST)
        # println("  i=", i)

        xs = PB.SIMDutils.vgatherind(ST, a, i)
        PB.SIMDutils.vscatterind!(xs, b, i)

    end
    return nothing
end

function copy_fillindex!(b, a)
  
    iter = eachindex(a)
    
    for i in iter
    
        b[i] = a[i]

    end
    return nothing
end




@testset "scatter gather" begin

a = zeros(7)
a .= 42 .+ (1:length(a))

b = similar(a)
fill!(b, 0.0)
copy_simd_fillindex!(b, a, PB.SIMDutils.FP64P4)
@test b == a

b = similar(a)
fill!(b, 0.0)
copy_scalar_fillindex!(b, a)
@test b == a

b = similar(a)
fill!(b, 0.0)
copy_fillindex!(b, a)
@test b == a

a = zeros(Float32, 27)
a .= 42 .+ (1:length(a))

b = similar(a)
fill!(b, 0.0)
copy_simd_fillindex!(b, a, PB.SIMDutils.FP32P8)
@test b == a

end