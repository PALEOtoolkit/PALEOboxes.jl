module IteratorUtils


"""
    named_tuple_ref(keys, eltype) 

Construct a "mutable NamedTuple" with `keys` where `values` are of type Ref{eltype}.

Access values as nt.x[] etc
"""
function named_tuple_ref(keys, eltype)
    return NamedTuple{keys}([Ref{eltype}() for i in 1:length(keys)])
end

"""
    foreach_tuple(f, t1::Tuple)
    foreach_tuple(f, t1::Tuple, t2::Tuple)

Call `f(t1[n])` for each element `n` of `t1::Tuple`.
Call `f(t1[n], t2[n])` for each element `n` of `t1::Tuple`, `t2::Tuple`.

# Implementation
Recursively generates inline code for speed and type stability.

NB: as of Julia 1.6, slows down and allocates if `length(t1) > 32` due to Julia limitation.
See [`foreach_longtuple`](@ref) for this.

See https://github.com/JuliaLang/julia/issues/31869
https://github.com/JuliaLang/julia/blob/master/base/tuple.jl (map implementation)
"""
foreach_tuple(f, t1::Tuple{}) = ()
foreach_tuple(f, t1::Tuple{Any, }) =
    (@Base._inline_meta; f(t1[1]); nothing)
foreach_tuple(f, t1::Tuple) =
    (@Base._inline_meta; f(t1[1]); foreach_tuple(f, Base.tail(t1)); nothing)

foreach_tuple(f, t1::Tuple{}, t2::Tuple{}) = ()
foreach_tuple(f, t1::Tuple{Any, }, t2::Tuple{Any,}) =
    (@Base._inline_meta; f(t1[1], t2[1]); nothing)
foreach_tuple(f, t1::Tuple, t2::Tuple) =
    (@Base._inline_meta; f(t1[1], t2[1]); foreach_tuple(f, Base.tail(t1), Base.tail(t2)); nothing)

"""
    foreach_longtuple(f, t1::Tuple, t2, ... tm)
    foreach_longtuple_p(f, t1::Tuple, t2, ... tm, p)

Call `f(t1[n], t2[n], ... tm[n])` or `f(t1[n], t2[n], ... tm[n], p)` for each element
`n` of `t1::Tuple`, `t2`, ... `tm`.

# Implementation
Uses `@generated` to generate unrolled code for speed and type stability without length restrictions on `tm`.

See https://discourse.julialang.org/t/manually-unroll-operations-with-objects-of-tuple/11604
"""
@generated function foreach_longtuple(f, t1::Tuple)
    ex = quote ; end  # empty expression
    for j=1:fieldcount(t1)
        push!(ex.args, quote; f(t1[$j]); end)
    end
    push!(ex.args, quote; return nothing; end)

    return ex
end

@generated function foreach_longtuple(f, t1::Tuple, t2)
    ex = quote ; end  # empty expression
    for j=1:fieldcount(t1)
        push!(ex.args, quote; f(t1[$j], t2[$j]); end)
    end
    push!(ex.args, quote; return nothing; end)

    return ex
end

@generated function foreach_longtuple(f, t1::Tuple, t2, t3)
    ex = quote ; end  # empty expression
    for j=1:fieldcount(t1)
        push!(ex.args, quote; f(t1[$j], t2[$j], t3[$j]); end)
    end
    push!(ex.args, quote; return nothing; end)

    return ex
end

@generated function foreach_longtuple(f, t1::Tuple, t2, t3, t4)
    ex = quote ; end  # empty expression
    for j=1:fieldcount(t1)
        push!(ex.args, quote; f(t1[$j], t2[$j], t3[$j], t4[$j]); end)
    end
    push!(ex.args, quote; return nothing; end)

    return ex
end

@generated function foreach_longtuple(f, t1::Tuple, t2, t3, t4, t5)
    ex = quote ; end  # empty expression
    for j=1:fieldcount(t1)
        push!(ex.args, quote; f(t1[$j], t2[$j], t3[$j], t4[$j], t5[$j]); end)
    end
    push!(ex.args, quote; return nothing; end)

    return ex
end

@generated function foreach_longtuple_p(f, t1::Tuple, p)
    ex = quote ; end  # empty expression
    for j=1:fieldcount(t1)
        push!(ex.args, quote; f(t1[$j], p); end)
    end
    push!(ex.args, quote; return nothing; end)

    return ex
end

@generated function foreach_longtuple_p(f, t1::Tuple, t2, p)
    ex = quote ; end  # empty expression
    for j=1:fieldcount(t1)
        push!(ex.args, quote; f(t1[$j], t2[$j], p); end)
    end
    push!(ex.args, quote; return nothing; end)

    return ex
end

@generated function foreach_longtuple_p(f, t1::Tuple, t2, t3, p)
    ex = quote ; end  # empty expression
    for j=1:fieldcount(t1)
        push!(ex.args, quote; f(t1[$j], t2[$j], t3[$j], p); end)
    end
    push!(ex.args, quote; return nothing; end)

    return ex
end

@generated function foreach_longtuple_p(f, t1::Tuple, t2, t3, t4, p)
    ex = quote ; end  # empty expression
    for j=1:fieldcount(t1)
        push!(ex.args, quote; f(t1[$j], t2[$j], t3[$j], t4[$j], p); end)
    end
    push!(ex.args, quote; return nothing; end)

    return ex
end

@generated function foreach_longtuple_p(f, t1::Tuple, t2, t3, t4, t5, p)
    ex = quote ; end  # empty expression
    for j=1:fieldcount(t1)
        push!(ex.args, quote; f(t1[$j], t2[$j], t3[$j], t4[$j], t5[$j], p); end)
    end
    push!(ex.args, quote; return nothing; end)

    return ex
end

"""
    reduce_longtuple(f, rinit, t1::Tuple, t2, ... tm) -> r

Call `r += f(r, t1[n], t2[n], ... tm[n])`  for each element
`n` of `t1::Tuple`, `t2`, ... `tm`. Initial value of `r = rinit`

# Implementation
Uses `@generated` to generate unrolled code for speed and type stability without length restrictions on `tm`.

See https://discourse.julialang.org/t/manually-unroll-operations-with-objects-of-tuple/11604
"""

@generated function reduce_longtuple(f, rinit, t1::Tuple)
    ex = quote ; end  # empty expression
    for j=1:fieldcount(t1)
        push!(ex.args, quote; rinit = f(rinit, t1[$j]); end)
    end
    push!(ex.args, quote; return rinit; end)

    return ex
end

@generated function reduce_longtuple(f, rinit, t1::Tuple, t2)
    ex = quote ; end  # empty expression
    for j=1:fieldcount(t1)
        push!(ex.args, quote; rinit = f(rinit, t1[$j], t2[$j]); end)
    end
    push!(ex.args, quote; return rinit; end)

    return ex
end


end # module
