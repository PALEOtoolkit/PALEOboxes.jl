module IteratorUtils

"""
    check_lengths_equal(it1, it2)
    check_lengths_equal(it1, it2, it3)
    check_lengths_equal(it1, it2, it3, it4)
    check_lengths_equal(it1, it2, it3, it4, it5)


Error if iterables it1 ... itn do not have the same length (length(t1) == length(t2) == ...) 
"""
@inline check_lengths_equal(it1) = (length(it1); true)
@inline check_lengths_equal(it1, it2; errmsg="lengths differ") = length(it1) == length(it2) || 
    throw(ArgumentError(errmsg))
@inline check_lengths_equal(it1, it2, it3; errmsg="lengths differ") = length(it1) == length(it2) == length(it3) || 
    throw(ArgumentError(errmsg))
@inline check_lengths_equal(it1, it2, it3, it4; errmsg="lengths differ") = length(it1) == length(it2) == length(it3) == length(it4) || 
    throw(ArgumentError(errmsg))
@inline check_lengths_equal(it1, it2, it3, it4, it5; errmsg="lengths differ") = length(it1) == length(it2) == length(it3) == length(it4) == length(it5) || 
    throw(ArgumentError(errmsg))


"""
    zipstrict(iters...; errmsg="iterables lengths differ")

`Base.zip` with additional check that lengths of iterables are equal.   
"""
@inline zipstrict(iters...; errmsg="iterables lengths differ") = (check_lengths_equal(iters...; errmsg=errmsg); zip(iters...))


"""
    named_tuple_ref(keys, eltype) 

Construct a "mutable NamedTuple" with `keys` where `values` are of type Ref{eltype}.

Access values as nt.x[] etc
"""
function named_tuple_ref(keys, eltype)
    return NamedTuple{keys}([Ref{eltype}() for i in 1:length(keys)])
end

"""
    foreach_tuple(f, t1::Tuple; errmsg="iterables lengths differ")
    foreach_tuple(f, t1::Tuple, t2::Tuple; errmsg="iterables lengths differ")

Call `f(t1[n])` for each element `n` of `t1::Tuple`.
Call `f(t1[n], t2[n])` for each element `n` of `t1::Tuple`, `t2::Tuple`.

# Implementation
Recursively generates inline code for speed and type stability.

NB: as of Julia 1.6, slows down and allocates if `length(t1) > 32` due to Julia limitation.
See [`foreach_longtuple`](@ref) for this.

See https://github.com/JuliaLang/julia/issues/31869
https://github.com/JuliaLang/julia/blob/master/base/tuple.jl (map implementation)
"""
@inline foreach_tuple(f::F, tuples...; errmsg="iterables lengths differ") where{F} = 
    (check_lengths_equal(tuples...; errmsg=errmsg); foreach_tuple_unchecked(f, tuples...))

foreach_tuple_unchecked(f, t1::Tuple{}) = ()
foreach_tuple_unchecked(f, t1::Tuple{Any, }) =
    (@Base._inline_meta; f(t1[1]); nothing)
foreach_tuple_unchecked(f, t1::Tuple) =
    (@Base._inline_meta; f(t1[1]); foreach_tuple_unchecked(f, Base.tail(t1)); nothing)

foreach_tuple_unchecked(f, t1::Tuple{}, t2::Tuple{}) = ()
foreach_tuple_unchecked(f, t1::Tuple{Any, }, t2::Tuple{Any,}) =
    (@Base._inline_meta; f(t1[1], t2[1]); nothing)
foreach_tuple_unchecked(f, t1::Tuple, t2::Tuple) =
    (@Base._inline_meta; f(t1[1], t2[1]); foreach_tuple_unchecked(f, Base.tail(t1), Base.tail(t2)); nothing)

"""
    foreach_longtuple(f, t1::Tuple, t2, ... tm; errmsg="iterables lengths differ")
    foreach_longtuple_p(f, t1::Tuple, t2, ... tm, p; errmsg="iterables lengths differ")

Call `f(t1[n], t2[n], ... tm[n])` or `f(t1[n], t2[n], ... tm[n], p)` for each element
`n` of `t1::Tuple`, `t2`, ... `tm`.

# Implementation
Uses `@generated` to generate unrolled code for speed and type stability without length restrictions on `tm`.

See https://discourse.julialang.org/t/manually-unroll-operations-with-objects-of-tuple/11604
"""
@inline foreach_longtuple(f::F, tuples...; errmsg="iterables lengths differ") where{F} = 
    (check_lengths_equal(tuples...; errmsg=errmsg); foreach_longtuple_unchecked(f, tuples...))

@generated function foreach_longtuple_unchecked(f, t1::Tuple)
    ex = quote ; end  # empty expression
    for j=1:fieldcount(t1)
        push!(ex.args, quote; f(t1[$j]); end)
    end
    push!(ex.args, quote; return nothing; end)

    return ex
end

@generated function foreach_longtuple_unchecked(f, t1::Tuple, t2)
    ex = quote ; end  # empty expression
    for j=1:fieldcount(t1)
        push!(ex.args, quote; f(t1[$j], t2[$j]); end)
    end
    push!(ex.args, quote; return nothing; end)

    return ex
end

@generated function foreach_longtuple_unchecked(f, t1::Tuple, t2, t3)
    ex = quote ; end  # empty expression
    for j=1:fieldcount(t1)
        push!(ex.args, quote; f(t1[$j], t2[$j], t3[$j]); end)
    end
    push!(ex.args, quote; return nothing; end)

    return ex
end

@generated function foreach_longtuple_unchecked(f, t1::Tuple, t2, t3, t4)
    ex = quote ; end  # empty expression
    for j=1:fieldcount(t1)
        push!(ex.args, quote; f(t1[$j], t2[$j], t3[$j], t4[$j]); end)
    end
    push!(ex.args, quote; return nothing; end)

    return ex
end

@generated function foreach_longtuple_unchecked(f, t1::Tuple, t2, t3, t4, t5)
    ex = quote ; end  # empty expression
    for j=1:fieldcount(t1)
        push!(ex.args, quote; f(t1[$j], t2[$j], t3[$j], t4[$j], t5[$j]); end)
    end
    push!(ex.args, quote; return nothing; end)

    return ex
end

@inline foreach_longtuple_p(f::F, t1, p; errmsg="iterables lengths differ") where{F} =
    foreach_longtuple_unchecked_p(f, t1, p)
@generated function foreach_longtuple_unchecked_p(f, t1::Tuple, p)
    ex = quote ; end  # empty expression
    for j=1:fieldcount(t1)
        push!(ex.args, quote; f(t1[$j], p); end)
    end
    push!(ex.args, quote; return nothing; end)

    return ex
end

@inline foreach_longtuple_p(f::F, t1, t2, p; errmsg="iterables lengths differ") where{F} =
    foreach_longtuple_unchecked_p(f, t1, t2, p)
@generated function foreach_longtuple_unchecked_p(f, t1::Tuple, t2, p)
    ex = quote ; end  # empty expression
    for j=1:fieldcount(t1)
        push!(ex.args, quote; f(t1[$j], t2[$j], p); end)
    end
    push!(ex.args, quote; return nothing; end)

    return ex
end

@inline foreach_longtuple_p(f::F, t1, t2, t3, p; errmsg="iterables lengths differ") where{F} =
    foreach_longtuple_unchecked_p(f, t1, t2, t3, p)
@generated function foreach_longtuple_unchecked_p(f, t1::Tuple, t2, t3, p)
    ex = quote ; end  # empty expression
    for j=1:fieldcount(t1)
        push!(ex.args, quote; f(t1[$j], t2[$j], t3[$j], p); end)
    end
    push!(ex.args, quote; return nothing; end)

    return ex
end

@inline foreach_longtuple_p(f::F, t1, t2, t3, t4, p; errmsg="iterables lengths differ") where{F} =
    foreach_longtuple_unchecked_p(f, t1, t2, t3, t4, p)
@generated function foreach_longtuple_unchecked_p(f, t1::Tuple, t2, t3, t4, p)
    ex = quote ; end  # empty expression
    for j=1:fieldcount(t1)
        push!(ex.args, quote; f(t1[$j], t2[$j], t3[$j], t4[$j], p); end)
    end
    push!(ex.args, quote; return nothing; end)

    return ex
end

@inline foreach_longtuple_p(f::F, t1, t2, t3, t4, t5, p; errmsg="iterables lengths differ") where{F} =
    foreach_longtuple_unchecked_p(f, t1, t2, t3, t4, t5, p)
@generated function foreach_longtuple_unchecked_p(f, t1::Tuple, t2, t3, t4, t5, p)
    ex = quote ; end  # empty expression
    for j=1:fieldcount(t1)
        push!(ex.args, quote; f(t1[$j], t2[$j], t3[$j], t4[$j], t5[$j], p); end)
    end
    push!(ex.args, quote; return nothing; end)

    return ex
end

"""
    reduce_longtuple(f, rinit, t1::Tuple, t2, ... tm; errmsg="iterables lengths differ") -> r

Call `r += f(r, t1[n], t2[n], ... tm[n])`  for each element
`n` of `t1::Tuple`, `t2`, ... `tm`. Initial value of `r = rinit`

# Implementation
Uses `@generated` to generate unrolled code for speed and type stability without length restrictions on `tm`.

See https://discourse.julialang.org/t/manually-unroll-operations-with-objects-of-tuple/11604
"""
@inline reduce_longtuple(f::F, rinit, tuples...; errmsg="iterables lengths differ") where{F} = 
    (check_lengths_equal(tuples...; errmsg=errmsg); reduce_longtuple_unchecked(f, rinit, tuples...))

@generated function reduce_longtuple_unchecked(f, rinit, t1::Tuple)
    ex = quote ; end  # empty expression
    for j=1:fieldcount(t1)
        push!(ex.args, quote; rinit = f(rinit, t1[$j]); end)
    end
    push!(ex.args, quote; return rinit; end)

    return ex
end

@generated function reduce_longtuple_unchecked(f, rinit, t1::Tuple, t2)
    ex = quote ; end  # empty expression
    for j=1:fieldcount(t1)
        push!(ex.args, quote; rinit = f(rinit, t1[$j], t2[$j]); end)
    end
    push!(ex.args, quote; return rinit; end)

    return ex
end


end # module
