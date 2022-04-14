
# Tests for Tuple iteration and dispatch 
# see https://github.com/JuliaLang/julia/issues/31869
# https://github.com/JuliaLang/julia/blob/master/base/tuple.jl (map implementation)
#
# (looking for a simple way to dispatch a Tuple of similar Reactions (do_ in Model.jl) 
# or Variables (RateStoich.jl, several other places) that uses the Julia type
# system to eliminate dynamic dispatch overhead, both speed and spurious allocations)


using BenchmarkTools


foreach_test(f, t::Tuple{})     = ()
# foreach_test(f, t::Tuple{Any, }) = (f(t[1]); nothing)
# foreach_test(f, t::Tuple)    = (f(t[1]); foreach_test(f, Base.tail(t)); nothing)
foreach_test(f, t::Tuple{Any, }) = (@Base._inline_meta; f(t[1]); nothing)
foreach_test(f, t::Tuple)    = (@Base._inline_meta; f(t[1]); foreach_test(f, Base.tail(t)); nothing)

# benchmarks demonstrate Base.tail hits some arbitrary limit at 32 ??

t2 = Tuple(ones(32))
@btime foreach_test(sin, $t2) # no allocations

t2 = Tuple(ones(33))
@btime foreach_test(sin, $t2) # allocates

fprint(a) = println("a=$a typeof(a)=$(typeof(a))")

foreach_test2(f, t1::Tuple{}, t2::Tuple{})     = ()
# foreach_test(f, t::Tuple{Any, }) = (f(t[1]); nothing)
# foreach_test(f, t::Tuple)    = (f(t[1]); foreach_test(f, Base.tail(t)); nothing)
foreach_test2(f, t1::Tuple{Any, }, t2::Tuple{Any,}) = (@Base._inline_meta; f(t1[1], t2[1]); nothing)
foreach_test2(f, t1::Tuple, t2::Tuple)    = (@Base._inline_meta; f(t1[1], t2[1]); foreach_test2(f, Base.tail(t1), Base.tail(t2)); nothing)

t2 = Tuple(ones(32))
@btime foreach_test2(min, $t2, $t2 ) # no allocations

t2 = Tuple(ones(33))
@btime foreach_test2(min, $t2, $t2) # allocates