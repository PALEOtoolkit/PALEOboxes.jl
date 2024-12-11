"""
    value_ad(x::ADT) -> x

Get scalar value from variable `x` (discarding any AD derivatives).

This can be used to exclude `x` from automatic differentation and hence a Jacobian calculation,
eg where `x` is a small contribution but would make the Jacobian much denser.

Model code should implement this for any AD types used, eg

    value_ad(x::SparsityTracing.ADval) = SparsityTracing.value(x)
    value_ad(x::ForwardDiff.Dual) = ForwardDiff.value(x)

"""
value_ad(x) = x

# get scalar or ad from variable x, as specified by first argument
value_ad(::Type{T}, x::T) where {T} = x        # pass through AD
value_ad(::Type{T}, x::Float64) where {T} = x  # pass through Float64
value_ad(::Type{Float64}, x)  = value_ad(x)    # strip AD
value_ad(::Type{Float64}, x::Float64)  = x     # avoid method ambiguity

"""
    zero_ad(x...) -> 0.0*x

Provide a zero of type of `x` (or type of `x1*x2*...` if given multiple arguments), retaining AD dependency information.

Workaround to enable use of conditional logic whilst retaining dependency information for tracing Jacobian sparsity pattern.
"""
zero_ad(x) = 0.0*x
zero_ad(x1, x2) = 0.0*x1*x2
zero_ad(x1, x2, x3) = 0.0*x1*x2*x3
zero_ad(x1, x2, x3, x4) = 0.0*x1*x2*x3*x4

"""
    smoothstepcubic(x, xedge, xwidth) -> y

Smoothed step function over width `xwidth` at location `xedge`.

Provides a way of smoothing a discontinuous step function so that the
derivative exists, to help numerical solution methods that require a
Jacobian (eg stiff ODE solvers).

Uses a cubic function in range `xedge +/- xwidth/2`, so first derivative
is continuous, higher derivatives are not.

Uses [`zero_ad`](@ref) to retain AD dependency information for tracing Jacobian sparsity pattern.

Returns:
- 0.0 for x < (xedge - xwidth/2)
- 1.0 for x > (xedge + xwidth/2)
- a smoothed step for xedge-xwidth/2 < x < xedge+xwidth/2
"""
function smoothstepcubic(x, xedge, xwidth)
    # rescale to 0 < xs < 1
    xs = (x - xedge + 0.5*xwidth)/xwidth
    # xs smoothly steps from 0 to 1 over interval 0 < xs < 1
    if xs > 1.0
        return one(x) + zero_ad(x)  
    elseif xs < 0.0
        return zero_ad(x)
    else
        return xs*xs*(3 - 2 * xs)
    end
end