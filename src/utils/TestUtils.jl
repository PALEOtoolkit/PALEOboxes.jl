"""
    TestUtils

Utility functions for testing and benchmarking
"""

module TestUtils

using BenchmarkTools
using InteractiveUtils

import PALEOboxes as PB

"repeatedly call `dispatch_methodlist` `N` times to generate a long-duration run for profiling"
function profile_dispatchlist(dispatchlist, N, deltat::Float64=0.0)
    for i in 1:N
        PB.dispatch_methodlist(dispatchlist..., deltat)
    end

    return nothing
end


"""
    bench_method(
        dispatchlist, domainname="", reactionname="", methodname="";
        deltat::Float64=0.0,
        use_time=false,
        use_btime=true, 
        do_code_llvm=false,
        do_code_native=true,
        do_code_warntype=false
    )
    
Benchmark or dump llvm code etc for specified method(s).
    
Iterates through `dispatchlist`, if names match then benchmark method else just call method normally.
"""
function bench_method(
    dispatchlist, domainname="", reactionname="", methodname="";
    deltat::Float64=0.0,
    call_all=true,
    use_time=false,
    use_btime=true, 
    do_code_llvm=false,
    do_code_native=false,
    do_code_warntype=false,
)

    # fns, methods, vardatas, crs = dispatchlist
    # iterate through dispatchlist, if names match then benchmark method else just call method
    for (method, vardata, cr) in PB.IteratorUtils.zipstrict(
                                        dispatchlist.methods,
                                        dispatchlist.vardatas, 
                                        dispatchlist.cellranges)
        if (isempty(domainname) || method[].domain.name == domainname) &&
           (isempty(reactionname) || method[].reaction.name == reactionname) &&
           (isempty(methodname) || method[].name == methodname)

            println(PB.fullname(method[]), ":")
            if use_time
                println("    @time:")
                @time PB.call_method(method[], vardata[], cr, deltat)
            end
            if use_btime
                println("    @btime:")
                @btime PB.call_method($method[], $vardata[], $cr, $deltat)
            end
            if do_code_llvm
                println("    @code_llvm:")
                b = IOBuffer()
                PB.call_method_codefn(IOContext(b, :color=>true), code_llvm, method[], vardata[], cr, deltat)
                println(String(take!(b)))                
            end
            if do_code_native
                println("    @code_native:")                
                b = IOBuffer()
                PB.call_method_codefn(IOContext(b, :color=>true), code_native, method[], vardata[], cr, deltat)
                println(String(take!(b)))
            end
            if do_code_warntype
                println("    @code_warntype:")
                b = IOBuffer()
                PB.call_method_codefn(IOContext(b, :color=>true), code_warntype, method[], vardata[], cr, deltat)
                println(String(take!(b)))
                # @code_warntype PB.call_method(method, vardata[], cr, deltat)
            end
        elseif call_all
            PB.call_method(method[], vardata[], cr, deltat)
        end
        
    end

    return nothing
end


"""
    bench_model(model, modeldata; bench_whole=true, domainname="", reactionname="", methodname="", [; kwargs...])
    
Optionally benchmark whole model, and then call [`bench_method`](@ref) to benchmark selected methods

# Keywords
- `bench_whole=true`:  true to benchmark whole model (using `dispatch_methodlist` and `do_deriv`)
- `domainname`, `reactionname`, `methodname`, `kwargs`: passed through to [`bench_method`](@ref)

# Example
Run @code_warntype on a single method:
```julia
julia> PB.TestUtils.bench_model(run.model, modeldata; bench_whole=false, reactionname="atmtransport", methodname="calc_gas_floormrbc",
  use_btime=false, do_code_warntype=true)
````
"""
function bench_model(
    model, modeldata;
    bench_whole=true,
    domainname="",
    reactionname="",
    methodname="",
    kwargs...
)

    dispatchlists = modeldata.dispatchlists_all

    if bench_whole
        println()
        println("@btime  PB.dispatch_list initialize")
        @btime PB.dispatch_methodlist($dispatchlists.list_initialize)

        println()
        println("@btime  PB.dispatch_list do")
        @btime PB.dispatch_methodlist($dispatchlists.list_do)

        println()
        println("@btime do_deriv")
        @btime PB.do_deriv($dispatchlists)
    end

    println()
    println("bench_method initialize")
    PB.TestUtils.bench_method(dispatchlists.list_initialize, domainname, reactionname, methodname; kwargs...) 

    println()
    println("bench_method do")
    PB.TestUtils.bench_method(dispatchlists.list_do, domainname, reactionname, methodname; kwargs...)
    
    return nothing
end


function check_round(a, b, sigdigits; message="")
    lhs = round(a, sigdigits=sigdigits)
    rhs = round(b, sigdigits=sigdigits)
    result = lhs == rhs
    if !result
        println("check_round (", message, ") ",a, " != ",b, " sigdigits=", sigdigits, " round(a) ", lhs, " != round(b) ", rhs)
    end
    return result
end

# macro to validate output
macro check_true(test_expression)
    message = string("check_true failed: ", test_expression)
    :($(esc(test_expression)) || error($message))
end

end
