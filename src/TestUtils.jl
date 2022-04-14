"""
    TestUtils

Utility functions for testing and benchmarking
"""

module TestUtils

using BenchmarkTools
using InteractiveUtils

import PALEOboxes as PB



function profile_dispatchlist(dispatchlist, N, deltat::Float64=0.0)
    for i in 1:N
        PB.dispatch_methodlist(dispatchlist..., deltat)
    end

    return nothing
end


"""
    bench_method(dispatchlist, domainname="", reactionname="", methodname="";
            deltat::Float64=0.0,
            use_time=false,
            use_btime=true, 
            code_llvm=false,
            code_native=true,
            code_warntype=false)
    
Benchmark or dump llvm code etc for specified method(s).
    
Iterates through `dispatchlist`, if names match then benchmark method else just call method normally.
"""
function bench_method(dispatchlist, domainname="", reactionname="", methodname="";
    deltat::Float64=0.0,
    call_all=true,
    use_time=false,
    use_btime=true, 
    do_code_llvm=false,
    do_code_native=false,
    do_code_warntype=false)

    # fns, methods, vardatas, crs = dispatchlist
    # iterate through dispatchlist, if names match then benchmark method else just call method
    for (fn, method, vardata, cr) in zip(
                                        dispatchlist.methodfns, 
                                        dispatchlist.methods,
                                        dispatchlist.vardatas, 
                                        dispatchlist.cellranges)
        if (isempty(domainname) || method.domain.name == domainname) &&
           (isempty(reactionname) || method.reaction.name == reactionname) &&
           (isempty(methodname) || method.name == methodname)

            println(PB.fullname(method), ":")
            if use_time
                println("    @time:")
                @time fn(method, vardata[], cr, deltat)
            elseif use_btime
                println("    @btime:")
                @btime $fn($method, $vardata[], $cr, $deltat)
            elseif do_code_llvm
                println("    @code_llvm:")
                @code_llvm fn(method, vardata[], cr, deltat)
            elseif do_code_native
                println("    @code_native:")
                @code_native fn(method, vardata[], cr, deltat)
            elseif do_code_warntype
                println("    @code_warntype:")
                # b = IOBuffer()
                # code_warntype(IOContext(b, :color=>true), 
                #    fn, (typeof(method), typeof(vardata[]), typeof(cr), typeof(deltat))
                #    )
                # println(String(take!(b)))
                @code_warntype fn(method, vardata[], cr, deltat)
            end
        elseif call_all
            fn(method, vardata[], cr, deltat)
        end
        
    end

    return nothing
end



"""
    bench_model(model)
    
Benchmark whole model.

# Example
Run @code_warntype on a single method:
```julia
julia> PB.TestUtils.bench_model(run.model, bench_whole=false, reactionname="atmtransport", methodname="calc_gas_floormrbc",
  use_btime=false, do_code_warntype=true)
````
"""
function bench_model(
    model;
    modeldata=nothing,
    tforce=0.0,
    bench_whole=true,
    domainname="",
    reactionname="",
    methodname="",
    kwargs...
)
    if isnothing(modeldata)
        modeldata =  PB.create_modeldata(model)
        
        PB.allocate_variables!(model, modeldata)
        PB.set_default_solver_view!(model, modeldata)

        PB.check_ready(model, modeldata)

        PB.initialize_reactiondata!(model, modeldata)
            
        PB.dispatch_setup(model, :norm_value, modeldata)
        PB.copy_norm!(modeldata.solver_view_all)
        PB.dispatch_setup(model, :initial_value, modeldata)        
    end

    dispatchlists = modeldata.dispatchlists_all
    PB.set_tforce!(modeldata.solver_view_all, tforce)

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
    
    return dispatchlists
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
