# code from https://github.com/cdsousa/SingleDispatchJuliaVsCXX
# demonstrating Julia dynamic dispatch overhead

using BenchmarkTools

abstract type AbstrType end
struct ConcrType1 <: AbstrType; x::Int; end
struct ConcrType2 <: AbstrType; x::Int; end
@noinline f(a) = a.x

const n = 1_000_000

const arrconcr = [ConcrType1(i) for i=1:n]
const arrabstr = AbstrType[rand(Bool) ? ConcrType1(i) : ConcrType2(i) for i=1:n]
const arrunion = Union{ConcrType1, ConcrType2}[rand(Bool) ? ConcrType1(i) : ConcrType2(i) for i=1:n]

println()

println("concr")
sum_arrconcr = 0::Int
@btime for i=1:n
    @inbounds $sum_arrconcr += f($arrconcr[i])
end

println("abstr")
sum_arrabstr = 0::Int
@btime for i=1:n
    @inbounds $sum_arrabstr += f($arrabstr[i])
end

println("manual dispatch")
# manual dispatch
sum_arrabstr2 = 0::Int
@btime for i=1:n
    @inbounds a = $arrabstr[i]
    T = typeof(a)
    if T === ConcrType1
        $sum_arrabstr2 += f(a::ConcrType1)
    elseif T === ConcrType2
        $sum_arrabstr2 += f(a::ConcrType2)
    else
        $sum_arrabstr2 += f(a)
    end
end

println("union")
sum_arrunion = 0::Int
@btime for i=1:n
    @inbounds $sum_arrunion += f($arrunion[i])
end
