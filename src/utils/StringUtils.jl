"""
    parse_number_name(nname::AbstractString) -> (number, name)

Parse a string of form "-1*A"

"""
function parse_number_name(
    nname::AbstractString; 
    sep=['*', ' '], 
    number_first=true,
    io=IOBuffer(),
    errmsg="invalid field in nname, not of form number*name: ",
)

    # parse multiplier
    svmn = split(nname, sep, keepempty=false)
    mult = nothing
    if length(svmn) == 1
        mult, name = 1, svmn[1]
    elseif length(svmn) == 2
        if !number_first
            tmp = svmn[1]
            svmn[1] = svmn[2]
            svmn[2] = tmp
        end
        mult = tryparse(Int64, svmn[1])
        if isnothing(mult)
            mult = tryparse(Float64, svmn[1])
        end
        name = svmn[2]
    end

    !isnothing(mult) || infoerror(io, errmsg*nname)
  
    return (mult, name)
end

parse_name_to_power_number(nname::AbstractString) = parse_number_name(nname; sep=['^', ' '], number_first=false)
