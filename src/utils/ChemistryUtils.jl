module ChemistryUtils

import PALEOboxes as PB
import OrderedCollections

"""
    parse_chemical_formula(formula::AbstractString) -> element_counts::OrderedDict{Symbol, Float64}

Parse a chemical formula into a Dict of element counts

# Examples
```jldoctest; setup = :(import PALEOboxes as PB)
julia> PB.ChemistryUtils.parse_chemical_formula("HO2NO2")
OrderedCollections.OrderedDict{Symbol, Float64} with 3 entries:
  :H => 1.0
  :N => 1.0
  :O => 4.0
```
"""
function parse_chemical_formula(formula::AbstractString)
    elementcounts = OrderedCollections.OrderedDict{Symbol, Float64}()
    for m in eachmatch(r"([A-Z][a-z]*)(\d*)", formula)
        elementstr, countstr = m.captures
        element = Symbol(elementstr)
        count = isempty(countstr) ? 1 : parse(Int, countstr)
        elementcounts[element] = get(elementcounts, element, 0) + count
    end

    sort!(elementcounts)

    return elementcounts
end



"""
    calc_molar_mass(element_counts::AbstractDict) -> molar_mass::Float64
    calc_molar_mass(element_counts::AbstractVector{<:Pair}) -> molar_mass::Float64
    calc_molar_mass(element_counts::NamedTuple) -> molar_mass::Float64

Add up element masses to get molar mass (g)

# Examples
```jldoctest; setup = :(import PALEOboxes as PB)
julia> PB.ChemistryUtils.calc_molar_mass([:C=>1, :O=>2])
44.009
```
"""
function calc_molar_mass(element_counts)

    molar_mass = 0.0
    for (element, count) in element_counts
        element isa Symbol || error("element $element is not a Symbol eg :O, :C")
        haskey(PB.Constants.STANDARD_ATOMIC_WEIGHTS, element) || error("unknown element $element")
        count isa Number || error("element count $count is not a Number")
        molar_mass += count*PB.Constants.STANDARD_ATOMIC_WEIGHTS[element]
    end

    return molar_mass
end

calc_molar_mass(element_counts::NamedTuple) = calc_molar_mass(pairs(element_counts))

end # module
