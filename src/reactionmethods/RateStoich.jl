import Infiltrator
"""
    RateStoich(
        ratevartemplate, stoich_statevarname;
        deltavarname_eta=nothing, prcessname="", sms_prefix="", sms_suffix="_sms"
    ) -> RateStoich

Calculate fluxes for a biogeochemical reaction given rate, stoichiometry, and optionally isotope eta.

Add to a Reaction using [`create_ratestoich_method`](@ref) and [`add_method_do!`](@ref).  

A Property Variable should be set to provide the reaction rate (often this is implemented by another method of the same Reaction).
This method will then link to that (using the local and link names supplied by `ratevartemplate`)
and calculate the appropriate product rates, omitting products that are not present (`VariableReaction` not linked)
in the `Model` configuration. Metadata for use when analysing model output should be added to the rate variable using [`add_rate_stoichiometry!`](@ref),
in the usual case where this Variable is supplied as `ratevartemplate` this will happen automatically.

# Arguments:
- `ratevartemplate::Union{VarPropT, VarDepT}`: used to define the rate variable local and link names.
- `stoich_statevarname`: collection of Tuple(stoichiometry, name) eg ((-2.0, "O2"), (-1.0,"H2S::Isotope"), (1.0, "SO4::Isotope"))
- `deltavarname_eta`: optional tuple of variable delta + eta ("SO4\\_delta", -30.0) or ("SO4\\_delta", rj.pars.delta). If a Parameter is supplied, this is read in `do_react_ratestoich` to allow modification.
- `processname::String`: optional tag to identify the type of reaction in model output
- `add_rate_stoichiometry=true`: `true` to add call [`add_rate_stoichiometry!`](@ref) to add metadata to `ratevartemplate`.

# Examples:
Create a `RateStoich` representing the reaction   2 O2 + H2S -> H2SO4
 ```jldoctest; setup = :(import PALEOboxes)
julia> myratevar = PALEOboxes.VarProp("myrate", "mol yr-1", "a rate");

julia> rs = PALEOboxes.RateStoich(myratevar, ((-2.0, "O2"), (-1.0,"H2S"), (1.0, "SO4")));

julia> rs.stoich_statevarname
((-2.0, "O2"), (-1.0, "H2S"), (1.0, "SO4"))
```
"""
mutable struct RateStoich
    "template to define rate variable name, units, etc"
    ratevartemplate::Union{VarPropT, VarDepT}
    "label for output analysis (ratevar :rate_processname attribute)"
    processname::String
    "Tuple of (stoich, statevarname) eg ((-2.0, \"O2\"), (-1.0, \"H2S\"), (1.0, \"SO4\"))"
    stoich_statevarname::Tuple
    sms_prefix::String
    sms_suffix::String
    "Tuple (delta_varname, eta)"
    deltavarname_eta
    isotope_data::Type

    "construct a new RateStoich"
    function RateStoich(
        ratevartemplate::Union{VarPropT, VarDepT}, stoich_statevarname;
        processname="",
        sms_prefix="",
        sms_suffix="_sms",
        deltavarname_eta=nothing,
        add_rate_stoichiometry=true,
    )
        rs = new(
            ratevartemplate,
            processname,
            stoich_statevarname,
            sms_prefix,
            sms_suffix,
            deltavarname_eta,
            UndefinedData,
        )

        if add_rate_stoichiometry
            add_rate_stoichiometry!(ratevartemplate, rs)
        end

        return rs
    end

end

"""
    get_stoich_statevarname(reactantnames, prodnames) -> stoich_statevarname

Convert collections of reactant and product names (which may be duplicated) to reaction stoichiometry
in form suitable to create a [`RateStoich`](@ref).

Returns Tuple of (stoich, statevarname) eg ((-2.0, \"O2\"), (-1.0, \"H2S\"), (1.0, \"SO4\"))"
"""
function get_stoich_statevarname(reactantnames, prodnames)
    stoich = Dict{String, Float64}()
    for s in reactantnames
        stoich[s] = get(stoich, s, 0) - 1
    end

    for s in prodnames
        stoich[s] = get(stoich, s, 0) + 1
    end

    stoich_statevarname = Tuple((n, s) for (s, n) in stoich)
    return stoich_statevarname
end


function add_method_do!(reaction::AbstractReaction, ratestoich::RateStoich; kwargs...)
    method = create_ratestoich_method(reaction, ratestoich; kwargs...)

    add_method_do!(reaction, method)

    return method
end

"""
    create_ratestoich_method(reaction::AbstractReaction, ratestoich::RateStoich; isotope_data=ScalarData)
        -> ReactionMethod

Create method (see [`RateStoich`](@ref)).
"""
function create_ratestoich_method(
    @nospecialize(reaction::AbstractReaction), ratestoich::RateStoich;
    isotope_data=ScalarData,
    ignore_unlinked=true,
    operatorID=reaction.operatorID,
    domain=reaction.domain
)
    ratevarlocalname = ratestoich.ratevartemplate.localname # used for diagnostic info
    space = get_attribute(ratestoich.ratevartemplate, :space)

    # if isotopes in use, create _delta variable
    do_isotopes = (isotope_data <: AbstractIsotopeScalar) && !isnothing(ratestoich.deltavarname_eta)
    if do_isotopes
        ratestoich.isotope_data = isotope_data

        delta_varname, tmpeta = ratestoich.deltavarname_eta
        if tmpeta isa Parameter
            etapar = tmpeta
        else
            etapar = ParDouble("eta", tmpeta, units="per mil", description="constant eta")
        end

        delta_var = VarDep(delta_varname, "", "generated by RateStoich.rate="*ratevarlocalname,
                attributes=(:space=>space,))
    else
        ratestoich.isotope_data = ScalarData
        etapar = nothing
        delta_var = nothing
    end
    isotope_dict = Dict("Isotope" => ratestoich.isotope_data)

    # iterate through the flux variables in stoich_statevarname and create if necessary
    vars_sms = VariableReaction[]
    vars_stoich = Float64[]    
    for (stoich, statevarnamefull) in ratestoich.stoich_statevarname
        # parse out ::Isotope, substituting with ratestoich.isotope_data
        statevarname, disotope = split_nameisotope(statevarnamefull, isotope_dict)

        # construct name for flux variable
        smsname = ratestoich.sms_prefix*statevarname*ratestoich.sms_suffix
        smsname = ignore_unlinked ? "("*smsname*")" : smsname

        # create variable        
        push!(
            vars_sms,
            VarContrib(smsname, "", "generated by RateStoich rate="*ratevarlocalname,
                attributes=(:field_data=>disotope, :space=>space,))
        )

        push!(vars_stoich, stoich)
    end

    p = (
            etapar,
            vars_stoich,
            ratestoich,
        )

    return ReactionMethod(
        do_react_ratestoich,
        reaction,
        "RateStoich_"*ratevarlocalname,
        (
            VarList_single(VarDep(ratestoich.ratevartemplate)),
            isnothing(delta_var) ? VarList_nothing() : VarList_single(delta_var),
            VarList_tuple(vars_sms),
        ),
        p,
        operatorID,
        domain
    )
end

"""
    add_rate_stoichiometry!(ratevar::VarPropT, ratestoich::RateStoich)

Add metadata to rate variable `ratevar` for use when analysing model output.

Only needs to be called explicitly if `RateStoich` was supplied with a VarDep that links to the rate variable,
not the rate variable itself.
    
Adds Variable attributes:
- `rate_processname::String = ratestoich.processname`
- `rate_species::Vector{String}` reactants + products from `ratestoich.stoich_statevarname`
- `rate_stoichiometry::Vector{Float64}` reaction stoichiometry from `ratestoich.stoich_statevarname` 
"""
function add_rate_stoichiometry!(
    ratevar, ratestoich::RateStoich,
)
    isa(ratevar, VarPropT) || @warn "add_rate_stoichiometry! ratevar $(fullname(ratevar)) is not a Property Variable, metadata will be ignored"

    set_attribute!(ratevar, :rate_processname, ratestoich.processname; allow_create=true)

    rate_species = String[species for (stoich, species) in ratestoich.stoich_statevarname]
    set_attribute!(ratevar, :rate_species, rate_species; allow_create=true)

    rate_stoichiometry = Float64[stoich for (stoich, species) in ratestoich.stoich_statevarname]
    set_attribute!(ratevar, :rate_stoichiometry, rate_stoichiometry; allow_create=true)

    return nothing
end


"""
    do_react_ratestoich(m::ReactionMethod, vardata, cellrange, deltat)

Calculate rates for reaction products defined by a [`RateStoich`](@ref).
`ratestoich_accessor` should be a NamedTuple generated by
[`create_accessors`](@ref).
"""
function do_react_ratestoich(m::ReactionMethod, (rate, delta, sms_data),  cellrange, deltat)
    (etapar, sms_stoich, _) = m.p

    IteratorUtils.foreach_longtuple_p(_do_sms, sms_data, sms_stoich, (cellrange, etapar, rate, delta))

    return nothing
end

function _do_sms(sms_data, stoich, (cellrange, etapar, rate, delta))
    isotype = eltype(sms_data)
    # @Infiltrator.infiltrate

    if isotype <: AbstractIsotopeScalar
        eta = etapar[]
        @inbounds for idx in cellrange.indices
            sms_data[idx] += isotope_totaldelta(isotype, stoich*rate[idx], delta[idx] + eta)
        end
    else
        @inbounds for idx in cellrange.indices
            sms_data[idx] += stoich*rate[idx]
        end
    end
    return nothing
end
 # ignore unlinked Variables
function _do_sms(sms_data::Nothing, stoich, p)
    return nothing
end
