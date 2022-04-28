"""
    Fluxes

Create and manage collections of Target and Contrib Variables used to define fluxes within or between Domains.

# Conventions for defining biogeochemical fluxes within Domains

- Particulate organic matter (including CaCO3) with stoichiometry Corg:N:P:Ccarb is usually transferred as:
  `flux_list` `["P", "N", "Corg::CIsotope", "Ccarb::CIsotope"]`
  with a prefix indicating the function (eg `export_`).
- Solute fluxes are usually transferred as:
  `flux_list`  `["DIC::CIsotope", "TAlk", "Ca", "P", "O2", "SO4::SIsotope", "H2S::SIsotope", "CH4::CIsotope"]`
  (where these names match the Reservoir names for the solutes), with a prefix indicating the function
  (usually just `flux_` or `soluteflux_`).

# Conventions for defining global fluxes between modules 

The model configuration `.yaml` file should create 
a `Domain` for each global flux, containing one or more [`Fluxes.ReactionFluxTarget`](@ref).  Fluxes are then transferred
(copied) by adding a [`Fluxes.ReactionFluxTransfer`](@ref) to each destination Domain.

Naming conventions for Earth system fluxes:

|Domain name            | target prefix | flux_list (illustrative, add as needed) |
|:----------------------|:---------------|:--------------------|
|fluxAtoLand            | flux_          | ["CO2::CIsotope", "O2"] |
|fluxRtoOcean           | flux_          | ["DIC::CIsotope", "TAlk", "Ca", "P", "SO4::SIsotope"]|
|fluxOceanBurial        | flux_          | ["Corg::CIsotope", "Ccarb::CIsotope", "Porg", "Pauth", "PFe", "P", "GYP::SIsotope", "PYR::SIsotope"]|
|fluxSedCrusttoAOcean   | flux_          | ["C::CIsotope", "S::SIsotope", "Redox"]|
|fluxLandtoSedCrust     | flux_          |["Ccarb::CIsotope", "Corg::CIsotope", "GYP::SIsotope", "PYR::SIsotope"]|
|fluxOceanfloor         | particulateflux_|["P", "N", "Corg::CIsotope", "Ccarb::CIsotope"] |
|                       | soluteflux_   | ["DIC::CIsotope", "TAlk", "Ca", "P", "O2", "SO4::SIsotope", "H2S::SIsotope", "CH4::CIsotope"] |
|fluxAtmtoOceansurface  | flux_         | ["CO2::CIsotope", "CH4::CIsotope", "O2"] |

"""
module Fluxes

import PALEOboxes as PB
import SparseArrays
import LinearAlgebra # for I

import Infiltrator # Julia debugger

using ..DocStrings
"""
    FluxContrib(
        fluxprefix::AbstractString, flux_list;
        isotope_data::Dict,
        [, space=PB.CellSpace][, alloptional=true]
    ) -> vartuple::NamedTuple

Create a `NamedTuple` of `VarContrib` from `flux_prefix.*flux_list`, needed for a Reaction to write to a flux coupler.
Entries in `flux_list` of form `fluxname::isotope` are parsed to look up isotope type from `isotope_data` Dict.

For each `fluxname` in `flux_list`, the generated `VarContrib` `vartuple.fluxname` have
`link_namestr=\$(fluxprefix)\$(fluxname)"` and 
`localname="\$(fluxprefix)\$(fluxname)` (with any `.` in `fluxprefix` substituted to `_`).

# Example:
```julia
julia> fluxRtoOcean = PB.Fluxes.FluxContrib("fluxRtoOcean.flux_", ["P", "SO4::SIsotope"], isotope_data=Dict("SIsotope"=>PB.IsotopeLinear), space=PB.ScalarSpace);

julia> fluxRtoOcean.P
PALEOboxes.VariableReaction{PALEOboxes.VT_ReactContributor}
  localname='fluxRtoOcean_flux_P'
  link_name='fluxRtoOcean.flux_P[unlinked]'
  attributes=Dict{Symbol, Any}(:description => "flux P", :space => PALEOboxes.ScalarSpace, :standard_variable => false, :norm_value => 1.0, :data_dims => String[], :long_name => "", :vfunction => PALEOboxes.VF_Undefined, :field_data => PALEOboxes.ScalarData, :initial_delta => 0.0, :units => "mol yr-1"…)

julia> fluxRtoOcean.SO4
PALEOboxes.VariableReaction{PALEOboxes.VT_ReactContributor}
  localname='fluxRtoOcean_flux_SO4'
  link_name='fluxRtoOcean.flux_SO4[unlinked]'
  attributes=Dict{Symbol, Any}(:description => "flux SO4", :space => PALEOboxes.ScalarSpace, :standard_variable => false, :norm_value => 1.0, :data_dims => String[], :long_name => "", :vfunction => PALEOboxes.VF_Undefined, :field_data => PALEOboxes.IsotopeLinear, :initial_delta => 0.0, :units => "mol yr-1"…)
```
"""
FluxContrib(fluxprefix::AbstractString, flux_list; alloptional=true, kwargs...) = 
    _FluxVars(PB.VarContrib, fluxprefix, flux_list; alloptional=alloptional, kwargs...)
FluxContribScalar(fluxprefix::AbstractString, flux_list; kwargs...) = 
    _FluxVars(PB.VarContrib, fluxprefix, flux_list; space=PB.ScalarSpace, kwargs...)

"""
    FluxTarget(fluxprefix::AbstractString, flux_list;
        isotope_data::Dict, 
        [, space=PB.CellSpace]
    ) -> vartuple::NamedTuple

Create a `NamedTuple` of `VarTarget` from `flux_list`, needed for a Reaction to define a flux coupler.
See [`FluxContrib`](@ref) for arguments.
"""
FluxTarget(fluxprefix::AbstractString, flux_list; kwargs...) = 
        _FluxVars(PB.VarTarget, fluxprefix, flux_list; kwargs...)
FluxTargetScalar(fluxprefix::AbstractString, flux_list; kwargs...) = 
        _FluxVars(PB.VarTarget, fluxprefix, flux_list; space=PB.ScalarSpace, kwargs...)
    
"Create a NamedTuple of Variables using `varctorfn`"
function _FluxVars(
    varctorfn, fluxprefix::AbstractString, flux_list; 
    isotope_data::Dict,
    localname_prefix=nothing,
    alloptional=false,
    description="flux",
    space=PB.CellSpace,
)

    (fluxdomain, fluxsubdomain, fluxnameprefix) = PB.split_link_name(fluxprefix)

    if isnothing(localname_prefix)
        localname_prefix = PB.combine_link_name(fluxdomain, fluxsubdomain, fluxnameprefix, sep="_")        
    end

    field_names = []
    variables = []

    for fluxnameisotopestr in flux_list
        optional, fluxnameisotope = PB.strip_brackets(fluxnameisotopestr)
        fluxname, IsotopeType = PB.split_nameisotope(fluxnameisotope, isotope_data)
        push!(field_names, Symbol(fluxname))
        attributes = (:space=>space, :field_data=>IsotopeType) 

        localname = localname_prefix*fluxname

        link_namestr = PB.combine_link_name(fluxdomain, fluxsubdomain, fluxnameprefix*fluxname)
        if alloptional || optional
            link_namestr = "("*link_namestr*")"
        end

        push!(
            variables, 
            varctorfn(
                localname,  "mol yr-1", "$description $fluxname", 
                link_namestr=link_namestr, 
                attributes=attributes
            )
        )
    end

    return NamedTuple{Tuple(field_names)}(variables)
end



"""
    ReactionFluxTarget

Provides either a target for fluxes from `fluxlist` in an input `Domain` or a constant stub, optionally calculates totals.

For each `fluxname` in `fluxlist`, creates:
- if `const_stub==false`, an input `VarTarget` `pars.target_prefix.v*"flux_"*fluxname`
  (or if `const_stub==true`, a constant `VarProp`).
- if `flux_totals` is `true`, a total `VarPropScalar` `target_prefix.v*"flux_total_"*fluxname`.

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionFluxTarget{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParStringVec("fluxlist", ["example"], 
            description="available fluxes"),

        PB.ParString("target_prefix", "flux_",
            description="target names will be \"<target_prefix><fluxname>\",             
                target total names will be \"<target_prefix>total_<fluxname>\"
                (where <fluxname> is an entry in fluxlist)"),

        PB.ParBool("flux_totals",     false, 
            description="true to calculate flux totals (as \"<target_prefix>flux_total_<fluxname>\")"),

        PB.ParBool("const_stub", false,
            description="true to provide constant flux defined by :initial_value, :initial_delta attributes "*
                        "on input Target Variables \"<target_prefix>flux_<fluxname>\""),
    )
end


function PB.register_methods!(rj::ReactionFluxTarget)
    # Add flux targets and transfers
   
    target_variables = []
    total_localnames = []

    target_prefix = rj.pars.target_prefix.v
    if !isempty(target_prefix) && target_prefix[end] != '_'
        target_prefix = target_prefix*"_"
        @warn "Reaction $(PB.fullname(rj)) bodging target_prefix without trailing underscore $(rj.pars.target_prefix.v) -> $target_prefix"
    end
 
    for fluxnameisotope in rj.pars.fluxlist.v
        fluxname, IsotopeType = PB.split_nameisotope(fluxnameisotope, rj.external_parameters)

        targetname = target_prefix*fluxname
        
        if rj.pars.const_stub.v
            target_var = PB.VarPropStateIndep(targetname, "mol yr-1", "constant flux", attributes=(:field_data=>IsotopeType,))
        else
            target_var = PB.VarTarget(targetname, "mol yr-1", "flux input", attributes=(:field_data=>IsotopeType,))
        end
        if rj.pars.flux_totals.v
            PB.set_attribute!(target_var, :calc_total, true, allow_create=true)
        end

        push!(target_variables, target_var) 

        push!(total_localnames, target_prefix*"total_"*fluxname)
  
    end
        
    if rj.pars.const_stub.v
        # method to set constant value
        # supply filterfn to initialize supplied Variables even though they aren't state variables
        PB.add_method_setup_initialvalue_vars_default!(rj, target_variables, filterfn = v->true) 
    else
        # target_variables not used by us, but must appear in a method to be linked and created
        PB.add_method_do_nothing!(rj, target_variables)
    end

    if rj.pars.flux_totals.v
        totals_method = PB.create_totals_method(rj, target_variables, total_localnames=total_localnames)        
        PB.add_method_do!(rj, totals_method)
    end

    PB.add_method_initialize_zero_vars_default!(rj)

    return nothing
end




"""
    ReactionFluxTransfer

Copy fluxes, optionally using transfer matrices to define cell-cell mappings if input and output fluxes 
are in a different Domain.

There are three common cases:
- Transfer within a Domain, using `transfer_matrix` `Identity`
- Transfer from a `fluxXXX` Domain to the Domain hosting the `ReactionFluxTransfer`, using `transfer_matrix` to
  define a mapping if these Domains are of different sizes.  
- Transfer from a `fluxXXX` Domain to the boundary cells of an interior Domain (eg `ocean.oceansurface`), 
  where the Domain hosting the `ReactionFluxTransfer` is the corresponding boundary Domain (eg `oceansurface`).

# Parameters
$(PARS)
"""
Base.@kwdef mutable struct ReactionFluxTransfer{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParStringVec("fluxlist", String[], 
            description="available fluxes"),

        PB.ParString("input_fluxes", "[inputdomain.][inputsubdomain.]flux_\$fluxname\$",
            description="string to match to find input flux Variables. "*
                "These Variables will local names \"input_\$fluxname\$\"."),
        PB.ParString("output_fluxes", "[outputdomain.][outputsubdomain.]\$fluxname\$_sms", 
            description="string to use to generate output flux Variables where \$fluxname\$ is substituted from input_fluxes. "*
                "These Variables will local names \"output_\$fluxname\$\", and are optional (ie ignored if not linked). "*
                "Usually the output Domain is the Domain hosting the ReactionFluxTransfer (in general, it must be a Domain "*
                "of the same size as the Domain hosting the ReactionFluxTransfer)."),

        PB.ParDouble("transfer_multiplier", 1.0, 
            description="scalar multiplier for transfer"),
        PB.ParString("transfer_matrix", "Identity", 
            allowed_values=["Identity", "Distribute", "Custom"],
            description="matrix defining input Domain (length n) -> output Domain (length m) mapping: 
                'Identity' copies fluxes directly, requires m == n; 
                'Distribute' uniformly distributes fluxes from input->output: sums over n and them distributes fraction 1/m evenly to each m;
                'Custom' requires transfer matrix to be supplied"),
    )

    input_domain = nothing

    transfer_matrix_tr::SparseArrays.SparseMatrixCSC{Float64, Int64} = SparseArrays.spzeros(0, 0)
end

PB.register_methods!(rj::ReactionFluxTransfer) = nothing

function PB.register_dynamic_methods!(rj::ReactionFluxTransfer, model::PB.Model)
    
    (input_domain_name, input_subdomain_name, input_var_matchname) =
        PB.split_link_name(rj.pars.input_fluxes.v)

    # get Domain for input fluxes 
    rj.input_domain = isempty(input_domain_name) ? rj.domain : PB.get_domain(model, input_domain_name)   
    !isnothing(rj.input_domain) ||
        error("$(PB.fullname(rj)) no Domain $input_domain_name") 

    input_var_matchname_split = split(input_var_matchname, "\$fluxname\$")
    length(input_var_matchname_split) == 2 ||
        error("$(PB.fullname(rj)) invalid input_fluxes: must contain \$fluxname\$ $(rj.pars.input_fluxes.v)")

    # Regular expression so $fluxname$ matches any alphanumeric characters    
    input_var_regex = input_var_matchname_split[1]*r"(\p{Xan}+)"*input_var_matchname_split[2]

    (output_domain_name, output_subdomain_name, output_var_patternname) =
        PB.split_link_name(rj.pars.output_fluxes.v)
    occursin("\$fluxname\$", output_var_patternname) ||
        error("$(PB.fullname(rj)) invalid output_fluxes: must contain \$fluxname\$ $(rj.pars.output_fluxes.v)")

    flux_names, var_inputs, var_outputs = [], [], []

    @debug "register_dynamic_methods! $(PB.fullname(rj)):"
    for input_var in PB.get_variables(rj.input_domain)
        m = match(input_var_regex, input_var.name)
        #   bodge in additional check as regex seems wrong when $fluxname is at end of string, 
        #   "prefix_$fluxname" apparently will also match "prefix_$fluxname_otherstuff"
        if !isnothing(m) && (input_var.name == input_var_matchname_split[1]*m.captures[1]*input_var_matchname_split[2])
            fluxname = m.captures[1]
            push!(flux_names, fluxname)

            input_namestr = PB.combine_link_name(
                input_domain_name, 
                input_subdomain_name,
                input_var.name
            )
            input_attrb = Tuple(atnm=>PB.get_attribute(input_var, atnm) for atnm in (:field_data, :data_dims, :space))
           
            push!(var_inputs, 
                PB.VarDep("input_"*fluxname, "mol yr-1", "flux input",
                    link_namestr=input_namestr,                  
                    attributes=input_attrb,
                )
            )

            output_namestr = "("*PB.combine_link_name(
                output_domain_name, 
                output_subdomain_name, 
                replace(output_var_patternname, "\$fluxname\$" => fluxname)
            )*")"
            push!(var_outputs, 
                PB.VarContrib("output_"*fluxname, "mol yr-1", "flux output",
                    link_namestr=output_namestr,
                    attributes=input_attrb,
                )
            )
        else
            @debug "  ignoring $(PB.fullname(input_var))"
        end
    end
     
    PB.add_method_do!(
        rj, 
        do_transfer,
        (PB.VarList_components(var_inputs), PB.VarList_components(var_outputs, allow_unlinked=true)),
        p = flux_names,
        preparefn = prepare_do_transfer
    )
   
    return nothing
end

"remove unlinked Variables and create transfer_matrix"
function prepare_do_transfer(m::PB.ReactionMethod, (input_vardata, output_vardata))

    @info "prepare_do_transfer: ReactionMethod $(PB.fullname(m))"

    # write diagnostic output
    @info "  active fluxes:"
    var_inputs, var_outputs = PB.get_variables.(m.vars)
    flux_names = m.p
    nactive = 0
    for (var_input, var_output, fluxname) in zip(var_inputs, var_outputs, flux_names)
        if isnothing(var_output.linkvar)
            @info "    $(rpad(fluxname, 20))   $(rpad(PB.fullname(var_input.linkvar),40)) output not linked"
        else
            @info "    $(rpad(fluxname, 20))   $(rpad(PB.fullname(var_input.linkvar),40)) -> $(PB.fullname(var_output.linkvar))"
            nactive += 1
        end
    end

    # remove unlinked Variables
    input_vardata_active = []
    output_vardata_active = []
    for (ainput, aoutput) in zip(input_vardata, output_vardata)
        if !isnothing(aoutput)
            push!(input_vardata_active, ainput)
            push!(output_vardata_active, aoutput)
        end
    end

    if !isempty(input_vardata_active)
        # figure out length of input and output arrays
        input_length = 0
        output_length = 0
        for (ainput, aoutput) in zip(input_vardata_active, output_vardata_active)
            
            if iszero(input_length)
                input_length = prod(size(ainput))::Int
            end                
            input_length == prod(size(ainput)) || 
                error("$(PB.fullname(m)) input Variables are not all the same size")
                                                            
            if iszero(output_length)
                output_length = prod(size(aoutput))::Int
            end
                
            output_length == prod(size(aoutput)) || 
                error("$(PB.fullname(m)) output Variables are not all the same size")           
        end
              
        # define transpose of transfer matrix (for sparse CSR efficiency, so output = input * transfer_matrix_tr)
        rj = m.reaction
        if rj.pars.transfer_matrix.v == "Custom"
            # will be supplied (not yet implemented)
            @error "$(PB.fullname(m)) not implemented  Custom transfer matrix * $(rj.pars.transfer_multiplier.v)"
        elseif rj.pars.transfer_matrix.v == "Identity"
            input_length == output_length || 
                error("$(PB.fullname(m)) Identity transfer_matrix but Variable lengths don't match (input, output) $(input_length) != $(output_length)")
            rj.transfer_matrix_tr = SparseArrays.sparse(LinearAlgebra.I, input_length, output_length).*rj.pars.transfer_multiplier.v
            @info "    using Identity transfer matrix * $(rj.pars.transfer_multiplier.v)"
        elseif rj.pars.transfer_matrix.v == "Distribute"
            rj.transfer_matrix_tr = SparseArrays.sparse(ones(input_length, output_length)./output_length.*rj.pars.transfer_multiplier.v)
            @info "    using Distribute transfer matrix size (input, output) $(size(rj.transfer_matrix_tr)) * $(rj.pars.transfer_multiplier.v)"
        else
            error("$(PB.fullname(m)) unknown transfer_matrix='$(rj.pars.transfer_matrix.v)'")
        end
    else
        @warn "ReactionMethod $(PB.fullname(m)) no active fluxes found"
    end
     
    rt = ([v for v in input_vardata_active], [v for v in output_vardata_active]) # recreate Vectors to make them typed
    
    return rt
end


function do_transfer(
    m::PB.ReactionMethod, 
    (input_vardata, output_vardata), 
    cellrange::PB.AbstractCellRange,
    deltat
)
    tm_tr = m.reaction.transfer_matrix_tr

    # _transfer(acinput::Nothing, acoutput::Nothing, tm_tr, cr) = nothing

    function _transfer(acinput, acoutput, tm_tr, cr)
        if length(acoutput) == 1
            # Scalar Variables, (we have already checked size(tm_tr, 2) == 1)
            @inbounds for idx in SparseArrays.nzrange(tm_tr, 1)
                i = tm_tr.rowval[idx]
                acoutput[1] += acinput[i]*tm_tr.nzval[idx]
            end
        else
            for j in cr.indices
                @inbounds for idx in SparseArrays.nzrange(tm_tr, j)
                    i = tm_tr.rowval[idx]
                    acoutput[j] += acinput[i]*tm_tr.nzval[idx]
                end
            end
        end

        return nothing
    end

    for (acinput, acoutput) in zip(input_vardata, output_vardata)
        _transfer(acinput, acoutput, tm_tr, cellrange)
    end
   
    return nothing
end

end
