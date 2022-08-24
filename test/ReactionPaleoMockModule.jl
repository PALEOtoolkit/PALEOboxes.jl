module ReactionPaleoMockModule

import PALEOboxes as PB

# using Infiltrator # Julia debugger

"Define variables, parameters, and any additional reaction state"
Base.@kwdef mutable struct ReactionPaleoMock{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDouble( "par1", 0.0, units="m-3", 
            description="test parameter double"),
        PB.ParInt(    "parint", 2, 
            description="test parameter int"),
        PB.ParBool(   "parbool", true, 
            description="test parameter bool"),
        PB.ParString( "parstring", "a_string", 
            description="test parameter string"),
    )

end


function do_stateandeqb(m::PB.ReactionMethod, pars, (vars, ), cellrange::PB.AbstractCellRange, deltat)

    vars.scalar_prop[] = vars.scalar_dep[]  # 25nS, 0 allocations
 
    return nothing
end

function do_react(m::PB.ReactionMethod, pars, (vars, ), cellrange::PB.AbstractCellRange, deltat)
 
    @inbounds  for i in cellrange.indices
        vars.phy[i] = pars.par1[]
    end

    return nothing
end


function PB.register_methods!(rj::ReactionPaleoMock)
    println("ReactionPaleoMock.register_methods! rj=", rj)

    println("  parint=", rj.pars.parint[])
    println("  parbool=", rj.pars.parbool[])
    println("  parstring=", rj.pars.parstring[])
  
    vars_stateandeqb = [
        PB.VarDepScalar("scalar_dep", "mol m-3", "a scalar dependency"),
        PB.VarPropScalar("%reaction%scalar_prop", "mol m-3", "a scalar property")
    ]

    PB.add_method_do!(rj, do_stateandeqb, (PB.VarList_namedtuple(vars_stateandeqb),))

    # default vector Variable
    vars_react  = [
        PB.VarProp(   "phy", "mol m-3", "phytoplankton concentration")
    ]
    
    PB.add_method_do!(rj, do_react, (PB.VarList_namedtuple(vars_react),))
    
    return nothing
end

# contiguous range is faster than using indices array
# @btime ReactionPaleoMockModule.do_react_kernel(reaction.var_phy._data, 1:1000, reaction.pars.par1.value)
#  168.262 ns (1 allocation: 16 bytes)
function  do_react_kernel(v,  indices, p)
    for i in indices
        @inbounds v[i] = p
    end
    return nothing
end


end # module
