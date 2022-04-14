module ReactionRateStoichMockModule

import PALEOboxes as PB

"Mock Reaction to test RateStoich"
Base.@kwdef mutable struct ReactionRateStoichMock <: PB.AbstractReaction
    base::PB.ReactionBase

    # a biogeochemically implausible fractionating H2S oxidation 
    # H2S + 2O2 -> SO4,   with -10 per mil fractionation
    ratestoich = PB.RateStoich(
        PB.VarProp("myrate", "mol yr-1", "a rate"),
        ((-2.0, "O2"), (-1.0,"H2S::Isotope"), (1.0, "SO4::Isotope")),
        deltavarname_eta=("H2S_delta", -10.0),
        processname="redox"
    )

end

function PB.register_methods!(rj::ReactionRateStoichMock)

    PB.add_method_do!(rj, do_rate, (PB.VarList_single(rj.ratestoich.ratevartemplate),))
    
    PB.add_method_do!(rj, rj.ratestoich, isotope_data=PB.IsotopeLinear)

    return nothing
end

function do_rate(m, (rate, ), cellrange::PB.AbstractCellRange, deltat)

    rate .= 42.0
  
    return nothing
end

"Install create_reaction when module imported"
function __init__()
    PB.add_reaction_factory(ReactionRateStoichMock)
    return nothing
end


end # module
