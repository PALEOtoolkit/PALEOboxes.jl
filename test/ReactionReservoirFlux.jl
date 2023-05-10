module ReactionRateStoichMockModule

import PALEOboxes as PB

"Mock Reaction to test calculating a flux that depends on a reservoir"
Base.@kwdef mutable struct ReactionReservoirFlux <: PB.AbstractReaction
    base::PB.ReactionBase

end

function PB.register_methods!(rj::ReactionReservoirFlux)

    vars = [
        PB.VarDepScalar("R", "mol", "a scalar reservoir"),
        PB.VarContribScalar("R_sms", "mol yr-1", "a scalar reservoir source-sink"),
    ]
    PB.add_method_do!(rj, do_reservoirflux, (PB.VarList_namedtuple(vars),))

    return nothing
end

function do_reservoirflux(m, (vars, ), cellrange::PB.AbstractCellRange, deltat)

    flux = 1.0*vars.R[] # first order rate
    vars.R_sms[] -= flux
  
    return nothing
end


end # module
