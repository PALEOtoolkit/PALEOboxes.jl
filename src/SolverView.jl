"""
    SolverView

Provides a view on the whole or some part of the Model for a numerical solver.
Contains [`VariableAggregator`](@ref)s for a subset of spatial locations 
(Domains, indices within spatial Domains) and Variables, with:
- ODE paired `stateexplicit` (S) and `stateexplicit_deriv` (dS/dt), where dS/dt = F(S).
- Implicit-ODE paired `total` (T) and `total_deriv` (dT/dt), where dT(S)/dt = F(T(S)) with
  `total` a function of explicit and implicit state Variables `statexplicit` and `state` (S).
- Algebraic `constraint`s (C), where C(S) = 0 with C a function of explicit and implicit state Variables `statexplicit` and `state` (S),

The number of `total` + number of `constraint` Variables must equal the number of implicit `state` Variables.

Optional access methods provide an `ODE/DAE solver` view with composite `statevar` and `statevar_sms`,
where:

- `statevar` is a concatenation of `stateexplicit` and `state` ([`set_statevar!`](@ref))
- `statevar_sms` is a concatenation of `stateexplicit_deriv`, `total_deriv`, `constraints` ([`get_statevar_sms!`](@ref))
"""
mutable struct SolverView{
    T, 
    VA1 <: VariableAggregator, VA2 <: VariableAggregator, VA3 <: VariableAggregator, VA4 <: VariableAggregator, VA5 <: VariableAggregator, VA6 <: VariableAggregator,
    VA7 <: VariableAggregatorNamed
}
    stateexplicit::VA1
    stateexplicit_deriv::VA2
    stateexplicit_norm::Vector{Float64}
    
    total::VA3
    total_deriv::VA4
    total_norm::Vector{Float64}

    constraints::VA5
    constraints_norm::Vector{Float64}

    state::VA6
    state_norm::Vector{Float64}

    hostdep::VA7

    function SolverView(
        etype::Type{T};  # eltype(modeldata)
        stateexplicit::A1, stateexplicit_deriv::A2, total::A3, total_deriv::A4, constraints::A5, state::A6, hostdep::A7
    ) where {T, A1, A2, A3, A4, A5, A6, A7}
        return new{T, A1, A2, A3, A4, A5, A6, A7}(
            stateexplicit,
            stateexplicit_deriv,
            Vector{Float64}(),
            total,
            total_deriv,
            Vector{Float64}(),
            constraints,
            Vector{Float64}(),
            state,
            Vector{Float64}(),
            hostdep
        )
    end
end

Base.eltype(::Type{SolverView{T, A1, A2, A3, A4, A5, A6, A7}}) where {T, A1, A2, A3, A4, A5, A6, A7} = T

"compact form"
function Base.show(io::IO, sv::SolverView)
    print(
        io,
        "SolverView{$(eltype(sv)), ...}",
        "(stateexplicit $(length(sv.stateexplicit)), "*
        "total $(length(sv.total)), "*
        "constraints $(length(sv.constraints)), "*
        "state $(length(sv.state)), "*
        "hostdep $(length(sv.hostdep)))"
    )
    return nothing
end
"multiline form"
function Base.show(io::IO, m::MIME"text/plain", sv::SolverView)
    println(io, "SolverView{$(eltype(sv)), ...}:")
    for f in (:stateexplicit, :stateexplicit_deriv, :total, :total_deriv, :constraints, :state, :hostdep)
        print(io, String(f)*":  ")
        va = getfield(sv, f)
        show(io, m, va)
    end
    return nothing
end


"get number of algebraic constraints (0 if problem is an ODE)"
function num_algebraic_constraints(sv::SolverView)
    return length(sv.constraints)
end

"get number of implicit total variables with time derivatives"
function num_total(sv::SolverView)
    return length(sv.total)
end

"""
    set_statevar!(sv::SolverView, u)

Set combined stateexplicit, state variables from u
"""
function set_statevar!(sv::SolverView, u::AbstractVector)
   
    l_ts = copyto!(sv.stateexplicit, u, sof=1)
    copyto!(sv.state, u, sof=l_ts+1)      
   
    return nothing
end

"statevar += a*u"
function add_statevar!(sv::SolverView, a, u::AbstractVector)
   
    l_ts = length(sv.stateexplicit)   
    add_data!(sv.stateexplicit, a, view(u, 1:l_ts))
    add_data!(sv.state, a, view(u, (l_ts+1):length(u)))      
   
    return nothing
end

function get_statevar!(u, sv::SolverView)
    l_ts = copyto!(u, sv.stateexplicit, dof=1)
    copyto!(u, sv.state, dof=1+l_ts)
    return nothing
end



function state_vars_isdifferential(sv::SolverView)
    isdifferential = trues(length(sv.stateexplicit) + length(sv.total) + length(sv.constraints))
    isdifferential[end-length(sv.constraints)+1:end] .= false
    return isdifferential
end

"""
    get_statevar_sms!(du, sv::SolverView)

Get combined derivatives and constraints, eg for an ODE solver
"""
function get_statevar_sms!(du::AbstractVector, sv::SolverView)
    l_ts = copyto!(du, sv.stateexplicit_deriv, dof=1)

    # TODO get total_deriv  
    # only used for initialisation, doesn't really make sense for an ODE solver ??
    l_ti = copyto!(du, sv.total_deriv, dof=l_ts+1)

    # add constraints to derivative (for ODE solvers that handle this using mass_matrix = 0)
    copyto!(du, sv.constraints, dof=l_ts+l_ti+1)

    return nothing
end

function get_statevar_sms(sv::SolverView)
    du = Vector{eltype(sv)}(undef, length(sv.stateexplicit_deriv) + length(sv.total_deriv) + length(sv.constraints))
    
    get_statevar_sms!(du, sv)

    return du
end

function get_statevar(sv::SolverView)
    return vcat(get_data(sv.stateexplicit),
                get_data(sv.state))
end

function get_statevar_norm(sv::SolverView)
    return vcat(sv.stateexplicit_norm,
                sv.state_norm)
end

function get_statevar_sms_norm(sv::SolverView)
    # use stateexplicit and total as estimates of the time derivative normalisation value
    return vcat(sv.stateexplicit_norm, sv.total_norm, sv.constraints_norm)
end

"copy norm values from state variable etc data"
function copy_norm!(sv::SolverView)
    sv.stateexplicit_norm   = value_ad.(get_data(sv.stateexplicit))
    sv.state_norm           = value_ad.(get_data(sv.state))
    sv.total_norm           = value_ad.(get_data(sv.total))
    sv.constraints_norm     = value_ad.(get_data(sv.constraints))
    return nothing
end

"copy norm values back into state variable etc data"
function uncopy_norm!(sv::SolverView)
    copyto!(sv.stateexplicit, sv.stateexplicit_norm)
    copyto!(sv.state, sv.state_norm)
    copyto!(sv.total, sv.total_norm)
    copyto!(sv.constraints, sv.constraints_norm)
   
    return nothing
end

"""
    set_tforce!(sv::SolverView, t)

Set `global.tforce` model time (if defined) from t.
"""
set_tforce!(sv::SolverView, t) = set_tforce!(sv.hostdep, t, allow_missing=true)
