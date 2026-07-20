"""Material data used by the two-dimensional plane-strain RLM solver."""
Base.@kwdef struct RLMMaterialConfig
    E::Float64 = 25_840.0
    nu::Float64 = 0.18
    G_c::Float64 = 0.65
    ell::Float64 = 10.0
    kappa::Float64 = 1.0e-6
    viscosity::Float64 = 0.05
end

"""Mesh input and quadrature choices for a Q1 quadrilateral discretization."""
Base.@kwdef struct RLMMeshConfig
    path::String = "l_shape.msh"
    quadrature_order::Int = 2
end

"""Dimensionless, continuous piecewise-linear history on physical time."""
struct RLMPiecewiseLinearHistory
    times::Vector{Float64}
    values::Vector{Float64}
end

function RLMPiecewiseLinearHistory(times::AbstractVector{<:Real}, values::AbstractVector{<:Real})
    length(times) == length(values) || throw(ArgumentError("history times and values must have equal length"))
    length(times) >= 2 || throw(ArgumentError("a history requires at least two nodes"))
    result = RLMPiecewiseLinearHistory(Float64.(times), Float64.(values))
    all(isfinite, result.times) || throw(ArgumentError("history times must be finite"))
    all(isfinite, result.values) || throw(ArgumentError("history values must be finite"))
    all(diff(result.times) .> 0.0) || throw(ArgumentError("history times must be strictly increasing"))
    return result
end

@inline function history_value(history::RLMPiecewiseLinearHistory, time::Real)
    t = Float64(time)
    first(history.times) <= t <= last(history.times) || throw(ArgumentError(
        "time $t lies outside history interval [$(first(history.times)), $(last(history.times))]",
    ))
    index = searchsortedlast(history.times, t)
    index == length(history.times) && return history.values[end]
    t0, t1 = history.times[index], history.times[index + 1]
    v0, v1 = history.values[index], history.values[index + 1]
    return v0 + (t - t0) * (v1 - v0) / (t1 - t0)
end

"""Fixed support, reference amplitudes, and independent physical-time histories."""
Base.@kwdef struct RLMLoadConfig
    fixed_boundary::String = "top"
    loaded_boundary::String = "right"
    component::Int = 2
    overlap_policy::Symbol = :loaded
    displacement_amplitude::Float64 = -0.01
    displacement_history::RLMPiecewiseLinearHistory = RLMPiecewiseLinearHistory([0.0, 1.0], [0.0, 1.0])
    initial_damage::Float64 = 0.0
    body_force::NTuple{2, Float64} = (0.0, 0.0)
    body_force_history::RLMPiecewiseLinearHistory = RLMPiecewiseLinearHistory([0.0, 1.0], [0.0, 0.0])
    traction_boundary::Union{Nothing, String} = nothing
    traction::NTuple{2, Float64} = (0.0, 0.0)
    traction_history::RLMPiecewiseLinearHistory = RLMPiecewiseLinearHistory([0.0, 1.0], [0.0, 0.0])
end

"""Fixed physical-time BDF1 grid and RLM energy parameter."""
Base.@kwdef struct RLMTimeConfig
    final_time::Float64 = 1.0
    dt::Float64 = 1.0e-3
    alpha::Float64 = 1.0
end

"""All floating-point decisions made by the first-stage RLM implementation."""
Base.@kwdef struct RLMToleranceConfig
    principal_zero_abs::Float64 = 1.0e-14
    principal_zero_rel::Float64 = 1.0e-12
    repeated_eigen_abs::Float64 = 1.0e-14
    repeated_eigen_rel::Float64 = 1.0e-12

    discriminant_abs::Float64 = 0.0
    discriminant_rel::Float64 = 1.0e-12
    coefficient_abs::Float64 = eps(Float64)
    coefficient_rel::Float64 = 1.0e-14
    positive_root::Float64 = 0.0
    duplicate_root_rel::Float64 = 1.0e-12
    scalar_denominator_epsilon::Float64 = eps(Float64)
    scalar_residual::Float64 = 1.0e-10

    c1_abs::Float64 = 1.0e-12
    c1_rel::Float64 = 1.0e-10
    branch_identity::Float64 = 1.0e-9
    energy_balance_abs::Float64 = 1.0e-10
    energy_balance_rel::Float64 = 1.0e-8

end

"""Filesystem output controls. Diagnostics are always retained in memory."""
Base.@kwdef struct RLMOutputConfig
    directory::String = "data/sims/rlm_bdf1"
    write_csv::Bool = true
    write_vtk::Bool = true
    vtk_every_time_step::Int = 1
    verbose::Bool = true
end

"""Complete, explicit configuration for Miehe--RLM-PE--BDF1."""
Base.@kwdef struct RLMConfig
    material::RLMMaterialConfig = RLMMaterialConfig()
    mesh::RLMMeshConfig = RLMMeshConfig()
    load::RLMLoadConfig = RLMLoadConfig()
    time::RLMTimeConfig = RLMTimeConfig()
    tolerances::RLMToleranceConfig = RLMToleranceConfig()
    output::RLMOutputConfig = RLMOutputConfig()
end

"""Accepted RLM physical state at one real-time instant."""
mutable struct RLMState
    u::Vector{Float64}
    d::Vector{Float64}
    q::Float64
    P::Float64
end

Base.copy(state::RLMState) = RLMState(copy(state.u), copy(state.d), state.q, state.P)

"""One row of the mandatory RLM diagnostics."""
Base.@kwdef mutable struct RLMDiagnostic
    step::Int
    time::Float64
    dt::Float64
    displacement_factor::Float64
    body_force_factor::Float64
    traction_factor::Float64
    displacement_rate::Float64
    body_force_rate::Float64
    traction_rate::Float64
    displacement::Float64
    accepted::Bool
    status::String
    raw_energy::Float64 = NaN
    internal_energy::Float64 = NaN
    proxy_energy::Float64 = NaN
    elastic_energy::Float64 = NaN
    positive_elastic_energy::Float64 = NaN
    negative_elastic_energy::Float64 = NaN
    fracture_energy::Float64 = NaN
    nonlinear_energy::Float64 = NaN
    predicted_energy::Float64 = NaN
    prediction_gap::Float64 = NaN
    proxy_gap::Float64 = NaN

    q::Float64 = NaN
    q_minus_one::Float64 = NaN
    c0::Float64 = NaN
    c1::Float64 = NaN
    A::Float64 = NaN
    B::Float64 = NaN
    C::Float64 = NaN
    discriminant::Float64 = NaN
    discriminant_used::Float64 = NaN
    scalar_residual::Float64 = NaN
    reaction_force::Float64 = NaN
    external_work::Float64 = NaN
    cumulative_external_work::Float64 = NaN
    viscous_dissipation::Float64 = NaN
    numerical_dissipation::Float64 = NaN
    cumulative_viscous_dissipation::Float64 = NaN
    cumulative_numerical_dissipation::Float64 = NaN

    phase_increment::Float64 = NaN
    phase_relative_increment::Float64 = NaN
    phase_equilibrium_residual::Float64 = NaN
    healing::Float64 = NaN
    min_d::Float64 = NaN
    max_d::Float64 = NaN
    energy_balance_residual::Float64 = NaN
end

"""Transactional failure: no state has been committed when this exception is made."""
struct RLMStepFailure <: Exception
    code::Symbol
    message::String
    data::NamedTuple
end

function Base.showerror(io::IO, err::RLMStepFailure)
    print(io, "RLM step failed [", err.code, "]: ", err.message)
end

"""Four affine branches and all scalar quantities for an accepted candidate step."""
struct RLMTrial
    u_a::Vector{Float64}
    u_b::Vector{Float64}
    d_a::Vector{Float64}
    d_b::Vector{Float64}
    u_star::Vector{Float64}
    d_star::Vector{Float64}
    u::Vector{Float64}
    d::Vector{Float64}
    n_u::Vector{Float64}
    n_d::Vector{Float64}
    P::Float64
    c0::Float64
    c1::Float64
    A::Float64
    B::Float64
    C::Float64
    discriminant::Float64
    discriminant_used::Float64
    q::Float64
    scalar_residual::Float64
    raw_energy::Float64
    proxy_energy::Float64
    elastic_energy::Float64
    positive_elastic_energy::Float64
    negative_elastic_energy::Float64
    fracture_energy::Float64
    nonlinear_energy::Float64
    reaction_force::Float64
    internal_energy::Float64
    external_work::Float64
    viscous_dissipation::Float64
    numerical_dissipation::Float64
    energy_balance_residual::Float64
    phase_increment::Float64
    phase_relative_increment::Float64
    phase_equilibrium_residual::Float64
    healing::Float64
    min_d::Float64
    max_d::Float64
end

"""Final status and complete accepted/failure diagnostic history."""
struct RLMResult{P}
    success::Bool
    completed::Bool
    message::String
    state::RLMState
    diagnostics::Vector{RLMDiagnostic}
    problem::P
end

@inline lame_lambda(mat::RLMMaterialConfig) =
    mat.E * mat.nu / ((1.0 + mat.nu) * (1.0 - 2.0 * mat.nu))

@inline lame_mu(mat::RLMMaterialConfig) = mat.E / (2.0 * (1.0 + mat.nu))

function validate_config(config::RLMConfig)
    mat = config.material
    mesh = config.mesh
    load = config.load
    time = config.time
    tol = config.tolerances
    output = config.output

    all(isfinite, (mat.E, mat.nu, mat.G_c, mat.ell, mat.kappa, mat.viscosity)) ||
        throw(ArgumentError("all material parameters must be finite"))
    mat.E > 0.0 || throw(ArgumentError("material.E must be positive"))
    (-1.0 < mat.nu < 0.5) || throw(ArgumentError("plane-strain material.nu must lie in (-1, 0.5)"))
    mat.G_c > 0.0 || throw(ArgumentError("material.G_c must be positive"))
    mat.ell > 0.0 || throw(ArgumentError("material.ell must be positive"))
    (0.0 < mat.kappa < 1.0) || throw(ArgumentError("material.kappa must lie in (0, 1)"))
    mat.viscosity > 0.0 || throw(ArgumentError("material.viscosity η must be positive"))

    mesh.quadrature_order >= 2 || throw(ArgumentError("Q1 RLM quadrature_order must be at least 2"))
    load.component in (1, 2) || throw(ArgumentError("load.component must be 1 or 2"))
    load.overlap_policy in (:loaded, :fixed) ||
        throw(ArgumentError("load.overlap_policy must be :loaded or :fixed"))
    isfinite(load.displacement_amplitude) ||
        throw(ArgumentError("load.displacement_amplitude must be finite"))
    isfinite(load.initial_damage) || throw(ArgumentError("load.initial_damage must be finite"))
    all(isfinite, load.body_force) || throw(ArgumentError("load.body_force must be finite"))
    all(isfinite, load.traction) || throw(ArgumentError("load.traction must be finite"))

    all(isfinite, (time.final_time, time.dt, time.alpha)) ||
        throw(ArgumentError("time.final_time, time.dt and time.alpha must be finite"))
    time.final_time > 0.0 || throw(ArgumentError("time.final_time must be positive"))
    time.dt > 0.0 || throw(ArgumentError("time.dt must be positive"))
    time.alpha > 0.0 || throw(ArgumentError("time.alpha must be a positive energy"))
    nsteps = time.final_time / time.dt
    isapprox(nsteps, round(nsteps); atol = 1.0e-10, rtol = 1.0e-12) ||
        throw(ArgumentError("time.final_time must be an integer multiple of time.dt"))
    histories = (load.displacement_history, load.body_force_history, load.traction_history)
    for history in histories
        first(history.times) <= 0.0 && last(history.times) >= time.final_time ||
            throw(ArgumentError("each load history must cover [0, time.final_time]"))
        history_value(history, 0.0) == 0.0 ||
            throw(ArgumentError("each load history must start from zero at t=0"))
    end

    for name in fieldnames(RLMToleranceConfig)
        value = getfield(tol, name)
        value >= 0.0 || throw(ArgumentError("tolerances.$name must be nonnegative"))
        isfinite(value) || throw(ArgumentError("tolerances.$name must be finite"))
    end
    output.vtk_every_time_step > 0 ||
        throw(ArgumentError("output.vtk_every_time_step must be positive"))
    return config
end
