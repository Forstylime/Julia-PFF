"""Material data used by the two-dimensional plane-strain RLM solver."""
Base.@kwdef struct RLMMaterialConfig
    E::Float64 = 25_840.0
    nu::Float64 = 0.18
    G_c::Float64 = 0.65
    ell::Float64 = 10.0
    kappa::Float64 = 1.0e-6
    mobility::Float64 = 20.0
end

"""Mesh input and quadrature choices for a Q1 quadrilateral discretization."""
Base.@kwdef struct RLMMeshConfig
    path::String = "l_shape.msh"
    quadrature_order::Int = 2
end

"""Fixed support, displacement loading, and fixed external loads."""
Base.@kwdef struct RLMLoadConfig
    fixed_boundary::String = "top"
    loaded_boundary::String = "right"
    component::Int = 2
    overlap_policy::Symbol = :loaded
    final_displacement::Float64 = -0.01
    load_steps::Int = 20
    initial_damage::Float64 = 0.0
    body_force::NTuple{2, Float64} = (0.0, 0.0)
    traction_boundary::Union{Nothing, String} = nothing
    traction::NTuple{2, Float64} = (0.0, 0.0)
end

"""BDF1 step size, fixed RLM energy parameter, and relaxation controls."""
Base.@kwdef struct RLMTimeConfig
    dt::Float64 = 1.0e-3
    alpha::Float64 = 1.0
    relaxation_mode::Symbol = :to_tolerance
    min_relax_steps::Int = 1
    max_relax_steps::Int = 100
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

    phase::Float64 = 1.0e-6
    q::Float64 = 1.0e-6
end

"""Filesystem output controls. Diagnostics are always retained in memory."""
Base.@kwdef struct RLMOutputConfig
    directory::String = "data/sims/rlm_bdf1"
    write_csv::Bool = true
    write_vtk::Bool = true
    vtk_every_load_step::Int = 1
    verbose::Bool = true
end

"""Complete, explicit configuration for Miehe--RLM-PC--BDF1."""
Base.@kwdef struct RLMConfig
    material::RLMMaterialConfig = RLMMaterialConfig()
    mesh::RLMMeshConfig = RLMMeshConfig()
    load::RLMLoadConfig = RLMLoadConfig()
    time::RLMTimeConfig = RLMTimeConfig()
    tolerances::RLMToleranceConfig = RLMToleranceConfig()
    output::RLMOutputConfig = RLMOutputConfig()
end

"""Accepted RLM history at one fixed-load relaxation instant."""
mutable struct RLMState
    u::Vector{Float64}
    d::Vector{Float64}
    q::Float64
    P::Float64
end

Base.copy(state::RLMState) = RLMState(copy(state.u), copy(state.d), state.q, state.P)

"""One row of the mandatory RLM diagnostics."""
Base.@kwdef struct RLMDiagnostic
    load_step::Int
    relax_step::Int
    load_fraction::Float64
    displacement::Float64
    accepted::Bool
    status::String

    raw_energy::Float64 = NaN
    proxy_energy::Float64 = NaN
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

    phase_increment::Float64 = NaN
    phase_relative_increment::Float64 = NaN
    healing::Float64 = NaN
    min_d::Float64 = NaN
    max_d::Float64 = NaN
    energy_balance_residual::Float64 = NaN
    reset_jump::Float64 = NaN
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
    energy_balance_residual::Float64
    phase_increment::Float64
    phase_relative_increment::Float64
    healing::Float64
    min_d::Float64
    max_d::Float64
end

"""Final status and complete accepted/failure diagnostic history."""
struct RLMResult{P}
    success::Bool
    converged::Bool
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

    all(isfinite, (mat.E, mat.nu, mat.G_c, mat.ell, mat.kappa, mat.mobility)) ||
        throw(ArgumentError("all material parameters must be finite"))
    mat.E > 0.0 || throw(ArgumentError("material.E must be positive"))
    (-1.0 < mat.nu < 0.5) || throw(ArgumentError("plane-strain material.nu must lie in (-1, 0.5)"))
    mat.G_c > 0.0 || throw(ArgumentError("material.G_c must be positive"))
    mat.ell > 0.0 || throw(ArgumentError("material.ell must be positive"))
    (0.0 < mat.kappa < 1.0) || throw(ArgumentError("material.kappa must lie in (0, 1)"))
    mat.mobility > 0.0 || throw(ArgumentError("material.mobility M must be positive"))

    mesh.quadrature_order >= 2 || throw(ArgumentError("Q1 RLM quadrature_order must be at least 2"))
    load.component in (1, 2) || throw(ArgumentError("load.component must be 1 or 2"))
    load.overlap_policy in (:loaded, :fixed) ||
        throw(ArgumentError("load.overlap_policy must be :loaded or :fixed"))
    load.load_steps > 0 || throw(ArgumentError("load.load_steps must be positive"))
    isfinite(load.final_displacement) ||
        throw(ArgumentError("load.final_displacement must be finite"))
    isfinite(load.initial_damage) || throw(ArgumentError("load.initial_damage must be finite"))
    all(isfinite, load.body_force) || throw(ArgumentError("load.body_force must be finite"))
    all(isfinite, load.traction) || throw(ArgumentError("load.traction must be finite"))

    all(isfinite, (time.dt, time.alpha)) ||
        throw(ArgumentError("time.dt and time.alpha must be finite"))
    time.dt > 0.0 || throw(ArgumentError("time.dt must be positive"))
    time.alpha > 0.0 || throw(ArgumentError("time.alpha must be a positive energy"))
    time.relaxation_mode in (:to_tolerance, :fixed_steps) || throw(ArgumentError(
        "time.relaxation_mode must be :to_tolerance or :fixed_steps",
    ))
    time.min_relax_steps >= 1 || throw(ArgumentError("time.min_relax_steps must be at least one"))
    time.max_relax_steps >= time.min_relax_steps ||
        throw(ArgumentError("time.max_relax_steps must not be smaller than min_relax_steps"))

    for name in fieldnames(RLMToleranceConfig)
        value = getfield(tol, name)
        value >= 0.0 || throw(ArgumentError("tolerances.$name must be nonnegative"))
        isfinite(value) || throw(ArgumentError("tolerances.$name must be finite"))
    end
    output.vtk_every_load_step > 0 ||
        throw(ArgumentError("output.vtk_every_load_step must be positive"))
    return config
end
