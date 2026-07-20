@inline function _named_value(data::NamedTuple, name::Symbol, default = NaN)
    return haskey(data, name) ? getfield(data, name) : default
end

function _throw_step_failure(code, message; kwargs...)
    throw(RLMStepFailure(code, message, (; kwargs...)))
end

function _solve_u_baseline(problem::RLMProblem)
    rhs_u_a = copy(problem.f_ext)
    apply_rhs!(problem.rhs_data_u, rhs_u_a, problem.ch_u, false)
    u_a = problem.factor_u \ rhs_u_a
    apply!(u_a, problem.ch_u)
    return u_a
end

function _solve_affine_branches(problem::RLMProblem, state::RLMState, n_u, n_d, u_a)
    viscosity_over_dt = problem.config.material.viscosity / problem.config.time.dt
    rhs_u_b = -n_u
    apply_rhs!(problem.rhs_data_u, rhs_u_b, problem.ch_u, true)
    u_b = problem.factor_u \ rhs_u_b
    apply_zero!(u_b, problem.ch_u)
    d_a = problem.factor_d \ (viscosity_over_dt * (problem.M_d * state.d))
    d_b = problem.factor_d \ (-n_d)
    return u_a, u_b, d_a, d_b
end

function _history_step_data(load::RLMLoadConfig, t_old::Real, t_new::Real)
    dt = Float64(t_new - t_old)
    dt > 0.0 || throw(ArgumentError("time step must be positive"))
    d_old, d_new = history_value(load.displacement_history, t_old), history_value(load.displacement_history, t_new)
    b_old, b_new = history_value(load.body_force_history, t_old), history_value(load.body_force_history, t_new)
    tr_old, tr_new = history_value(load.traction_history, t_old), history_value(load.traction_history, t_new)
    return (; displacement_factor = d_new, body_force_factor = b_new, traction_factor = tr_new,
        displacement_rate = (d_new - d_old) / dt, body_force_rate = (b_new - b_old) / dt,
        traction_rate = (tr_new - tr_old) / dt)
end

function _dirichlet_increment(problem::RLMProblem, t_old::Real, t_new::Real)
    u_old = zeros(ndofs(problem.dh_u))
    update!(problem.ch_u, t_old)
    apply!(u_old, problem.ch_u)
    u_new = zeros(ndofs(problem.dh_u))
    update!(problem.ch_u, t_new)
    apply!(u_new, problem.ch_u)
    return u_new - u_old
end

"""
    compute_rlm_bdf1_trial(problem, state, t_next)

Compute one real-time BDF1 candidate from `t_next - config.time.dt` to `t_next`.
The accepted `state` is never mutated; the problem's constraint and external-force
data are advanced to `t_next` so that a caller may commit the returned candidate.
"""
function compute_rlm_bdf1_trial(problem::RLMProblem, state::RLMState, t_next::Real)
    config = problem.config
    t_new = Float64(t_next)
    t_old = t_new - config.time.dt
    step_data = _history_step_data(config.load, t_old, t_new)
    update!(problem.ch_u, t_new)
    update_rlm_external_force!(problem, t_new)

    tolerances = config.tolerances
    alpha = config.time.alpha
    n_u = zeros(ndofs(problem.dh_u))
    n_d = zeros(ndofs(problem.dh_d))
    assemble_rlm_nonlinear_forces!(n_u, n_d, problem, state.u, state.d)
    u_a = _solve_u_baseline(problem)
    u_a, u_b, d_a, d_b = _solve_affine_branches(problem, state, n_u, n_d, u_a)

    u_star = u_a + u_b
    d_star = d_a + d_b
    P_next = rlm_nonlinear_energy(problem, u_star, d_star)
    c0 = dot(n_u, u_a - state.u) + dot(n_d, d_a - state.d)
    c1 = dot(n_u, u_b) + dot(n_d, d_b)
    c1_from_branches = -dot(u_b, problem.K_u * u_b) - dot(d_b, problem.K_d * d_b)
    branch_scale = abs(c1) + abs(c1_from_branches) + eps(Float64)
    branch_residual = abs(c1 - c1_from_branches) / branch_scale
    branch_residual <= tolerances.branch_identity || _throw_step_failure(
        :c1_branch_identity, "c1 does not match the negative affine-branch quadratic form";
        P = P_next, c0, c1, c1_from_branches, branch_residual,
    )
    c1_tolerance = tolerances.c1_abs + tolerances.c1_rel * max(abs(c1), abs(c1_from_branches))
    c1 <= c1_tolerance || _throw_step_failure(
        :positive_c1, "c1=$c1 is significantly positive (tolerance=$c1_tolerance)"; P = P_next, c0, c1,
    )

    A = alpha - c1
    B = -c0
    C = P_next - state.P - alpha * state.q^2
    root = solve_rlm_quadratic(A, B, C, state.q, tolerances)
    root.success || _throw_step_failure(root.code, root.message;
        P = P_next, c0, c1, A, B, C, discriminant = root.discriminant,
        discriminant_used = root.discriminant_used, scalar_residual = root.residual,
    )

    q_next = root.q
    u_next = u_a + q_next * u_b
    d_next = d_a + q_next * d_b
    apply!(u_next, problem.ch_u)

    Ku_u = problem.K_u * u_next
    Kat2_d = problem.K_AT2 * d_next
    delta_u = u_next - state.u
    delta_d = d_next - state.d
    viscosity_over_dt = config.material.viscosity / config.time.dt
    nonlinear_actual = rlm_nonlinear_energy(problem, u_next, d_next)
    positive_elastic_energy, negative_elastic_energy = rlm_elastic_split_energies(problem, u_next, d_next)
    elastic_energy = positive_elastic_energy + negative_elastic_energy
    fracture_energy = 0.5 * dot(d_next, Kat2_d)
    internal_energy = elastic_energy + fracture_energy
    raw_energy = internal_energy - dot(problem.f_ext, u_next)
    proxy_energy = 0.5 * dot(u_next, Ku_u) + fracture_energy + P_next + alpha * (q_next^2 - 1.0)
    proxy_old = 0.5 * dot(state.u, problem.K_u * state.u) +
                0.5 * dot(state.d, problem.K_AT2 * state.d) + state.P + alpha * (state.q^2 - 1.0)
    viscous_dissipation = viscosity_over_dt * dot(delta_d, problem.M_d * delta_d)
    numerical_dissipation = 0.5 * dot(delta_u, problem.K_u * delta_u) +
                            0.5 * dot(delta_d, problem.K_AT2 * delta_d)
    mechanical_residual = Ku_u + q_next * n_u - problem.f_ext
    dirichlet_increment = _dirichlet_increment(problem, t_old, t_new)
    external_work = dot(problem.f_ext, delta_u) + dot(mechanical_residual, dirichlet_increment)
    balance_raw = proxy_energy - proxy_old - external_work + viscous_dissipation + numerical_dissipation
    balance_scale = abs(proxy_energy) + abs(proxy_old) + abs(external_work) +
                    abs(viscous_dissipation) + abs(numerical_dissipation) + eps(Float64)
    balance_residual = abs(balance_raw) / balance_scale
    balance_tolerance = tolerances.energy_balance_abs + tolerances.energy_balance_rel * balance_scale
    abs(balance_raw) <= balance_tolerance || _throw_step_failure(
        :energy_balance, "BDF1 work-energy identity residual $(abs(balance_raw)) exceeds $balance_tolerance";
        P = P_next, c0, c1, A, B, C, discriminant = root.discriminant,
        discriminant_used = root.discriminant_used, scalar_residual = root.residual,
        energy_balance_residual = balance_residual,
    )

    phase_residual_vector = viscosity_over_dt * (problem.M_d * delta_d) + Kat2_d + q_next * n_d
    phase_residual_scale = viscosity_over_dt * norm(problem.M_d * delta_d) + norm(Kat2_d) +
                           abs(q_next) * norm(n_d) + eps(Float64)
    phase_equilibrium_residual = norm(phase_residual_vector) / phase_residual_scale
    reaction_force = sum(mechanical_residual[problem.loaded_component_dofs])
    phase_increment, phase_relative_increment, healing, min_d, max_d = phase_field_metrics(problem, d_next, state.d)
    return RLMTrial(u_a, u_b, d_a, d_b, u_star, d_star, u_next, d_next, n_u, n_d,
        P_next, c0, c1, A, B, C, root.discriminant, root.discriminant_used, q_next,
        root.residual, raw_energy, proxy_energy, elastic_energy, positive_elastic_energy, negative_elastic_energy, fracture_energy,
        nonlinear_actual, reaction_force, internal_energy, external_work, viscous_dissipation,
        numerical_dissipation, balance_residual, phase_increment, phase_relative_increment,
        phase_equilibrium_residual, healing, min_d, max_d)
end

function _diagnostic(problem, state, trial, step::Int, time::Float64, cumulative_work,
    cumulative_viscous, cumulative_numerical; status = "accepted", accepted = true)
    load = problem.config.load
    data = if step == 0
        (; displacement_factor = history_value(load.displacement_history, 0.0),
            body_force_factor = history_value(load.body_force_history, 0.0),
            traction_factor = history_value(load.traction_history, 0.0),
            displacement_rate = 0.0, body_force_rate = 0.0, traction_rate = 0.0)
    else
        _history_step_data(load, time - problem.config.time.dt, time)
    end
    return RLMDiagnostic(step = step, time = time, dt = step == 0 ? 0.0 : problem.config.time.dt,
        displacement_factor = data.displacement_factor, body_force_factor = data.body_force_factor,
        traction_factor = data.traction_factor, displacement_rate = step == 0 ? 0.0 : data.displacement_rate,
        body_force_rate = step == 0 ? 0.0 : data.body_force_rate,
        traction_rate = step == 0 ? 0.0 : data.traction_rate,
        displacement = load.displacement_amplitude * data.displacement_factor,
        accepted = accepted, status = status, raw_energy = trial.raw_energy,
        internal_energy = trial.internal_energy, proxy_energy = trial.proxy_energy,
        elastic_energy = trial.elastic_energy, positive_elastic_energy = trial.positive_elastic_energy,
        negative_elastic_energy = trial.negative_elastic_energy, fracture_energy = trial.fracture_energy,
        nonlinear_energy = trial.nonlinear_energy, predicted_energy = trial.P,
        prediction_gap = trial.P - trial.nonlinear_energy, proxy_gap = abs(trial.proxy_energy - trial.internal_energy),
        q = trial.q, q_minus_one = trial.q - 1.0, c0 = trial.c0, c1 = trial.c1, A = trial.A,
        B = trial.B, C = trial.C, discriminant = trial.discriminant,
        discriminant_used = trial.discriminant_used, scalar_residual = trial.scalar_residual,
        reaction_force = trial.reaction_force, external_work = trial.external_work,
        cumulative_external_work = cumulative_work, viscous_dissipation = trial.viscous_dissipation,
        numerical_dissipation = trial.numerical_dissipation,
        cumulative_viscous_dissipation = cumulative_viscous,
        cumulative_numerical_dissipation = cumulative_numerical,
        phase_increment = trial.phase_increment, phase_relative_increment = trial.phase_relative_increment,
        phase_equilibrium_residual = trial.phase_equilibrium_residual, healing = trial.healing,
        min_d = trial.min_d, max_d = trial.max_d, energy_balance_residual = trial.energy_balance_residual)
end

function _initial_diagnostic(problem, state)
    nonlinear = state.P
    positive, negative = rlm_elastic_split_energies(problem, state.u, state.d)
    elastic = positive + negative
    fracture = 0.5 * dot(state.d, problem.K_AT2 * state.d)
    internal = elastic + fracture
    trial = RLMTrial(state.u, zeros(length(state.u)), state.d, zeros(length(state.d)), state.u, state.d,
        state.u, state.d, zeros(length(state.u)), zeros(length(state.d)), state.P, NaN, NaN, NaN, NaN, NaN,
        NaN, NaN, state.q, NaN, internal, internal, elastic, positive, negative, fracture, nonlinear, 0.0, internal, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, minimum(state.d), maximum(state.d))
    return _diagnostic(problem, state, trial, 0, 0.0, 0.0, 0.0, 0.0; status = "initial")
end

function _failure_diagnostic(problem, state, err, step, time, cumulative_work, cumulative_viscous, cumulative_numerical)
    nonlinear = rlm_nonlinear_energy(problem, state.u, state.d)
    positive, negative = rlm_elastic_split_energies(problem, state.u, state.d)
    elastic = positive + negative
    fracture = 0.5 * dot(state.d, problem.K_AT2 * state.d)
    internal = elastic + fracture
    data = err.data
    load_data = _history_step_data(problem.config.load, time - problem.config.time.dt, time)
    return RLMDiagnostic(step = step, time = time, dt = problem.config.time.dt,
        displacement_factor = load_data.displacement_factor, body_force_factor = load_data.body_force_factor,
        traction_factor = load_data.traction_factor, displacement_rate = load_data.displacement_rate,
        body_force_rate = load_data.body_force_rate, traction_rate = load_data.traction_rate,
        displacement = problem.config.load.displacement_amplitude * load_data.displacement_factor,
        accepted = false, status = "rollback:$(err.code):$(err.message)", raw_energy = internal - dot(problem.f_ext, state.u),
        internal_energy = internal, proxy_energy = internal, elastic_energy = elastic,
        positive_elastic_energy = positive, negative_elastic_energy = negative, fracture_energy = fracture,
        nonlinear_energy = nonlinear, predicted_energy = _named_value(data, :P),
        prediction_gap = _named_value(data, :P) - nonlinear, proxy_gap = 0.0, q = state.q,
        q_minus_one = state.q - 1.0, c0 = _named_value(data, :c0), c1 = _named_value(data, :c1),
        A = _named_value(data, :A), B = _named_value(data, :B), C = _named_value(data, :C),
        discriminant = _named_value(data, :discriminant), discriminant_used = _named_value(data, :discriminant_used),
        scalar_residual = _named_value(data, :scalar_residual, Inf), cumulative_external_work = cumulative_work,
        cumulative_viscous_dissipation = cumulative_viscous, cumulative_numerical_dissipation = cumulative_numerical,
        min_d = minimum(state.d), max_d = maximum(state.d), energy_balance_residual = _named_value(data, :energy_balance_residual))
end

_csv_value(value) = value isa AbstractString ? "\"" * replace(value, "\"" => "\"\"") * "\"" : string(value)

function write_rlm_diagnostics(path::AbstractString, diagnostics::Vector{RLMDiagnostic})
    mkpath(dirname(path)); names = fieldnames(RLMDiagnostic)
    open(path, "w") do io
        println(io, join(string.(names), ','))
        for diagnostic in diagnostics
            println(io, join((_csv_value(getfield(diagnostic, name)) for name in names), ','))
        end
    end
    return path
end

function write_rlm_time_history(path::AbstractString, diagnostics::Vector{RLMDiagnostic})
    mkpath(dirname(path))
    names = (:step, :time, :dt, :displacement, :reaction_force, :internal_energy, :proxy_energy,
        :external_work, :cumulative_external_work, :viscous_dissipation, :numerical_dissipation,
        :cumulative_viscous_dissipation, :cumulative_numerical_dissipation, :q, :min_d, :max_d,
        :energy_balance_residual)
    open(path, "w") do io
        println(io, join(string.(names), ','))
        for diagnostic in diagnostics
            diagnostic.accepted || continue
            println(io, join((_csv_value(getfield(diagnostic, name)) for name in names), ','))
        end
    end
    return path
end

function _write_rlm_vtk(problem, state, step)
    output = problem.config.output; mkpath(output.directory)
    path = joinpath(output.directory, "time_$(lpad(step, 4, '0'))")
    VTKGridFile(path, problem.dh_u) do vtk
        write_solution(vtk, problem.dh_u, state.u); write_solution(vtk, problem.dh_d, state.d)
    end
    return path
end

function _flush_outputs(problem, diagnostics)
    problem.config.output.write_csv || return nothing
    write_rlm_diagnostics(joinpath(problem.config.output.directory, "diagnostics.csv"), diagnostics)
    write_rlm_time_history(joinpath(problem.config.output.directory, "time_history.csv"), diagnostics)
    return nothing
end

"""Advance Miehe--RLM-PE--BDF1 once per physical-time interval, without load relaxation."""
function solve_rlm_bdf1(problem::RLMProblem)
    config = problem.config
    update!(problem.ch_u, 0.0); update_rlm_external_force!(problem, 0.0)
    state = RLMState(zeros(ndofs(problem.dh_u)), fill(config.load.initial_damage, ndofs(problem.dh_d)), 1.0, 0.0)
    apply!(state.u, problem.ch_u); state.P = rlm_nonlinear_energy(problem, state.u, state.d)
    diagnostics = RLMDiagnostic[_initial_diagnostic(problem, state)]
    cumulative_work = 0.0; cumulative_viscous = 0.0; cumulative_numerical = 0.0
    nsteps = round(Int, config.time.final_time / config.time.dt)
    for step in 1:nsteps
        # Use the configured endpoint exactly; repeated floating-point
        # multiplication can otherwise produce e.g. 1.2000000000000002.
        time = step == nsteps ? config.time.final_time : step * config.time.dt
        local trial
        try
            trial = compute_rlm_bdf1_trial(problem, state, time)
        catch err
            if err isa RLMStepFailure
                push!(diagnostics, _failure_diagnostic(problem, state, err, step, time, cumulative_work, cumulative_viscous, cumulative_numerical))
                _flush_outputs(problem, diagnostics)
                config.output.verbose && @error sprint(showerror, err)
                return RLMResult(false, false, sprint(showerror, err) * "; candidate rolled back", state, diagnostics, problem)
            end
            rethrow()
        end
        cumulative_work += trial.external_work
        cumulative_viscous += trial.viscous_dissipation
        cumulative_numerical += trial.numerical_dissipation
        state = RLMState(trial.u, trial.d, trial.q, trial.P)
        diagnostic = _diagnostic(problem, state, trial, step, time, cumulative_work, cumulative_viscous, cumulative_numerical)
        push!(diagnostics, diagnostic)
        config.output.verbose && println("time $time: q-1=$(diagnostic.q_minus_one), reaction=$(diagnostic.reaction_force), Δd=$(diagnostic.phase_increment), energy residual=$(diagnostic.energy_balance_residual)")
        if config.output.write_vtk && (step % config.output.vtk_every_time_step == 0 || step == nsteps)
            _write_rlm_vtk(problem, state, step)
        end
    end
    _flush_outputs(problem, diagnostics)
    return RLMResult(true, true, "all real-time steps completed", state, diagnostics, problem)
end
