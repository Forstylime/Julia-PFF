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
    u_old = zeros(ndofs(problem.dh_u)); update!(problem.ch_u, t_old); apply!(u_old, problem.ch_u)
    u_new = zeros(ndofs(problem.dh_u)); update!(problem.ch_u, t_new); apply!(u_new, problem.ch_u)
    return u_new - u_old
end

@inline _affine_state(u_a, u_b, d_a, d_b, q) = (u_a + q * u_b, d_a + q * d_b)

function _qm_initial_curvature(problem, u_a, u_b, d_a, d_b, phi_star)
    qm = problem.config.qm
    h = min(qm.finite_difference_step, 1.0 - qm.q_min, qm.q_max - 1.0)
    estimated = 0.0
    if h > 0.0
        u_minus, d_minus = _affine_state(u_a, u_b, d_a, d_b, 1.0 - h)
        u_plus, d_plus = _affine_state(u_a, u_b, d_a, d_b, 1.0 + h)
        estimated = max(0.0, (rlm_nonlinear_energy(problem, u_plus, d_plus) -
            2.0 * phi_star + rlm_nonlinear_energy(problem, u_minus, d_minus)) / h^2)
    end
    return max(qm.initial_curvature, estimated)
end

@inline function _qm_majorant(phi_star, g_star, curvature, q)
    return phi_star + g_star * (q - 1.0) + 0.5 * curvature * (q - 1.0)^2
end

function _next_curvature(curvature, scale, qm)
    seed = max(1.0e-12, sqrt(eps(Float64)) * max(1.0, scale))
    return curvature == 0.0 ? seed : curvature * qm.curvature_growth
end

"""Return QM-admissible roots, their reconstructed fields, and majorant margins."""
function _qm_admissible_candidates(problem, root, A, B, C, u_a, u_b, d_a, d_b,
    phi_star, g_star, curvature)
    qm, tol = problem.config.qm, problem.config.tolerances
    accepted = NamedTuple[]
    for q in root.candidates
        isfinite(q) && q > tol.positive_root || continue
        q >= qm.q_min && q <= qm.q_max || continue
        residual = scalar_equation_residual(A, B, C, q, tol.scalar_denominator_epsilon)
        residual <= tol.scalar_residual || continue
        u, d = _affine_state(u_a, u_b, d_a, d_b, q)
        phi = rlm_nonlinear_energy(problem, u, d)
        majorant = _qm_majorant(phi_star, g_star, curvature, q)
        margin = majorant - phi
        margin_tolerance = qm.majorant_abs + qm.majorant_rel * max(abs(majorant), abs(phi))
        margin >= -margin_tolerance || continue
        push!(accepted, (; q, u, d, phi, margin, residual))
    end
    return accepted
end

"""Compute one transactional BDF1--QM candidate without mutating `state`."""
function compute_rlm_bdf1_trial(problem::RLMProblem, state::RLMState, t_next::Real)
    config = problem.config
    t_new = Float64(t_next); t_old = t_new - config.time.dt
    _history_step_data(config.load, t_old, t_new)
    update!(problem.ch_u, t_new); update_rlm_external_force!(problem, t_new)

    tolerances, alpha = config.tolerances, config.time.alpha
    n_u = zeros(ndofs(problem.dh_u)); n_d = zeros(ndofs(problem.dh_d))
    assemble_rlm_nonlinear_forces!(n_u, n_d, problem, state.u, state.d)
    u_a = _solve_u_baseline(problem)
    u_a, u_b, d_a, d_b = _solve_affine_branches(problem, state, n_u, n_d, u_a)
    u_star, d_star = _affine_state(u_a, u_b, d_a, d_b, 1.0)
    phi_star = rlm_nonlinear_energy(problem, u_star, d_star)
    n_u_star = zeros(ndofs(problem.dh_u)); n_d_star = zeros(ndofs(problem.dh_d))
    assemble_rlm_nonlinear_forces!(n_u_star, n_d_star, problem, u_star, d_star)
    g_star = dot(n_u_star, u_b) + dot(n_d_star, d_b)
    phi_old = rlm_nonlinear_energy(problem, state.u, state.d)
    c0 = dot(n_u, u_a - state.u) + dot(n_d, d_a - state.d)
    c1 = dot(n_u, u_b) + dot(n_d, d_b)
    c1_from_branches = -dot(u_b, problem.K_u * u_b) - dot(d_b, problem.K_d * d_b)
    branch_scale = abs(c1) + abs(c1_from_branches) + eps(Float64)
    abs(c1 - c1_from_branches) / branch_scale <= tolerances.branch_identity || _throw_step_failure(
        :c1_branch_identity, "c1 does not match the negative affine-branch quadratic form";
        c0, c1, c1_from_branches,
    )
    c1_tolerance = tolerances.c1_abs + tolerances.c1_rel * max(abs(c1), abs(c1_from_branches))
    c1 <= c1_tolerance || _throw_step_failure(:positive_c1, "c1 is significantly positive"; c0, c1)

    curvature = _qm_initial_curvature(problem, u_a, u_b, d_a, d_b, phi_star)
    selected = nothing
    A = B = C = discriminant = discriminant_used = NaN
    for attempt in 0:config.qm.max_backtracks
        A = alpha + 0.5 * curvature - c1
        B = g_star - curvature - c0
        C = phi_star - g_star + 0.5 * curvature - phi_old - alpha * state.q^2
        root = solve_rlm_quadratic(A, B, C, state.q, tolerances)
        discriminant, discriminant_used = root.discriminant, root.discriminant_used
        candidates = _qm_admissible_candidates(problem, root, A, B, C, u_a, u_b, d_a, d_b,
            phi_star, g_star, curvature)
        if !isempty(candidates)
            selected = candidates[argmin(abs(candidate.q - state.q) for candidate in candidates)]
            break
        end
        attempt == config.qm.max_backtracks ||
            (curvature = _next_curvature(curvature, abs(phi_star) + abs(g_star) + abs(phi_old), config.qm))
    end
    selected === nothing && _throw_step_failure(
        :qm_no_admissible_root,
        "no scalar root satisfied the configured q interval, residual, and QM majorant after $(config.qm.max_backtracks) backtracks";
        c0, c1, A, B, C, phi_star, g_star, curvature, discriminant, discriminant_used,
    )

    q_next, u_next, d_next = selected.q, selected.u, selected.d
    apply!(u_next, problem.ch_u)
    Ku_u, Kat2_d = problem.K_u * u_next, problem.K_AT2 * d_next
    delta_u, delta_d = u_next - state.u, d_next - state.d
    viscosity_over_dt = config.material.viscosity / config.time.dt
    nonlinear_actual = selected.phi
    positive_elastic_energy, negative_elastic_energy = rlm_elastic_split_energies(problem, u_next, d_next)
    elastic_energy = positive_elastic_energy + negative_elastic_energy
    fracture_energy = 0.5 * dot(d_next, Kat2_d)
    internal_energy = elastic_energy + fracture_energy
    raw_energy = internal_energy - dot(problem.f_ext, u_next)
    relaxed_internal_energy = internal_energy + alpha * (q_next^2 - 1.0)
    relaxed_old = rlm_relaxed_internal_energy(problem, state.u, state.d, state.q)
    viscous_dissipation = viscosity_over_dt * dot(delta_d, problem.M_d * delta_d)
    numerical_dissipation = 0.5 * dot(delta_u, problem.K_u * delta_u) + 0.5 * dot(delta_d, problem.K_AT2 * delta_d)
    mechanical_residual = Ku_u + q_next * n_u - problem.f_ext
    external_work = dot(problem.f_ext, delta_u) + dot(mechanical_residual, _dirichlet_increment(problem, t_old, t_new))
    inequality_raw = relaxed_internal_energy - relaxed_old - external_work + viscous_dissipation + numerical_dissipation
    inequality_scale = abs(relaxed_internal_energy) + abs(relaxed_old) + abs(external_work) +
        abs(viscous_dissipation) + abs(numerical_dissipation) + eps(Float64)
    inequality_violation = max(0.0, inequality_raw) / inequality_scale
    inequality_tolerance = tolerances.energy_balance_abs + tolerances.energy_balance_rel * inequality_scale
    max(0.0, inequality_raw) <= inequality_tolerance || _throw_step_failure(
        :energy_inequality, "BDF1-QM work-energy inequality violation exceeds tolerance";
        c0, c1, A, B, C, phi_star, g_star, curvature, energy_inequality_violation = inequality_violation,
    )
    phase_residual_vector = viscosity_over_dt * (problem.M_d * delta_d) + Kat2_d + q_next * n_d
    phase_scale = viscosity_over_dt * norm(problem.M_d * delta_d) + norm(Kat2_d) + abs(q_next) * norm(n_d) + eps(Float64)
    phase_equilibrium_residual = norm(phase_residual_vector) / phase_scale
    # The frozen RLM equilibrium and the original Miehe constitutive residual
    # coincide only when q=1 and the frozen/nonlinear states coincide.
    reaction_rlm = sum(mechanical_residual[problem.loaded_component_dofs])
    n_u_phys = zeros(ndofs(problem.dh_u)); n_d_phys = zeros(ndofs(problem.dh_d))
    assemble_rlm_nonlinear_forces!(n_u_phys, n_d_phys, problem, u_next, d_next)
    physical_residual = Ku_u + n_u_phys - problem.f_ext
    reaction_phys = sum(physical_residual[problem.loaded_component_dofs])
    phase_increment, phase_relative_increment, healing, min_d, max_d = phase_field_metrics(problem, d_next, state.d)
    return RLMTrial(u_a, u_b, d_a, d_b, u_star, d_star, u_next, d_next, n_u, n_d,
        phi_star, g_star, curvature, selected.margin, c0, c1, A, B, C, discriminant, discriminant_used,
        q_next, selected.residual, raw_energy, relaxed_internal_energy, elastic_energy, positive_elastic_energy,
        negative_elastic_energy, fracture_energy, nonlinear_actual, reaction_rlm, reaction_phys, internal_energy, external_work,
        viscous_dissipation, numerical_dissipation, inequality_violation, phase_increment, phase_relative_increment,
        phase_equilibrium_residual, healing, min_d, max_d)
end

function _diagnostic(problem, trial, step::Int, time::Float64, cumulative_work, cumulative_viscous, cumulative_numerical; status = "accepted", accepted = true)
    load = problem.config.load
    data = step == 0 ? (; displacement_factor = history_value(load.displacement_history, 0.0), body_force_factor = history_value(load.body_force_history, 0.0), traction_factor = history_value(load.traction_history, 0.0), displacement_rate = 0.0, body_force_rate = 0.0, traction_rate = 0.0) : _history_step_data(load, time - problem.config.time.dt, time)
    return RLMDiagnostic(step = step, time = time, dt = step == 0 ? 0.0 : problem.config.time.dt,
        displacement_factor = data.displacement_factor, body_force_factor = data.body_force_factor, traction_factor = data.traction_factor,
        displacement_rate = data.displacement_rate, body_force_rate = data.body_force_rate, traction_rate = data.traction_rate,
        displacement = load.displacement_amplitude * data.displacement_factor, accepted = accepted, status = status,
        raw_energy = trial.raw_energy, internal_energy = trial.internal_energy, relaxed_internal_energy = trial.relaxed_internal_energy,
        elastic_energy = trial.elastic_energy, positive_elastic_energy = trial.positive_elastic_energy, negative_elastic_energy = trial.negative_elastic_energy,
        fracture_energy = trial.fracture_energy, nonlinear_energy = trial.nonlinear_energy, phi_star = trial.phi_star, g_star = trial.g_star,
        curvature = trial.curvature, majorant_margin = trial.majorant_margin, q = trial.q, q_minus_one = trial.q - 1.0,
        c0 = trial.c0, c1 = trial.c1, A = trial.A, B = trial.B, C = trial.C, discriminant = trial.discriminant,
        discriminant_used = trial.discriminant_used, scalar_residual = trial.scalar_residual,
        reaction_rlm = trial.reaction_rlm, reaction_phys = trial.reaction_phys,
        external_work = trial.external_work, cumulative_external_work = cumulative_work, viscous_dissipation = trial.viscous_dissipation,
        numerical_dissipation = trial.numerical_dissipation, cumulative_viscous_dissipation = cumulative_viscous,
        cumulative_numerical_dissipation = cumulative_numerical, phase_increment = trial.phase_increment,
        phase_relative_increment = trial.phase_relative_increment, phase_equilibrium_residual = trial.phase_equilibrium_residual,
        healing = trial.healing, min_d = trial.min_d, max_d = trial.max_d, energy_inequality_violation = trial.energy_inequality_violation)
end

function _initial_diagnostic(problem, state)
    positive, negative = rlm_elastic_split_energies(problem, state.u, state.d)
    elastic, fracture = positive + negative, 0.5 * dot(state.d, problem.K_AT2 * state.d)
    internal = elastic + fracture
    trial = RLMTrial(state.u, zeros(length(state.u)), state.d, zeros(length(state.d)), state.u, state.d, state.u, state.d,
        zeros(length(state.u)), zeros(length(state.d)), NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN,
        state.q, NaN, internal, internal + problem.config.time.alpha * (state.q^2 - 1.0), elastic, positive, negative,
        fracture, rlm_nonlinear_energy(problem, state.u, state.d), 0.0, 0.0, internal, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        minimum(state.d), maximum(state.d))
    return _diagnostic(problem, trial, 0, 0.0, 0.0, 0.0, 0.0; status = "initial")
end

function _failure_diagnostic(problem, state, err, step, time, cumulative_work, cumulative_viscous, cumulative_numerical)
    positive, negative = rlm_elastic_split_energies(problem, state.u, state.d)
    elastic, fracture = positive + negative, 0.5 * dot(state.d, problem.K_AT2 * state.d)
    internal = elastic + fracture; data = err.data
    load_data = _history_step_data(problem.config.load, time - problem.config.time.dt, time)
    return RLMDiagnostic(step = step, time = time, dt = problem.config.time.dt,
        displacement_factor = load_data.displacement_factor, body_force_factor = load_data.body_force_factor,
        traction_factor = load_data.traction_factor, displacement_rate = load_data.displacement_rate,
        body_force_rate = load_data.body_force_rate, traction_rate = load_data.traction_rate,
        displacement = problem.config.load.displacement_amplitude * load_data.displacement_factor, accepted = false,
        status = "rollback:$(err.code):$(err.message)", raw_energy = internal - dot(problem.f_ext, state.u), internal_energy = internal,
        relaxed_internal_energy = rlm_relaxed_internal_energy(problem, state.u, state.d, state.q), elastic_energy = elastic,
        positive_elastic_energy = positive, negative_elastic_energy = negative, fracture_energy = fracture,
        nonlinear_energy = rlm_nonlinear_energy(problem, state.u, state.d), q = state.q, q_minus_one = state.q - 1.0,
        c0 = _named_value(data, :c0), c1 = _named_value(data, :c1), A = _named_value(data, :A), B = _named_value(data, :B),
        C = _named_value(data, :C), phi_star = _named_value(data, :phi_star), g_star = _named_value(data, :g_star),
        curvature = _named_value(data, :curvature), discriminant = _named_value(data, :discriminant),
        discriminant_used = _named_value(data, :discriminant_used), energy_inequality_violation = _named_value(data, :energy_inequality_violation),
        cumulative_external_work = cumulative_work, cumulative_viscous_dissipation = cumulative_viscous,
        cumulative_numerical_dissipation = cumulative_numerical, min_d = minimum(state.d), max_d = maximum(state.d))
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
    names = (:step, :time, :dt, :displacement, :reaction_rlm, :reaction_phys, :internal_energy, :relaxed_internal_energy,
        :external_work, :cumulative_external_work, :viscous_dissipation, :numerical_dissipation,
        :cumulative_viscous_dissipation, :cumulative_numerical_dissipation, :q, :curvature, :majorant_margin,
        :min_d, :max_d, :energy_inequality_violation)
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

function _format_duration(seconds::Real)
    total = max(0, round(Int, seconds))
    hours, remainder = divrem(total, 3600)
    minutes, seconds_part = divrem(remainder, 60)
    return hours > 0 ? "$(hours)h $(lpad(minutes, 2, '0'))m" :
           minutes > 0 ? "$(minutes)m $(lpad(seconds_part, 2, '0'))s" : "$(seconds_part)s"
end

function _write_rlm_progress(io, step, nsteps, started_ns; final = false)
    elapsed = (time_ns() - started_ns) / 1.0e9
    fraction = step / nsteps
    width = 24
    filled = clamp(round(Int, width * fraction), 0, width)
    bar = repeat("█", filled) * repeat("░", width - filled)
    eta = step == 0 ? "--" : _format_duration(elapsed * (nsteps - step) / step)
    print(io, "\rRLM-BDF1-QM [$bar] $(lpad(round(Int, 100 * fraction), 3))% " *
        "($step/$nsteps) elapsed $(_format_duration(elapsed)), ETA $eta")
    final && println(io)
    flush(io)
end

"""Advance Miehe--RLM--BDF1-QM once per physical-time interval."""
function solve_rlm_bdf1(problem::RLMProblem)
    config = problem.config
    update!(problem.ch_u, 0.0); update_rlm_external_force!(problem, 0.0)
    state = RLMState(zeros(ndofs(problem.dh_u)), fill(config.load.initial_damage, ndofs(problem.dh_d)), 1.0)
    apply!(state.u, problem.ch_u)
    diagnostics = RLMDiagnostic[_initial_diagnostic(problem, state)]
    cumulative_work = 0.0; cumulative_viscous = 0.0; cumulative_numerical = 0.0
    nsteps = round(Int, config.time.final_time / config.time.dt)
    progress_started_ns = time_ns()
    progress_last_ns = progress_started_ns
    config.output.show_progress && _write_rlm_progress(stdout, 0, nsteps, progress_started_ns)
    for step in 1:nsteps
        time = step == nsteps ? config.time.final_time : step * config.time.dt
        local trial
        try
            trial = compute_rlm_bdf1_trial(problem, state, time)
        catch err
            if err isa RLMStepFailure
                push!(diagnostics, _failure_diagnostic(problem, state, err, step, time, cumulative_work, cumulative_viscous, cumulative_numerical))
                _flush_outputs(problem, diagnostics)
                config.output.show_progress && println(stdout)
                config.output.verbose && @error sprint(showerror, err)
                return RLMResult(false, false, sprint(showerror, err) * "; candidate rolled back", state, diagnostics, problem)
            end
            rethrow()
        end
        cumulative_work += trial.external_work; cumulative_viscous += trial.viscous_dissipation; cumulative_numerical += trial.numerical_dissipation
        state = RLMState(trial.u, trial.d, trial.q)
        diagnostic = _diagnostic(problem, trial, step, time, cumulative_work, cumulative_viscous, cumulative_numerical)
        push!(diagnostics, diagnostic)
        config.output.verbose && println("time $time: q-1=$(diagnostic.q_minus_one), reaction_RLM=$(diagnostic.reaction_rlm), reaction_phys=$(diagnostic.reaction_phys), Δd=$(diagnostic.phase_increment), QM violation=$(diagnostic.energy_inequality_violation)")
        if config.output.write_vtk && (step % config.output.vtk_every_time_step == 0 || step == nsteps)
            _write_rlm_vtk(problem, state, step)
        end
        now_ns = time_ns()
        if config.output.show_progress && (step == nsteps ||
            (now_ns - progress_last_ns) / 1.0e9 >= config.output.progress_refresh_seconds)
            _write_rlm_progress(stdout, step, nsteps, progress_started_ns; final = step == nsteps)
            progress_last_ns = now_ns
        end
    end
    _flush_outputs(problem, diagnostics)
    return RLMResult(true, true, "all real-time steps completed", state, diagnostics, problem)
end
