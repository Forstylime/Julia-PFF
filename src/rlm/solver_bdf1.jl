@inline function _named_value(data::NamedTuple, name::Symbol, default = NaN)
    return haskey(data, name) ? getfield(data, name) : default
end

function _throw_step_failure(code, message; kwargs...)
    throw(RLMStepFailure(code, message, (; kwargs...)))
end

function _solve_affine_branches(problem::RLMProblem, state::RLMState, n_u, n_d)
    inverse_M_dt = 1.0 / (
        problem.config.material.mobility * problem.config.time.dt
    )

    rhs_u_a = copy(problem.f_ext)
    apply_rhs!(problem.rhs_data_u, rhs_u_a, problem.ch_u, false)
    u_a = problem.factor_u \ rhs_u_a
    apply!(u_a, problem.ch_u)

    rhs_u_b = -n_u
    apply_rhs!(problem.rhs_data_u, rhs_u_b, problem.ch_u, true)
    u_b = problem.factor_u \ rhs_u_b
    apply_zero!(u_b, problem.ch_u)

    rhs_d_a = inverse_M_dt * (problem.M_d * state.d)
    d_a = problem.factor_d \ rhs_d_a
    d_b = problem.factor_d \ (-n_d)
    return u_a, u_b, d_a, d_b
end

"""
    compute_rlm_bdf1_trial(problem, state)

Compute one complete BDF1 candidate without mutating `state`. Any failed check
throws `RLMStepFailure`; callers can therefore roll back by simply not committing
the returned trial.
"""
function compute_rlm_bdf1_trial(problem::RLMProblem, state::RLMState)
    config = problem.config
    tolerances = config.tolerances
    alpha = config.time.alpha
    n_u = zeros(ndofs(problem.dh_u))
    n_d = zeros(ndofs(problem.dh_d))
    assemble_rlm_nonlinear_forces!(n_u, n_d, problem, state.u, state.d)
    u_a, u_b, d_a, d_b = _solve_affine_branches(problem, state, n_u, n_d)

    u_star = u_a + u_b
    d_star = d_a + d_b
    P_next = rlm_nonlinear_energy(problem, u_star, d_star)
    c0 = dot(n_u, u_a - state.u) + dot(n_d, d_a - state.d)
    c1 = dot(n_u, u_b) + dot(n_d, d_b)

    c1_from_branches = -dot(u_b, problem.K_u * u_b) -
                       dot(d_b, problem.K_d * d_b)
    branch_scale = abs(c1) + abs(c1_from_branches) + eps(Float64)
    branch_residual = abs(c1 - c1_from_branches) / branch_scale
    if branch_residual > tolerances.branch_identity
        _throw_step_failure(
            :c1_branch_identity,
            "c1 does not match the negative affine-branch quadratic form";
            P = P_next,
            c0,
            c1,
            c1_from_branches,
            branch_residual,
        )
    end
    c1_tolerance = tolerances.c1_abs +
                   tolerances.c1_rel * max(abs(c1), abs(c1_from_branches))
    if c1 > c1_tolerance
        _throw_step_failure(
            :positive_c1,
            "c1=$c1 is significantly positive (tolerance=$c1_tolerance)";
            P = P_next,
            c0,
            c1,
        )
    end

    # These are exactly the coefficients in section 12 of the source document.
    A = alpha - c1
    B = -c0
    C = P_next - state.P - alpha * state.q^2
    root = solve_rlm_quadratic(A, B, C, state.q, tolerances)
    if !root.success
        _throw_step_failure(
            root.code,
            root.message;
            P = P_next,
            c0,
            c1,
            A,
            B,
            C,
            discriminant = root.discriminant,
            discriminant_used = root.discriminant_used,
            scalar_residual = root.residual,
        )
    end

    q_next = root.q
    u_next = u_a + q_next * u_b
    d_next = d_a + q_next * d_b
    # Reassert the analytically known branch boundary values to remove solve roundoff.
    apply!(u_next, problem.ch_u)

    raw_energy = rlm_raw_energy(problem, u_next, d_next)
    proxy_energy = rlm_proxy_energy(problem, u_next, d_next, q_next, P_next)
    proxy_old = rlm_proxy_energy(problem, state.u, state.d, state.q, state.P)
    delta_u = u_next - state.u
    delta_d = d_next - state.d
    inverse_M_dt = 1.0 / (config.material.mobility * config.time.dt)
    dissipation = 0.5 * dot(delta_u, problem.K_u * delta_u) +
                  0.5 * dot(delta_d, problem.K_AT2 * delta_d) +
                  inverse_M_dt * dot(delta_d, problem.M_d * delta_d)
    balance_raw = proxy_energy - proxy_old + dissipation
    balance_scale = abs(proxy_energy) + abs(proxy_old) + abs(dissipation) + eps(Float64)
    balance_residual = abs(balance_raw) / balance_scale
    balance_tolerance = tolerances.energy_balance_abs +
                        tolerances.energy_balance_rel * balance_scale
    if abs(balance_raw) > balance_tolerance
        _throw_step_failure(
            :energy_balance,
            "BDF1 proxy-energy identity residual $(abs(balance_raw)) exceeds $balance_tolerance";
            P = P_next,
            c0,
            c1,
            A,
            B,
            C,
            discriminant = root.discriminant,
            discriminant_used = root.discriminant_used,
            scalar_residual = root.residual,
            energy_balance_residual = balance_residual,
        )
    end

    phase_increment, phase_relative_increment, healing, min_d, max_d =
        phase_field_metrics(problem, d_next, state.d)
    return RLMTrial(
        u_a,
        u_b,
        d_a,
        d_b,
        u_star,
        d_star,
        u_next,
        d_next,
        n_u,
        n_d,
        P_next,
        c0,
        c1,
        A,
        B,
        C,
        root.discriminant,
        root.discriminant_used,
        q_next,
        root.residual,
        raw_energy,
        proxy_energy,
        balance_residual,
        phase_increment,
        phase_relative_increment,
        healing,
        min_d,
        max_d,
    )
end

function _accepted_diagnostic(problem, trial, load_step, relax_step, load_fraction)
    nonlinear_actual = rlm_nonlinear_energy(problem, trial.u, trial.d)
    return RLMDiagnostic(
        load_step = load_step,
        relax_step = relax_step,
        load_fraction = load_fraction,
        displacement = load_fraction * problem.config.load.final_displacement,
        accepted = true,
        status = "accepted",
        raw_energy = trial.raw_energy,
        proxy_energy = trial.proxy_energy,
        nonlinear_energy = nonlinear_actual,
        predicted_energy = trial.P,
        prediction_gap = trial.P - nonlinear_actual,
        proxy_gap = abs(trial.proxy_energy - trial.raw_energy),
        q = trial.q,
        q_minus_one = trial.q - 1.0,
        c0 = trial.c0,
        c1 = trial.c1,
        A = trial.A,
        B = trial.B,
        C = trial.C,
        discriminant = trial.discriminant,
        discriminant_used = trial.discriminant_used,
        scalar_residual = trial.scalar_residual,
        phase_increment = trial.phase_increment,
        phase_relative_increment = trial.phase_relative_increment,
        healing = trial.healing,
        min_d = trial.min_d,
        max_d = trial.max_d,
        energy_balance_residual = trial.energy_balance_residual,
    )
end

function _reset_diagnostic(problem, state, load_step, load_fraction, reset_jump)
    nonlinear = rlm_nonlinear_energy(problem, state.u, state.d)
    raw = rlm_quadratic_energy(problem, state.u, state.d) + nonlinear
    proxy = rlm_proxy_energy(problem, state.u, state.d, state.q, state.P)
    _, _, _, min_d, max_d = phase_field_metrics(problem, state.d, state.d)
    return RLMDiagnostic(
        load_step = load_step,
        relax_step = 0,
        load_fraction = load_fraction,
        displacement = load_fraction * problem.config.load.final_displacement,
        accepted = true,
        status = load_step == 0 ? "initial" : "load_reset",
        raw_energy = raw,
        proxy_energy = proxy,
        nonlinear_energy = nonlinear,
        predicted_energy = state.P,
        prediction_gap = state.P - nonlinear,
        proxy_gap = abs(proxy - raw),
        q = state.q,
        q_minus_one = state.q - 1.0,
        phase_increment = 0.0,
        phase_relative_increment = 0.0,
        healing = 0.0,
        min_d = min_d,
        max_d = max_d,
        reset_jump = reset_jump,
    )
end

function _failure_diagnostic(problem, state, err, load_step, relax_step, load_fraction)
    raw = rlm_raw_energy(problem, state.u, state.d)
    proxy = rlm_proxy_energy(problem, state.u, state.d, state.q, state.P)
    nonlinear = rlm_nonlinear_energy(problem, state.u, state.d)
    _, _, _, min_d, max_d = phase_field_metrics(problem, state.d, state.d)
    data = err.data
    return RLMDiagnostic(
        load_step = load_step,
        relax_step = relax_step,
        load_fraction = load_fraction,
        displacement = load_fraction * problem.config.load.final_displacement,
        accepted = false,
        status = "rollback:$(err.code):$(err.message)",
        raw_energy = raw,
        proxy_energy = proxy,
        nonlinear_energy = nonlinear,
        predicted_energy = _named_value(data, :P),
        prediction_gap = _named_value(data, :P) - nonlinear,
        proxy_gap = abs(proxy - raw),
        q = state.q,
        q_minus_one = state.q - 1.0,
        c0 = _named_value(data, :c0),
        c1 = _named_value(data, :c1),
        A = _named_value(data, :A),
        B = _named_value(data, :B),
        C = _named_value(data, :C),
        discriminant = _named_value(data, :discriminant),
        discriminant_used = _named_value(data, :discriminant_used),
        scalar_residual = _named_value(data, :scalar_residual, Inf),
        phase_increment = 0.0,
        phase_relative_increment = 0.0,
        healing = 0.0,
        min_d = min_d,
        max_d = max_d,
        energy_balance_residual = _named_value(data, :energy_balance_residual),
    )
end

function _csv_value(value)
    if value isa AbstractString
        return "\"" * replace(value, "\"" => "\"\"") * "\""
    end
    return string(value)
end

function write_rlm_diagnostics(path::AbstractString, diagnostics::Vector{RLMDiagnostic})
    mkpath(dirname(path))
    names = fieldnames(RLMDiagnostic)
    open(path, "w") do io
        println(io, join(string.(names), ','))
        for diagnostic in diagnostics
            values = (_csv_value(getfield(diagnostic, name)) for name in names)
            println(io, join(values, ','))
        end
    end
    return path
end

function _write_rlm_vtk(problem, state, load_step)
    output = problem.config.output
    mkpath(output.directory)
    path = joinpath(output.directory, "load_$(lpad(load_step, 4, '0'))")
    VTKGridFile(path, problem.dh_u) do vtk
        write_solution(vtk, problem.dh_u, state.u)
        write_solution(vtk, problem.dh_d, state.d)
    end
    return path
end

function _flush_outputs(problem, diagnostics)
    output = problem.config.output
    output.write_csv && write_rlm_diagnostics(
        joinpath(output.directory, "diagnostics.csv"), diagnostics,
    )
    return nothing
end

"""
    solve_rlm_bdf1(problem)

Run outer displacement increments and inner fixed-load BDF1 relaxation. Load
resets follow section 19. A failed trial is reported and returned without
changing the last accepted `RLMState`.
"""
function solve_rlm_bdf1(problem::RLMProblem)
    config = problem.config
    output = config.output
    state = RLMState(
        zeros(ndofs(problem.dh_u)),
        fill(config.load.initial_damage, ndofs(problem.dh_d)),
        1.0,
        0.0,
    )
    apply!(state.u, problem.ch_u)
    state.P = rlm_nonlinear_energy(problem, state.u, state.d)
    diagnostics = RLMDiagnostic[_reset_diagnostic(problem, state, 0, 0.0, 0.0)]

    for load_step in 1:config.load.load_steps
        load_fraction = load_step / config.load.load_steps
        update!(problem.ch_u, load_fraction)

        # Construct U^{k,0} in the new affine space, then isolate the auxiliary
        # reset jump at fixed new physical fields.
        u_reset = copy(state.u)
        apply!(u_reset, problem.ch_u)
        proxy_before_reset = rlm_quadratic_energy(problem, u_reset, state.d) +
                             state.P + config.time.alpha * (state.q^2 - 1.0)
        P_reset = rlm_nonlinear_energy(problem, u_reset, state.d)
        state = RLMState(u_reset, copy(state.d), 1.0, P_reset)
        proxy_after_reset = rlm_proxy_energy(problem, state.u, state.d, state.q, state.P)
        reset_jump = proxy_after_reset - proxy_before_reset
        push!(diagnostics, _reset_diagnostic(
            problem, state, load_step, load_fraction, reset_jump,
        ))

        output.verbose && println(
            "load $load_step/$(config.load.load_steps), displacement=",
            load_fraction * config.load.final_displacement,
            ", reset_jump=", reset_jump,
        )

        converged_this_load = false
        for relax_step in 1:config.time.max_relax_steps
            local trial
            try
                trial = compute_rlm_bdf1_trial(problem, state)
            catch err
                if err isa RLMStepFailure
                    push!(diagnostics, _failure_diagnostic(
                        problem, state, err, load_step, relax_step, load_fraction,
                    ))
                    _flush_outputs(problem, diagnostics)
                    output.verbose && @error sprint(showerror, err)
                    return RLMResult(
                        false,
                        false,
                        sprint(showerror, err) * "; candidate rolled back",
                        state,
                        diagnostics,
                        problem,
                    )
                end
                rethrow()
            end

            # The only mutation point for an accepted BDF1 step.
            state = RLMState(trial.u, trial.d, trial.q, trial.P)
            diagnostic = _accepted_diagnostic(
                problem, trial, load_step, relax_step, load_fraction,
            )
            push!(diagnostics, diagnostic)
            output.verbose && println(
                "  relax $relax_step: q-1=$(diagnostic.q_minus_one), " *
                "c1=$(diagnostic.c1), D=$(diagnostic.discriminant), " *
                "rq=$(diagnostic.scalar_residual), Δd=$(diagnostic.phase_increment), " *
                "heal=$(diagnostic.healing)",
            )

            if relax_step >= config.time.min_relax_steps &&
               trial.phase_relative_increment < config.tolerances.phase &&
               abs(trial.q - 1.0) < config.tolerances.q
                converged_this_load = true
                break
            end
        end

        if !converged_this_load
            message = "load step $load_step did not converge in " *
                      "$(config.time.max_relax_steps) fixed-load relaxation steps"
            _flush_outputs(problem, diagnostics)
            output.verbose && @error message
            return RLMResult(false, false, message, state, diagnostics, problem)
        end

        if output.write_vtk &&
           (load_step % output.vtk_every_load_step == 0 ||
            load_step == config.load.load_steps)
            _write_rlm_vtk(problem, state, load_step)
        end
        _flush_outputs(problem, diagnostics)
    end

    message = "all load steps converged"
    return RLMResult(true, true, message, state, diagnostics, problem)
end
