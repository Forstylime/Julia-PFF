using Ferrite
using LinearAlgebra

const DEFAULT_STAGGERED_OUTPUT_INTERVAL = 100
const DEFAULT_MAX_NEWTON_ITERATIONS = 20

function validate_config(config::StaggeredConfig)
    time = config.time
    solver = config.solver
    output = config.output
    isfinite(time.total_time) || throw(ArgumentError("time.total_time must be finite"))
    (isnothing(time.dt) || isfinite(time.dt)) || throw(ArgumentError("time.dt must be finite"))
    all(isfinite, (solver.tolerance, solver.viscosity)) ||
        throw(ArgumentError("solver tolerance and viscosity must be finite"))
    time.n_steps > 0 || throw(ArgumentError("time.n_steps must be positive"))
    solver.tolerance > 0 || throw(ArgumentError("solver.tolerance must be positive"))
    solver.max_staggered_iterations > 0 ||
        throw(ArgumentError("solver.max_staggered_iterations must be positive"))
    solver.max_newton_iterations > 0 ||
        throw(ArgumentError("solver.max_newton_iterations must be positive"))
    solver.viscosity >= 0 || throw(ArgumentError("solver.viscosity must be non-negative"))
    output.vtk_interval > 0 || throw(ArgumentError("output.vtk_interval must be positive"))
    _staggered_time_grid(time.n_steps, time.total_time, time.time_points, time.dt)
    return config
end

function _staggered_time_grid(
    n_steps::Integer,
    total_time::Real,
    time_points::Union{Nothing,AbstractVector},
    dt::Union{Nothing,Real},
)
    if !isnothing(dt)
        isnothing(time_points) ||
            throw(ArgumentError("dt and time_points cannot be used together"))
        dt > 0 || throw(ArgumentError("dt must be positive"))
        return collect(0.0:float(dt):(n_steps * float(dt)))
    end
    if isnothing(time_points)
        total_time > 0 || throw(ArgumentError("total_time must be positive"))
        return collect(range(0.0, float(total_time); length = n_steps + 1))
    end
    length(time_points) == n_steps + 1 ||
        throw(ArgumentError("time_points must have n_steps + 1 entries"))
    times = Float64.(time_points)
    all(isfinite, times) || throw(ArgumentError("time_points must be finite"))
    all(diff(times) .> 0) || throw(ArgumentError("time_points must be strictly increasing"))
    return times
end

"""
    _staggered_load_function(times, load_history, ramp_time)

Build the displacement load-factor history.  `ramp_time` selects a two-stage
history: a linear ramp from zero to one followed by a constant hold.  A custom
`load_history` and `ramp_time` are mutually exclusive so that the prescribed
boundary condition is unambiguous.
"""
function _staggered_load_function(
    times::AbstractVector{<:Real},
    load_history::Union{Nothing,Function},
    ramp_time::Union{Nothing,Real},
)
    if !isnothing(ramp_time)
        isnothing(load_history) ||
            throw(ArgumentError("ramp_time and load_history cannot be used together"))
        isfinite(ramp_time) || throw(ArgumentError("ramp_time must be finite"))
        0.0 < ramp_time < times[end] ||
            throw(ArgumentError("ramp_time must lie strictly between zero and the final time"))
        return t -> min(float(t) / float(ramp_time), 1.0)
    end
    return isnothing(load_history) ? (t -> t / times[end]) : load_history
end

_staggered_load_function(times, load::StaggeredLoadConfig) =
    _staggered_load_function(times, load.load_history, load.ramp_time)

function _staggered_output_directory(config::StaggeredConfig)
    output = config.output
    !isnothing(output.directory) && return output.directory
    suffix = config.solver.viscosity > 0 ? "staggered_bdf1" :
             (config.solver.enforce_irreversibility ? "staggered" : "staggered2")
    return joinpath("data", "sims", suffix)
end

function _solve_displacement!(
    stiffness,
    residual,
    setup::TensionSetup,
    displacement,
    damage,
    material::MaterialParams,
    displacement_values,
    damage_values;
    tolerance::Real,
    max_iterations::Integer,
)
    assemble_u!(
        stiffness,
        residual,
        setup.dh_u,
        setup.dh_d,
        displacement,
        damage,
        material,
        displacement_values,
        damage_values,
    )
    apply_zero!(stiffness, residual, setup.ch_u)

    residual_norm = norm(residual)
    iteration = 0
    while residual_norm > tolerance && iteration < max_iterations
        increment = stiffness \ (-residual)
        apply_zero!(increment, setup.ch_u)
        displacement .+= increment

        assemble_u!(
            stiffness,
            residual,
            setup.dh_u,
            setup.dh_d,
            displacement,
            damage,
            material,
            displacement_values,
            damage_values,
        )
        apply_zero!(stiffness, residual, setup.ch_u)
        residual_norm = norm(residual)
        iteration += 1
    end

    return (; iteration, residual_norm, converged = residual_norm <= tolerance)
end

function _solve_damage!(
    stiffness,
    force,
    setup::TensionSetup,
    damage,
    driving_force,
    material::MaterialParams,
    damage_values;
    eta::Real = 0.0,
    dt::Real = 1.0,
    d_old::Union{Nothing,AbstractVector} = nothing,
)
    assemble_d!(
        stiffness, force, setup.dh_d, driving_force, material, damage_values;
        eta, dt, d_old,
    )
    apply!(stiffness, force, setup.ch_d)
    damage .= stiffness \ force
    apply!(damage, setup.ch_d)
    return nothing
end

"""Return `∫ (eta/dt) (d_new-d_old)^2 dΩ`, the BDF1 viscous dissipation."""
function _viscous_dissipation(
    setup::TensionSetup,
    damage_new,
    damage_old,
    eta::Real,
    dt::Real,
    damage_values::CellValues,
)
    eta == 0 && return 0.0
    dissipation = 0.0
    for cell in CellIterator(setup.dh_d)
        reinit!(damage_values, cell)
        new_local = damage_new[celldofs(cell)]
        old_local = damage_old[celldofs(cell)]
        for q_point in 1:getnquadpoints(damage_values)
            increment = function_value(damage_values, q_point, new_local) -
                        function_value(damage_values, q_point, old_local)
            dissipation += eta / dt * increment^2 * getdetJdV(damage_values, q_point)
        end
    end
    return dissipation
end

function _write_staggered_vtk(output_directory, setup::TensionSetup, displacement, damage, step)
    path = joinpath(output_directory, "fracture_step_$step")
    VTKGridFile(path, setup.dh_u) do vtk
        write_solution(vtk, setup.dh_u, displacement)
        write_solution(vtk, setup.dh_d, damage)
    end
    return nothing
end

"""
    solve_staggered(setup, material, config)

Solve the first-order phase-field fracture model with the staggered scheme. At
each load step, the displacement and damage subproblems are solved alternately
until the damage increment satisfies `tol`.

The function returns a `StaggeredResult`. Its `u_final` and `d_final` fields
contain the raw final degree-of-freedom vectors.

`StaggeredConfig` groups physical time, displacement loading, nonlinear solver,
and output controls. A keyword-based compatibility wrapper is defined below.
"""
function solve_staggered(
    setup::TensionSetup,
    material::MaterialParams,
    config::StaggeredConfig,
)
    validate_config(config)
    time_config = config.time
    solver_config = config.solver
    output_config = config.output
    times = _staggered_time_grid(
        time_config.n_steps,
        time_config.total_time,
        time_config.time_points,
        time_config.dt,
    )
    load_at_time = _staggered_load_function(times, config.load)
    vtk_directory = _staggered_output_directory(config)
    output_config.write_vtk && mkpath(vtk_directory)

    right_boundary_dofs = get_right_dofs(setup.grid, setup.dh_u, setup.dir)
    n_displacement_dofs = ndofs(setup.dh_u)
    n_damage_dofs = ndofs(setup.dh_d)

    displacement = zeros(n_displacement_dofs)
    previous_displacement = similar(displacement)
    damage = zeros(n_damage_dofs)
    previous_damage = similar(damage)
    damage_before_iteration = similar(damage)

    apply!(displacement, setup.ch_u)
    apply!(damage, setup.ch_d)
    copyto!(previous_displacement, displacement)
    copyto!(previous_damage, damage)

    quadrature = QuadratureRule{RefQuadrilateral}(2)
    displacement_values = CellValues(quadrature, Lagrange{RefQuadrilateral,1}()^2)
    damage_values = CellValues(quadrature, Lagrange{RefQuadrilateral,1}())
    n_quadrature_points = getncells(setup.grid) * getnquadpoints(displacement_values)
    accepted_history = zeros(n_quadrature_points)
    trial_history = similar(accepted_history)
    current_driving_force = similar(accepted_history)

    displacement_stiffness = allocate_matrix(setup.dh_u)
    damage_stiffness = allocate_matrix(setup.dh_d)
    displacement_residual = zeros(n_displacement_dofs)
    damage_force = zeros(n_damage_dofs)

    displacements = Float64[0.0]
    reaction_forces = Float64[0.0]
    elastic_energies = Float64[0.0]
    surface_energies = Float64[0.0]
    diagnostics = StaggeredDiagnostic[]
    cumulative_viscous_dissipation = 0.0

    started_at = time()
    total_newton_iterations = 0
    output_config.verbose && println("开始 Staggered 交错求解，总步数: $(time_config.n_steps)")
    completed = true

    for step in 1:time_config.n_steps
        physical_time = times[step + 1]
        step_dt = physical_time - times[step]
        load_factor = float(load_at_time(physical_time))
        isfinite(load_factor) || throw(ArgumentError("load_history must return a finite factor"))
        imposed_displacement = load_factor * setup.final_displacement
        update!(setup.ch_u, load_factor)

        copyto!(displacement, previous_displacement)
        copyto!(damage, previous_damage)
        apply!(displacement, setup.ch_u)
        apply!(damage, setup.ch_d)

        output_config.verbose && println(
            "=== 物理时间步 $step / $(time_config.n_steps) | t = ",
            round(physical_time; digits = 6),
            " | 位移: ",
            round(imposed_displacement; digits = 5),
            " ===",
        )

        # These are the accepted n-state.  They are never modified inside the
        # staggered loop, which makes the BDF1 term exactly (d^(n+1)-d^n)/dt.
        d_old = copy(previous_damage)
        history_old = copy(accepted_history)
        staggered_converged = false
        staggered_iterations = 0
        last_damage_increment = Inf
        last_newton_iterations = 0
        last_displacement_residual = NaN
        for staggered_iteration in 1:solver_config.max_staggered_iterations
            staggered_iterations = staggered_iteration
            newton = _solve_displacement!(
                displacement_stiffness,
                displacement_residual,
                setup,
                displacement,
                damage,
                material,
                displacement_values,
                damage_values;
                tolerance = solver_config.tolerance,
                max_iterations = solver_config.max_newton_iterations,
            )
            total_newton_iterations += newton.iteration
            last_newton_iterations = newton.iteration
            last_displacement_residual = newton.residual_norm
            newton.converged || @warn(
                "位移子问题未在最大 Newton 迭代次数内收敛",
                load_step = step,
                staggered_iteration,
                residual_norm = newton.residual_norm,
            )

            # First calculate ψ⁺(u) without mutating accepted history.  The
            # trial maximum is committed only once this physical time step is
            # accepted; a failed step therefore rolls back every state field.
            compute_driving_force!(
                current_driving_force,
                setup.dh_u,
                displacement,
                material,
                displacement_values,
                false,
            )
            if solver_config.enforce_irreversibility
                trial_history .= max.(history_old, current_driving_force)
            else
                copyto!(trial_history, current_driving_force)
            end

            copyto!(damage_before_iteration, damage)
            _solve_damage!(
                damage_stiffness,
                damage_force,
                setup,
                damage,
                trial_history,
                material,
                damage_values,
                eta = solver_config.viscosity,
                dt = step_dt,
                d_old = d_old,
            )
            damage_increment = norm(damage - damage_before_iteration)
            last_damage_increment = damage_increment

            if output_config.verbose && staggered_iteration % 5 == 0
                @info "交错迭代" load_step = step staggered_iteration damage_increment
            end
            if damage_increment < solver_config.tolerance
                staggered_converged = true
                break
            end
        end

        if !staggered_converged
            # Transactional rollback: retain the last accepted time state.
            copyto!(displacement, previous_displacement)
            copyto!(damage, previous_damage)
            copyto!(accepted_history, history_old)
            completed = false
            @warn "物理时间步未收敛，已回滚并停止求解" time_step = step physical_time max_iterations = solver_config.max_staggered_iterations
            push!(diagnostics, StaggeredDiagnostic(
                step = step,
                time = physical_time,
                dt = step_dt,
                load_factor = load_factor,
                imposed_displacement = imposed_displacement,
                staggered_iterations = staggered_iterations,
                newton_iterations = last_newton_iterations,
                displacement_residual = last_displacement_residual,
                converged = false,
                damage_increment = last_damage_increment,
                damage_min = minimum(damage),
                damage_max = maximum(damage),
                viscous_dissipation = 0.0,
                cumulative_viscous_dissipation = cumulative_viscous_dissipation,
            ))
            break
        end

        copyto!(previous_displacement, displacement)
        copyto!(previous_damage, damage)
        copyto!(accepted_history, trial_history)
        viscous_dissipation = _viscous_dissipation(
            setup, damage, d_old, solver_config.viscosity, step_dt, damage_values,
        )
        cumulative_viscous_dissipation += viscous_dissipation
        push!(diagnostics, StaggeredDiagnostic(
            step = step,
            time = physical_time,
            dt = step_dt,
            load_factor = load_factor,
            imposed_displacement = imposed_displacement,
            staggered_iterations = staggered_iterations,
            newton_iterations = last_newton_iterations,
            displacement_residual = last_displacement_residual,
            converged = true,
            damage_increment = last_damage_increment,
            damage_min = minimum(damage),
            damage_max = maximum(damage),
            viscous_dissipation = viscous_dissipation,
            cumulative_viscous_dissipation = cumulative_viscous_dissipation,
        ))

        push!(
            elastic_energies,
            elastic_energy(
                setup.dh_u,
                setup.dh_d,
                displacement,
                damage,
                material,
                displacement_values,
                damage_values,
            ),
        )
        push!(surface_energies, surface_energy(setup.dh_d, damage, material, damage_values))

        reaction_force = compute_reaction_forces(
            right_boundary_dofs,
            displacement_stiffness,
            displacement_residual,
            setup.dh_u,
            setup.dh_d,
            displacement,
            damage,
            material,
            displacement_values,
            damage_values,
        )
        push!(displacements, imposed_displacement)
        push!(reaction_forces, reaction_force)

        if output_config.write_vtk && (step % output_config.vtk_interval == 0 || step == time_config.n_steps)
            _write_staggered_vtk(vtk_directory, setup, displacement, damage, step)
        end
    end

    if output_config.verbose
        output_config.write_vtk && println("仿真结束！VTK 文件保存在 $vtk_directory 目录下。")
        println("总 Newton 迭代次数: $total_newton_iterations")
        println("计算耗时: $(round(time() - started_at; digits = 2)) 秒")
    end

    message = completed ? "all staggered physical-time steps completed" :
              "staggered solver stopped after a non-converged physical-time step"
    return StaggeredResult(
        completed,
        completed,
        message,
        config,
        times[1:length(displacements)],
        displacements,
        reaction_forces,
        elastic_energies,
        surface_energies,
        diagnostics,
        copy(displacement),
        copy(damage),
        cumulative_viscous_dissipation,
        total_newton_iterations,
    )
end

"""
    solve_staggered(setup, material; kwargs...)

Compatibility wrapper for the original keyword-based API. New code should pass
a `StaggeredConfig` explicitly so time, loading, solver, and output settings
remain grouped at the call site. With `return_diagnostics=true`, this wrapper
returns `StaggeredResult`; otherwise it preserves the original four-vector
return value.
"""
function solve_staggered(
    setup::TensionSetup,
    material::MaterialParams;
    n_steps::Integer = 100,
    tol::Real = 1.0e-5,
    max_iter::Integer = 20,
    max_newton_iter::Integer = DEFAULT_MAX_NEWTON_ITERATIONS,
    enforce_irreversibility::Bool = true,
    eta::Real = 0.0,
    dt::Union{Nothing,Real} = nothing,
    total_time::Real = 1.0,
    time_points::Union{Nothing,AbstractVector} = nothing,
    load_history::Union{Nothing,Function} = nothing,
    ramp_time::Union{Nothing,Real} = nothing,
    return_diagnostics::Bool = false,
    output_directory::Union{Nothing,AbstractString} = nothing,
    write_vtk::Bool = true,
    vtk_interval::Integer = DEFAULT_STAGGERED_OUTPUT_INTERVAL,
    verbose::Bool = true,
)
    config = StaggeredConfig(
        time = StaggeredTimeConfig(
            n_steps = Int(n_steps),
            dt = isnothing(dt) ? nothing : Float64(dt),
            total_time = Float64(total_time),
            time_points = isnothing(time_points) ? nothing : Float64.(time_points),
        ),
        load = StaggeredLoadConfig(
            load_history = load_history,
            ramp_time = isnothing(ramp_time) ? nothing : Float64(ramp_time),
        ),
        solver = StaggeredSolverConfig(
            tolerance = Float64(tol),
            max_staggered_iterations = Int(max_iter),
            max_newton_iterations = Int(max_newton_iter),
            enforce_irreversibility = enforce_irreversibility,
            viscosity = Float64(eta),
        ),
        output = StaggeredOutputConfig(
            directory = isnothing(output_directory) ? nothing : String(output_directory),
            write_vtk = write_vtk,
            vtk_interval = Int(vtk_interval),
            verbose = verbose,
        ),
    )
    result = solve_staggered(setup, material, config)
    return return_diagnostics ? result :
           (result.displacements, result.reaction_forces, result.elastic_energies, result.surface_energies)
end
