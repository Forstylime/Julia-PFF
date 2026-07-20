using Ferrite
using LinearAlgebra

const DEFAULT_STAGGERED_OUTPUT_INTERVAL = 5
const DEFAULT_MAX_NEWTON_ITERATIONS = 10

function _validate_staggered_options(
    n_steps::Integer,
    tolerance::Real,
    max_staggered_iterations::Integer,
    max_newton_iterations::Integer,
    vtk_interval::Integer,
)
    n_steps > 0 || throw(ArgumentError("n_steps must be positive"))
    tolerance > 0 || throw(ArgumentError("tol must be positive"))
    max_staggered_iterations > 0 || throw(ArgumentError("max_iter must be positive"))
    max_newton_iterations > 0 || throw(ArgumentError("max_newton_iter must be positive"))
    vtk_interval > 0 || throw(ArgumentError("vtk_interval must be positive"))
    return nothing
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
    solve_staggered(setup, material; kwargs...)

Solve the first-order phase-field fracture model with the staggered scheme. At
each load step, the displacement and damage subproblems are solved alternately
until the damage increment satisfies `tol`.

The function returns `(displacements, reaction_forces, elastic_energies,
surface_energies)`, preserving the original public API.

# Keywords
- `n_steps=100`: number of displacement-controlled load steps.
- `tol=1e-5`: convergence tolerance for both inner solves.
- `max_iter=20`: maximum staggered iterations per load step.
- `max_newton_iter=10`: maximum Newton iterations per displacement solve.
- `enforce_irreversibility=true`: use the maximum tensile-energy history field.
- `eta=0`: phase-field viscosity. Positive values enable BDF1 real-time evolution.
- `dt=nothing`: physical BDF1 time step; mutually exclusive with `time_points`.
- `total_time=1`: end time used when `time_points` and `dt` are omitted.
- `time_points=nothing`: strictly increasing physical-time nodes (including zero).
- `load_history=nothing`: function of physical time returning the load factor.
- `return_diagnostics=false`: return a named result including time-step diagnostics.
- `output_directory=nothing`: VTK directory. A mode-specific default is used
  when omitted.
- `write_vtk=true`: enable VTK output.
- `vtk_interval=5`: write every N load steps, plus the final step.
- `verbose=true`: print load-step progress and the run summary.
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
    return_diagnostics::Bool = false,
    output_directory::Union{Nothing,AbstractString} = nothing,
    write_vtk::Bool = true,
    vtk_interval::Integer = DEFAULT_STAGGERED_OUTPUT_INTERVAL,
    verbose::Bool = true,
)
    _validate_staggered_options(n_steps, tol, max_iter, max_newton_iter, vtk_interval)
    eta >= 0 || throw(ArgumentError("eta must be non-negative"))
    times = _staggered_time_grid(n_steps, total_time, time_points, dt)
    load_at_time = isnothing(load_history) ? (t -> t / times[end]) : load_history

    vtk_directory = if isnothing(output_directory)
        suffix = eta > 0 ? "staggered_bdf1" :
                 (enforce_irreversibility ? "staggered" : "staggered2")
        joinpath("data", "sims", suffix)
    else
        String(output_directory)
    end
    write_vtk && mkpath(vtk_directory)

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
    diagnostics = NamedTuple[]
    cumulative_viscous_dissipation = 0.0

    started_at = time()
    total_newton_iterations = 0
    verbose && println("开始 Staggered 交错求解，总步数: $n_steps")

    for step in 1:n_steps
        physical_time = times[step + 1]
        dt = physical_time - times[step]
        load_factor = float(load_at_time(physical_time))
        isfinite(load_factor) || throw(ArgumentError("load_history must return a finite factor"))
        imposed_displacement = load_factor * setup.final_displacement
        update!(setup.ch_u, load_factor)

        copyto!(displacement, previous_displacement)
        copyto!(damage, previous_damage)
        apply!(displacement, setup.ch_u)
        apply!(damage, setup.ch_d)

        verbose && println(
            "=== 物理时间步 $step / $n_steps | t = ",
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
        for staggered_iteration in 1:max_iter
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
                tolerance = tol,
                max_iterations = max_newton_iter,
            )
            total_newton_iterations += newton.iteration
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
            if enforce_irreversibility
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
                eta = eta,
                dt = dt,
                d_old = d_old,
            )
            damage_increment = norm(damage - damage_before_iteration)
            last_damage_increment = damage_increment

            if verbose && staggered_iteration % 5 == 0
                @info "交错迭代" load_step = step staggered_iteration damage_increment
            end
            if damage_increment < tol
                staggered_converged = true
                break
            end
        end

        if !staggered_converged
            # Transactional rollback: retain the last accepted time state.
            copyto!(displacement, previous_displacement)
            copyto!(damage, previous_damage)
            copyto!(accepted_history, history_old)
            @warn "物理时间步未收敛，已回滚并停止求解" time_step = step physical_time max_iterations = max_iter
            push!(diagnostics, (
                step, time = physical_time, dt, load_factor, staggered_iterations,
                converged = false, damage_increment = last_damage_increment,
                damage_min = minimum(damage), damage_max = maximum(damage),
                viscous_dissipation = 0.0,
                cumulative_viscous_dissipation,
            ))
            break
        end

        copyto!(previous_displacement, displacement)
        copyto!(previous_damage, damage)
        copyto!(accepted_history, trial_history)
        viscous_dissipation = _viscous_dissipation(
            setup, damage, d_old, eta, dt, damage_values,
        )
        cumulative_viscous_dissipation += viscous_dissipation
        push!(diagnostics, (
            step, time = physical_time, dt, load_factor, staggered_iterations,
            converged = true, damage_increment = last_damage_increment,
            damage_min = minimum(damage), damage_max = maximum(damage),
            viscous_dissipation, cumulative_viscous_dissipation,
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

        if write_vtk && (step % vtk_interval == 0 || step == n_steps)
            _write_staggered_vtk(vtk_directory, setup, displacement, damage, step)
        end
    end

    if verbose
        write_vtk && println("仿真结束！VTK 文件保存在 $vtk_directory 目录下。")
        println("总 Newton 迭代次数: $total_newton_iterations")
        println("计算耗时: $(round(time() - started_at; digits = 2)) 秒")
    end

    if return_diagnostics
        return (
            displacements, reaction_forces, elastic_energies, surface_energies,
            times = times[1:length(displacements)], diagnostics,
            cumulative_viscous_dissipation,
        )
    end
    return displacements, reaction_forces, elastic_energies, surface_energies
end
