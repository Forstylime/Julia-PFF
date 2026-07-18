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
    damage_values,
)
    assemble_d!(stiffness, force, setup.dh_d, driving_force, material, damage_values)
    apply!(stiffness, force, setup.ch_d)
    damage .= stiffness \ force
    apply!(damage, setup.ch_d)
    return nothing
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
    output_directory::Union{Nothing,AbstractString} = nothing,
    write_vtk::Bool = true,
    vtk_interval::Integer = DEFAULT_STAGGERED_OUTPUT_INTERVAL,
    verbose::Bool = true,
)
    _validate_staggered_options(n_steps, tol, max_iter, max_newton_iter, vtk_interval)

    vtk_directory = if isnothing(output_directory)
        suffix = enforce_irreversibility ? "staggered" : "staggered2"
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
    driving_force = zeros(n_quadrature_points)

    displacement_stiffness = allocate_matrix(setup.dh_u)
    damage_stiffness = allocate_matrix(setup.dh_d)
    displacement_residual = zeros(n_displacement_dofs)
    damage_force = zeros(n_damage_dofs)

    displacements = Float64[0.0]
    reaction_forces = Float64[0.0]
    elastic_energies = Float64[0.0]
    surface_energies = Float64[0.0]

    started_at = time()
    total_newton_iterations = 0
    verbose && println("开始 Staggered 交错求解，总步数: $n_steps")

    for step in 1:n_steps
        load_factor = step / n_steps
        imposed_displacement = load_factor * setup.final_displacement
        update!(setup.ch_u, load_factor)

        copyto!(displacement, previous_displacement)
        copyto!(damage, previous_damage)
        apply!(displacement, setup.ch_u)
        apply!(damage, setup.ch_d)

        verbose && println(
            "=== 载荷步 $step / $n_steps | 位移: ",
            round(imposed_displacement; digits = 5),
            " ===",
        )

        staggered_converged = false
        for staggered_iteration in 1:max_iter
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

            compute_driving_force!(
                driving_force,
                setup.dh_u,
                displacement,
                material,
                displacement_values,
                enforce_irreversibility,
            )

            copyto!(damage_before_iteration, damage)
            _solve_damage!(
                damage_stiffness,
                damage_force,
                setup,
                damage,
                driving_force,
                material,
                damage_values,
            )
            damage_increment = norm(damage - damage_before_iteration)

            if verbose && staggered_iteration % 5 == 0
                @info "交错迭代" load_step = step staggered_iteration damage_increment
            end
            if damage_increment < tol
                staggered_converged = true
                break
            end
        end

        staggered_converged || @warn(
            "载荷步未在最大交错迭代次数内收敛",
            load_step = step,
            max_iterations = max_iter,
        )

        copyto!(previous_displacement, displacement)
        copyto!(previous_damage, damage)

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

    return displacements, reaction_forces, elastic_energies, surface_energies
end
