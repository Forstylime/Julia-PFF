"""
    compute_driving_force!(history, dh_u, displacement, material, cell_values,
                           enforce_irreversibility)

Update the tensile-energy driving force at every quadrature point. With
irreversibility enabled, `history` stores the largest value reached so far;
otherwise it stores the current value.
"""
function compute_driving_force!(
    history::AbstractVector,
    dh_u::DofHandler,
    displacement::AbstractVector,
    material::MaterialParams,
    cell_values::CellValues,
    enforce_irreversibility::Bool,
)
    quadrature_index = 1
    for cell in CellIterator(dh_u)
        reinit!(cell_values, cell)
        local_displacement = displacement[celldofs(cell)]

        for quadrature_point in 1:getnquadpoints(cell_values)
            strain = function_symmetric_gradient(
                cell_values,
                quadrature_point,
                local_displacement,
            )
            tensile_energy, _ = elastic_energy_densities(PlaneStrain(), strain, material)

            if enforce_irreversibility
                history[quadrature_index] = max(history[quadrature_index], tensile_energy)
            else
                history[quadrature_index] = tensile_energy
            end
            quadrature_index += 1
        end
    end
    return history
end

"""
    get_right_dofs(grid, dh_u, component; tol=1e-12)

Return the displacement degrees of freedom on the right boundary for one
component (`1` for x and `2` for y).
"""
function get_right_dofs(grid, dh_u, component::Integer; tol::Real = 1.0e-12)
    component in (1, 2) || throw(ArgumentError("component must be 1 (x) or 2 (y)"))

    node_dofs = zeros(Int, 2, getnnodes(grid))
    for cell_index in 1:getncells(grid)
        cell = getcells(grid, cell_index)
        cell_dofs = celldofs(dh_u, cell_index)
        for (local_node, node_index) in pairs(cell.nodes)
            offset = 2 * (local_node - 1)
            node_dofs[1, node_index] = cell_dofs[offset + 1]
            node_dofs[2, node_index] = cell_dofs[offset + 2]
        end
    end

    x_coordinates = [node.x[1] for node in grid.nodes]
    right_x = maximum(x_coordinates)
    right_nodes = findall(x -> isapprox(x, right_x; atol = tol), x_coordinates)
    right_dofs = [node_dofs[component, node] for node in right_nodes]
    return sort!(unique!(right_dofs))
end

"""
    compute_reaction_forces(reaction_dofs, stiffness, residual, ...)

Reassemble the unconstrained internal-force vector and return the total
reaction on `reaction_dofs`.
"""
function compute_reaction_forces(
    reaction_dofs,
    stiffness,
    residual,
    dh_u,
    dh_d,
    displacement,
    damage,
    material,
    displacement_values,
    damage_values,
)
    fill!(stiffness.nzval, 0.0)
    fill!(residual, 0.0)
    assemble_u!(
        stiffness,
        residual,
        dh_u,
        dh_d,
        displacement,
        damage,
        material,
        displacement_values,
        damage_values,
    )
    return sum(residual[dof] for dof in reaction_dofs)
end
