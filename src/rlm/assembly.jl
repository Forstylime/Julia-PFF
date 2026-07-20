struct RLMProblem{G,DHU,DHD,CHU,CVU,CVD,FVU,C,RHS,FU,FD}
    config::RLMConfig
    grid::G
    dh_u::DHU
    dh_d::DHD
    ch_u::CHU
    cellvalues_u::CVU
    cellvalues_d::CVD
    facetvalues_u::FVU
    C0::C
    K_u::SparseMatrixCSC{Float64, Int}
    M_d::SparseMatrixCSC{Float64, Int}
    K_AT2::SparseMatrixCSC{Float64, Int}
    K_d::SparseMatrixCSC{Float64, Int}
    K_u_constrained::SparseMatrixCSC{Float64, Int}
    f_ext::Vector{Float64}
    f_body_reference::Vector{Float64}
    f_traction_reference::Vector{Float64}
    loaded_component_dofs::Vector{Int}
    rhs_data_u::RHS
    factor_u::FU
    factor_d::FD
end

function _boundary_component_dofs(dh_u, grid, boundary::String, component::Int)
    facets = _require_facetset(grid, boundary)
    scalar_interpolation = Lagrange{RefQuadrilateral, 1}()
    local_facet_nodes = Ferrite.dirichlet_facetdof_indices(scalar_interpolation)
    dofs = Int[]
    for (cell_index, facet_index) in facets
        cell_dofs = celldofs(dh_u, cell_index)
        for local_node in local_facet_nodes[facet_index]
            push!(dofs, cell_dofs[2 * (local_node - 1) + component])
        end
    end
    sort!(unique!(dofs))
    isempty(dofs) && throw(ArgumentError(
        "loaded boundary '$boundary' has no component-$component displacement dofs",
    ))
    return dofs
end

function _require_facetset(grid, name::String)
    try
        return getfacetset(grid, name)
    catch err
        throw(ArgumentError("mesh does not contain required facet set '$name': $(sprint(showerror, err))"))
    end
end

function _validate_q1_quadrilateral_grid(grid)
    cells = getcells(grid)
    isempty(cells) && throw(ArgumentError("mesh contains no cells"))
    all(cell -> cell isa Quadrilateral, cells) || throw(ArgumentError(
        "the first-stage RLM solver supports only four-node Quadrilateral cells; " *
        "the imported mesh contains $(unique(typeof.(cells)))",
    ))
    return grid
end

function _create_rlm_dofhandlers(grid)
    interpolation_u = Lagrange{RefQuadrilateral, 1}()^2
    interpolation_d = Lagrange{RefQuadrilateral, 1}()

    dh_u = DofHandler(grid)
    add!(dh_u, :u, interpolation_u)
    close!(dh_u)

    dh_d = DofHandler(grid)
    add!(dh_d, :d, interpolation_d)
    close!(dh_d)
    return dh_u, dh_d, interpolation_u, interpolation_d
end

function _create_rlm_constraints(dh_u, grid, load::RLMLoadConfig)
    fixed = _require_facetset(grid, load.fixed_boundary)
    loaded = _require_facetset(grid, load.loaded_boundary)
    ch_u = ConstraintHandler(dh_u)
    displacement_amplitude = load.displacement_amplitude
    fixed_condition = Dirichlet(:u, fixed, (x, t) -> Vec(0.0, 0.0), [1, 2])
    loaded_condition = Dirichlet(
        :u,
        loaded,
        (x, time) -> displacement_amplitude * history_value(load.displacement_history, time),
        load.component,
    )
    # Ferrite resolves a shared corner degree of freedom with the last-added
    # condition. Make that otherwise ambiguous choice explicit in the config.
    if load.overlap_policy == :loaded
        add!(ch_u, fixed_condition)
        add!(ch_u, loaded_condition)
    else
        add!(ch_u, loaded_condition)
        add!(ch_u, fixed_condition)
    end
    close!(ch_u)
    update!(ch_u, 0.0)
    return ch_u
end

function _assemble_constant_matrices!(
    K_u,
    M_d,
    K_AT2,
    dh_u,
    dh_d,
    cellvalues_u,
    cellvalues_d,
    C0,
    material,
)
    assembler_u = start_assemble(K_u)
    assembler_mass = start_assemble(M_d)
    assembler_at2 = start_assemble(K_AT2)
    n_u = getnbasefunctions(cellvalues_u)
    n_d = getnbasefunctions(cellvalues_d)
    Kue = zeros(n_u, n_u)
    Mde = zeros(n_d, n_d)
    Kde = zeros(n_d, n_d)

    for (cell_u, cell_d) in zip(CellIterator(dh_u), CellIterator(dh_d))
        reinit!(cellvalues_u, cell_u)
        reinit!(cellvalues_d, cell_d)
        fill!(Kue, 0.0)
        fill!(Mde, 0.0)
        fill!(Kde, 0.0)

        for qp in 1:getnquadpoints(cellvalues_u)
            dOmega = getdetJdV(cellvalues_u, qp)
            for i in 1:n_u
                strain_i = shape_symmetric_gradient(cellvalues_u, qp, i)
                for j in 1:n_u
                    strain_j = shape_symmetric_gradient(cellvalues_u, qp, j)
                    Kue[i, j] += (strain_i ⊡ C0 ⊡ strain_j) * dOmega
                end
            end
        end

        for qp in 1:getnquadpoints(cellvalues_d)
            dOmega = getdetJdV(cellvalues_d, qp)
            for i in 1:n_d
                value_i = shape_value(cellvalues_d, qp, i)
                gradient_i = shape_gradient(cellvalues_d, qp, i)
                for j in 1:n_d
                    value_j = shape_value(cellvalues_d, qp, j)
                    gradient_j = shape_gradient(cellvalues_d, qp, j)
                    Mde[i, j] += value_i * value_j * dOmega
                    Kde[i, j] += (
                        material.G_c / material.ell * value_i * value_j +
                        material.G_c * material.ell * (gradient_i ⋅ gradient_j)
                    ) * dOmega
                end
            end
        end

        assemble!(assembler_u, celldofs(cell_u), Kue)
        assemble!(assembler_mass, celldofs(cell_d), Mde)
        assemble!(assembler_at2, celldofs(cell_d), Kde)
    end
    return nothing
end

function _assemble_body_force!(
    f_body,
    dh_u,
    grid,
    cellvalues_u,
    body_force::NTuple{2, Float64},
)
    fill!(f_body, 0.0)
    body = Vec(body_force)
    if body != zero(body)
        n_u = getnbasefunctions(cellvalues_u)
        fe = zeros(n_u)
        for cell in CellIterator(dh_u)
            reinit!(cellvalues_u, cell)
            fill!(fe, 0.0)
            for qp in 1:getnquadpoints(cellvalues_u)
                dOmega = getdetJdV(cellvalues_u, qp)
                for i in 1:n_u
                    test_value = shape_value(cellvalues_u, qp, i)
                    fe[i] += (body ⋅ test_value) * dOmega
                end
            end
            assemble!(f_body, celldofs(cell), fe)
        end
    end

    return f_body
end

function _assemble_traction_force!(
    f_traction,
    dh_u,
    grid,
    facetvalues_u,
    load::RLMLoadConfig,
)
    fill!(f_traction, 0.0)
    if load.traction_boundary !== nothing
        facets = _require_facetset(grid, load.traction_boundary)
        traction = Vec(load.traction)
        n_u = getnbasefunctions(facetvalues_u)
        fe = zeros(n_u)
        for facet in FacetIterator(dh_u, facets)
            reinit!(facetvalues_u, facet)
            fill!(fe, 0.0)
            for qp in 1:getnquadpoints(facetvalues_u)
                dGamma = getdetJdV(facetvalues_u, qp)
                for i in 1:n_u
                    test_value = shape_value(facetvalues_u, qp, i)
                    fe[i] += (traction ⋅ test_value) * dGamma
                end
            end
            assemble!(f_traction, celldofs(facet), fe)
        end
    elseif load.traction != (0.0, 0.0)
        throw(ArgumentError("a nonzero traction requires load.traction_boundary"))
    end
    return f_traction
end

"""Update physical-time external force coefficients without reassembling FEM integrals."""
function update_rlm_external_force!(problem::RLMProblem, time::Real)
    load = problem.config.load
    body_factor = history_value(load.body_force_history, time)
    traction_factor = history_value(load.traction_history, time)
    problem.f_ext .= body_factor .* problem.f_body_reference .+
                     traction_factor .* problem.f_traction_reference
    return problem.f_ext
end

"""
    build_rlm_problem(config; grid=nothing)

Build and factor the two constant BDF1 matrices. Tests may inject a Ferrite grid;
normal runs import `config.mesh.path` with FerriteGmsh.
"""
function build_rlm_problem(config::RLMConfig; grid = nothing)
    validate_config(config)
    if grid === nothing
        isfile(config.mesh.path) || throw(ArgumentError(
            "mesh file '$(config.mesh.path)' was not found; place l_shape.msh in the project root " *
            "or set RLMConfig.mesh.path explicitly",
        ))
        grid = FerriteGmsh.togrid(config.mesh.path)
    end
    _validate_q1_quadrilateral_grid(grid)

    dh_u, dh_d, interpolation_u, interpolation_d = _create_rlm_dofhandlers(grid)
    ch_u = _create_rlm_constraints(dh_u, grid, config.load)
    quadrature = QuadratureRule{RefQuadrilateral}(config.mesh.quadrature_order)
    facet_quadrature = FacetQuadratureRule{RefQuadrilateral}(config.mesh.quadrature_order)
    cellvalues_u = CellValues(quadrature, interpolation_u)
    cellvalues_d = CellValues(quadrature, interpolation_d)
    facetvalues_u = FacetValues(facet_quadrature, interpolation_u)

    K_u = allocate_matrix(dh_u)
    M_d = allocate_matrix(dh_d)
    K_AT2 = allocate_matrix(dh_d)
    C0 = plane_strain_elasticity_tensor(config.material)
    _assemble_constant_matrices!(
        K_u,
        M_d,
        K_AT2,
        dh_u,
        dh_d,
        cellvalues_u,
        cellvalues_d,
        C0,
        config.material,
    )

    viscosity_over_dt = config.material.viscosity / config.time.dt
    K_d = K_AT2 + viscosity_over_dt * M_d
    f_ext = zeros(ndofs(dh_u))
    f_body_reference = zeros(ndofs(dh_u))
    f_traction_reference = zeros(ndofs(dh_u))
    _assemble_body_force!(
        f_body_reference,
        dh_u,
        grid,
        cellvalues_u,
        config.load.body_force,
    )
    _assemble_traction_force!(
        f_traction_reference,
        dh_u,
        grid,
        facetvalues_u,
        config.load,
    )

    rhs_data_u = get_rhs_data(ch_u, K_u)
    K_u_constrained = copy(K_u)
    apply!(K_u_constrained, ch_u)
    factor_u = cholesky(Symmetric(K_u_constrained))
    factor_d = cholesky(Symmetric(K_d))
    loaded_component_dofs = _boundary_component_dofs(
        dh_u,
        grid,
        config.load.loaded_boundary,
        config.load.component,
    )

    return RLMProblem(
        config,
        grid,
        dh_u,
        dh_d,
        ch_u,
        cellvalues_u,
        cellvalues_d,
        facetvalues_u,
        C0,
        K_u,
        M_d,
        K_AT2,
        K_d,
        K_u_constrained,
        f_ext,
        f_body_reference,
        f_traction_reference,
        loaded_component_dofs,
        rhs_data_u,
        factor_u,
        factor_d,
    )
end

function assemble_rlm_nonlinear_forces!(
    n_u::Vector{Float64},
    n_d::Vector{Float64},
    problem::RLMProblem,
    u::Vector{Float64},
    d::Vector{Float64},
)
    fill!(n_u, 0.0)
    fill!(n_d, 0.0)
    cv_u = problem.cellvalues_u
    cv_d = problem.cellvalues_d
    material = problem.config.material
    tolerances = problem.config.tolerances
    n_base_u = getnbasefunctions(cv_u)
    n_base_d = getnbasefunctions(cv_d)
    f_u = zeros(n_base_u)
    f_d = zeros(n_base_d)
    u_local = zeros(n_base_u)
    d_local = zeros(n_base_d)

    for (cell_u, cell_d) in zip(CellIterator(problem.dh_u), CellIterator(problem.dh_d))
        reinit!(cv_u, cell_u)
        reinit!(cv_d, cell_d)
        fill!(f_u, 0.0)
        fill!(f_d, 0.0)
        dofs_u = celldofs(cell_u)
        dofs_d = celldofs(cell_d)
        @inbounds for i in eachindex(u_local)
            u_local[i] = u[dofs_u[i]]
        end
        @inbounds for i in eachindex(d_local)
            d_local[i] = d[dofs_d[i]]
        end

        for qp in 1:getnquadpoints(cv_u)
            dOmega = getdetJdV(cv_u, qp)
            strain = function_symmetric_gradient(cv_u, qp, u_local)
            damage = function_value(cv_d, qp, d_local)
            psi_plus, _, sigma_plus, _, _ =
                miehe_response_2d(strain, material, tolerances)
            g_minus_one = degradation(damage, material.kappa) - 1.0
            g_prime = degradation_derivative(damage, material.kappa)

            for i in 1:n_base_u
                test_strain = shape_symmetric_gradient(cv_u, qp, i)
                f_u[i] += g_minus_one * (sigma_plus ⊡ test_strain) * dOmega
            end
            for i in 1:n_base_d
                test_value = shape_value(cv_d, qp, i)
                f_d[i] += g_prime * psi_plus * test_value * dOmega
            end
        end
        assemble!(n_u, celldofs(cell_u), f_u)
        assemble!(n_d, celldofs(cell_d), f_d)
    end
    return n_u, n_d
end

function rlm_nonlinear_energy(
    problem::RLMProblem,
    u::Vector{Float64},
    d::Vector{Float64},
)
    energy = 0.0
    cv_u = problem.cellvalues_u
    cv_d = problem.cellvalues_d
    material = problem.config.material
    tolerances = problem.config.tolerances
    u_local = zeros(getnbasefunctions(cv_u))
    d_local = zeros(getnbasefunctions(cv_d))
    for (cell_u, cell_d) in zip(CellIterator(problem.dh_u), CellIterator(problem.dh_d))
        reinit!(cv_u, cell_u)
        reinit!(cv_d, cell_d)
        dofs_u = celldofs(cell_u)
        dofs_d = celldofs(cell_d)
        @inbounds for i in eachindex(u_local)
            u_local[i] = u[dofs_u[i]]
        end
        @inbounds for i in eachindex(d_local)
            d_local[i] = d[dofs_d[i]]
        end
        for qp in 1:getnquadpoints(cv_u)
            strain = function_symmetric_gradient(cv_u, qp, u_local)
            damage = function_value(cv_d, qp, d_local)
            psi_plus, _, _, _, _ = miehe_response_2d(strain, material, tolerances)
            energy += (degradation(damage, material.kappa) - 1.0) *
                      psi_plus * getdetJdV(cv_u, qp)
        end
    end
    return energy
end

"""Return degraded positive and undegraded negative elastic energies by quadrature."""
function rlm_elastic_split_energies(problem::RLMProblem, u::Vector{Float64}, d::Vector{Float64})
    positive = 0.0
    negative = 0.0
    cv_u, cv_d = problem.cellvalues_u, problem.cellvalues_d
    material, tolerances = problem.config.material, problem.config.tolerances
    u_local = zeros(getnbasefunctions(cv_u))
    d_local = zeros(getnbasefunctions(cv_d))
    for (cell_u, cell_d) in zip(CellIterator(problem.dh_u), CellIterator(problem.dh_d))
        reinit!(cv_u, cell_u); reinit!(cv_d, cell_d)
        dofs_u, dofs_d = celldofs(cell_u), celldofs(cell_d)
        @inbounds for i in eachindex(u_local)
            u_local[i] = u[dofs_u[i]]
        end
        @inbounds for i in eachindex(d_local)
            d_local[i] = d[dofs_d[i]]
        end
        for qp in 1:getnquadpoints(cv_u)
            strain = function_symmetric_gradient(cv_u, qp, u_local)
            damage = function_value(cv_d, qp, d_local)
            psi_plus, psi_minus, _, _, _ = miehe_response_2d(strain, material, tolerances)
            dOmega = getdetJdV(cv_u, qp)
            positive += degradation(damage, material.kappa) * psi_plus * dOmega
            negative += psi_minus * dOmega
        end
    end
    return positive, negative
end

@inline function rlm_quadratic_energy(problem::RLMProblem, u, d)
    return 0.5 * dot(u, problem.K_u * u) +
           0.5 * dot(d, problem.K_AT2 * d) -
           dot(problem.f_ext, u)
end

@inline function rlm_raw_energy(problem::RLMProblem, u, d)
    return rlm_quadratic_energy(problem, u, d) + rlm_nonlinear_energy(problem, u, d)
end

@inline function rlm_proxy_energy(problem::RLMProblem, u, d, q, P)
    return rlm_quadratic_energy(problem, u, d) + P +
           problem.config.time.alpha * (q^2 - 1.0)
end

function phase_field_metrics(problem::RLMProblem, d_new, d_old)
    cv_d = problem.cellvalues_d
    increment_squared = 0.0
    healing_squared = 0.0
    norm_new_squared = 0.0
    d_new_local = zeros(getnbasefunctions(cv_d))
    d_old_local = similar(d_new_local)
    # Q1 values are convex combinations of nodal values on the reference cell,
    # so nodal extrema report the discrete field bounds more faithfully than
    # sampling only the interior Gauss points.
    min_damage = minimum(d_new)
    max_damage = maximum(d_new)
    for cell in CellIterator(problem.dh_d)
        reinit!(cv_d, cell)
        dofs = celldofs(cell)
        @inbounds for i in eachindex(d_new_local)
            d_new_local[i] = d_new[dofs[i]]
            d_old_local[i] = d_old[dofs[i]]
        end
        for qp in 1:getnquadpoints(cv_d)
            dOmega = getdetJdV(cv_d, qp)
            damage_new = function_value(cv_d, qp, d_new_local)
            damage_old = function_value(cv_d, qp, d_old_local)
            increment = damage_new - damage_old
            increment_squared += increment^2 * dOmega
            healing_squared += min(increment, 0.0)^2 * dOmega
            norm_new_squared += damage_new^2 * dOmega
        end
    end
    increment_norm = sqrt(max(increment_squared, 0.0))
    healing_norm = sqrt(max(healing_squared, 0.0))
    new_norm = sqrt(max(norm_new_squared, 0.0))
    relative_increment = increment_norm / max(1.0, new_norm)
    return increment_norm, relative_increment, healing_norm, min_damage, max_damage
end
