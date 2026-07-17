using Test
using LinearAlgebra
using Ferrite
using PffSAV
using Tensors

const MAT = RLMMaterialConfig(
    E = 1_000.0,
    nu = 0.25,
    G_c = 1.0,
    ell = 0.2,
    kappa = 1.0e-6,
    mobility = 1.0,
)
const TOL = RLMToleranceConfig()

function small_config(; scalar_residual = 1.0e-10, phase = 1.0, q = 1.0, max_relax = 4)
    return RLMConfig(
        material = MAT,
        mesh = RLMMeshConfig(path = "unused", quadrature_order = 2),
        load = RLMLoadConfig(
            fixed_boundary = "left",
            loaded_boundary = "right",
            component = 1,
            final_displacement = 1.0e-4,
            load_steps = 2,
            initial_damage = 0.0,
        ),
        time = RLMTimeConfig(
            dt = 1.0e-3,
            alpha = 10.0,
            min_relax_steps = 1,
            max_relax_steps = max_relax,
        ),
        tolerances = RLMToleranceConfig(
            scalar_residual = scalar_residual,
            phase = phase,
            q = q,
        ),
        output = RLMOutputConfig(write_csv = false, write_vtk = false, verbose = false),
    )
end

small_grid() = generate_grid(
    Quadrilateral,
    (2, 2),
    Vec(0.0, 0.0),
    Vec(1.0, 1.0),
)

@testset "normalized degradation" begin
    kappa = MAT.kappa
    @test degradation(0.0, kappa) == 1.0
    @test degradation(1.0, kappa) == kappa
    d = 0.37
    h = 1.0e-7
    finite_difference = (degradation(d + h, kappa) - degradation(d - h, kappa)) / (2h)
    @test degradation_derivative(d, kappa) ≈ finite_difference rtol = 1.0e-8
end

@testset "robust plane-strain Miehe split" begin
    lambda = MAT.E * MAT.nu / ((1 + MAT.nu) * (1 - 2MAT.nu))
    mu = MAT.E / (2 * (1 + MAT.nu))
    C0 = PffSAV.plane_strain_elasticity_tensor(MAT)

    epsilon_tension = SymmetricTensor{2, 2, Float64}((0.01, 0.0, 0.0))
    pp, pm, sp, sm, split = miehe_response_2d(epsilon_tension, MAT, TOL)
    @test pp ≈ (0.5lambda + mu) * 0.01^2
    @test pm ≈ 0.0 atol = 1.0e-14
    @test split.epsilon_plus + split.epsilon_minus ≈ epsilon_tension
    @test sp + sm ≈ C0 ⊡ epsilon_tension

    epsilon_compression = SymmetricTensor{2, 2, Float64}((-0.01, 0.0, 0.0))
    pp, pm, sp, sm, split = miehe_response_2d(epsilon_compression, MAT, TOL)
    @test pp ≈ 0.0 atol = 1.0e-14
    @test pm ≈ (0.5lambda + mu) * 0.01^2
    @test split.epsilon_plus + split.epsilon_minus ≈ epsilon_compression
    @test sp + sm ≈ C0 ⊡ epsilon_compression

    epsilon_shear = SymmetricTensor{2, 2, Float64}((0.0, 0.005, 0.0))
    pp, pm, sp, sm, split = miehe_response_2d(epsilon_shear, MAT, TOL)
    @test pp ≈ mu * 0.005^2
    @test pm ≈ mu * 0.005^2
    @test split.principal_min ≈ -0.005
    @test split.principal_max ≈ 0.005

    epsilon_equal_tension = SymmetricTensor{2, 2, Float64}((0.002, 0.0, 0.002))
    pp, pm, _, _, _ = miehe_response_2d(epsilon_equal_tension, MAT, TOL)
    @test pp ≈ 0.5 * (epsilon_equal_tension ⊡ C0 ⊡ epsilon_equal_tension)
    @test pm ≈ 0.0 atol = 1.0e-14
    epsilon_equal_compression = -epsilon_equal_tension
    pp, pm, _, _, _ = miehe_response_2d(epsilon_equal_compression, MAT, TOL)
    @test pp ≈ 0.0 atol = 1.0e-14
    @test pm ≈ 0.5 * (epsilon_equal_compression ⊡ C0 ⊡ epsilon_equal_compression)

    epsilon_reference = SymmetricTensor{2, 2, Float64}((0.003, 0.001, -0.002))
    theta = 0.731
    rotation = [cos(theta) -sin(theta); sin(theta) cos(theta)]
    rotated_matrix = rotation * Matrix(epsilon_reference) * transpose(rotation)
    epsilon_rotated = SymmetricTensor{2, 2, Float64}((
        rotated_matrix[1, 1], rotated_matrix[1, 2], rotated_matrix[2, 2],
    ))
    pp_ref, pm_ref, _, _, _ = miehe_response_2d(epsilon_reference, MAT, TOL)
    pp_rot, pm_rot, _, _, _ = miehe_response_2d(epsilon_rotated, MAT, TOL)
    @test pp_rot ≈ pp_ref rtol = 1.0e-12
    @test pm_rot ≈ pm_ref rtol = 1.0e-12

    cases = (
        zero(epsilon_tension),
        SymmetricTensor{2, 2, Float64}((1.0e-6, 1.0e-18, 1.0e-6)),
        SymmetricTensor{2, 2, Float64}((1.0e-16, 0.0, -1.0e-16)),
        SymmetricTensor{2, 2, Float64}((0.002, 0.003, -0.004)),
    )
    for epsilon in cases
        pp, pm, sp, sm, split = miehe_response_2d(epsilon, MAT, TOL)
        unsplit = 0.5 * (epsilon ⊡ C0 ⊡ epsilon)
        @test split.epsilon_plus + split.epsilon_minus ≈ epsilon atol = 1.0e-14
        @test pp + pm ≈ unsplit atol = 1.0e-14 rtol = 1.0e-12
        @test sp + sm ≈ C0 ⊡ epsilon atol = 1.0e-12 rtol = 1.0e-12
    end
end

@testset "stable scalar equation" begin
    result = solve_rlm_quadratic(1.0, 0.0, -1.0, 1.0, TOL)
    @test result.success
    @test result.q ≈ 1.0
    negative_zero_B = solve_rlm_quadratic(1.0, -0.0, -1.0, 1.0, TOL)
    @test negative_zero_B.success
    @test negative_zero_B.q ≈ 1.0

    cancellation = solve_rlm_quadratic(1.0, 1.0e16, -1.0, 1.0, TOL)
    @test cancellation.success
    @test cancellation.q ≈ 1.0e-16 rtol = 1.0e-12

    two_positive = solve_rlm_quadratic(1.0, -3.0, 2.0, 1.1, TOL)
    @test two_positive.success
    @test two_positive.q == 1.0

    nearly_double = solve_rlm_quadratic(1.0, -2.0, 1.0 + 1.0e-14, 1.0, TOL)
    @test nearly_double.success
    @test nearly_double.discriminant < 0.0
    @test nearly_double.discriminant_used == 0.0

    linear = solve_rlm_quadratic(0.0, 2.0, -4.0, 1.0, TOL)
    @test linear.success
    @test linear.q == 2.0

    no_real = solve_rlm_quadratic(1.0, 0.0, 1.0, 1.0, TOL)
    @test !no_real.success
    @test no_real.code == :negative_discriminant

    no_positive = solve_rlm_quadratic(1.0, 2.0, 1.0, 1.0, TOL)
    @test !no_positive.success
    @test no_positive.code == :no_admissible_positive_root

    zero_residual_tolerance = RLMToleranceConfig(scalar_residual = 0.0)
    rejected_residual = solve_rlm_quadratic(1.0, 0.0, -2.0, 1.0, zero_residual_tolerance)
    @test !rejected_residual.success
    @test rejected_residual.code == :no_admissible_positive_root
    @test isfinite(rejected_residual.residual)
    @test rejected_residual.residual > 0.0
end

@testset "affine branches and one-step identities" begin
    problem = build_rlm_problem(small_config(); grid = small_grid())
    @test issymmetric(problem.K_u)
    @test issymmetric(problem.M_d)
    @test issymmetric(problem.K_AT2)
    @test isposdef(Symmetric(Matrix(problem.K_u_constrained)))
    @test isposdef(Symmetric(Matrix(problem.K_d)))

    update!(problem.ch_u, 1.0)
    state = RLMState(zeros(ndofs(problem.dh_u)), zeros(ndofs(problem.dh_d)), 1.0, 0.0)
    apply!(state.u, problem.ch_u)
    state.P = rlm_nonlinear_energy(problem, state.u, state.d)
    old_state = copy(state)
    trial = compute_rlm_bdf1_trial(problem, state)

    @test state.u == old_state.u
    @test state.d == old_state.d
    @test state.q == old_state.q
    @test state.P == old_state.P

    u_a_checked = copy(trial.u_a)
    apply!(u_a_checked, problem.ch_u)
    @test trial.u_a ≈ u_a_checked
    u_b_checked = copy(trial.u_b)
    apply_zero!(u_b_checked, problem.ch_u)
    @test trial.u_b ≈ u_b_checked

    c1_quadratic = -dot(trial.u_b, problem.K_u * trial.u_b) -
                   dot(trial.d_b, problem.K_d * trial.d_b)
    @test trial.c1 <= 1.0e-12
    @test trial.c1 ≈ c1_quadratic rtol = 1.0e-9 atol = 1.0e-12

    free = free_dofs(problem.ch_u)
    equilibrium = problem.K_u * trial.u + trial.q * trial.n_u - problem.f_ext
    @test norm(equilibrium[free]) <= 1.0e-9
    inverse_M_dt = 1 / (problem.config.material.mobility * problem.config.time.dt)
    phase_residual = inverse_M_dt * problem.M_d * (trial.d - state.d) +
                     problem.K_AT2 * trial.d + trial.q * trial.n_d
    @test norm(phase_residual) <= 1.0e-9
    @test scalar_equation_residual(
        trial.A, trial.B, trial.C, trial.q,
        problem.config.tolerances.scalar_denominator_epsilon,
    ) <= problem.config.tolerances.scalar_residual

    proxy_old = rlm_proxy_energy(problem, state.u, state.d, state.q, state.P)
    @test trial.proxy_energy <= proxy_old + 1.0e-11
    @test trial.energy_balance_residual <= problem.config.tolerances.energy_balance_rel
    @test trial.healing >= 0.0
end

@testset "fixed external-load functional" begin
    base = small_config()
    load = RLMLoadConfig(
        fixed_boundary = "left",
        loaded_boundary = "right",
        component = 1,
        final_displacement = 1.0e-4,
        load_steps = 1,
        body_force = (1.0, 2.0),
        traction_boundary = "top",
        traction = (0.0, 3.0),
    )
    config = RLMConfig(
        material = base.material,
        mesh = base.mesh,
        load = load,
        time = base.time,
        tolerances = base.tolerances,
        output = base.output,
    )
    problem = build_rlm_problem(config; grid = generate_grid(
        Quadrilateral, (1, 1), Vec(0.0, 0.0), Vec(1.0, 1.0),
    ))
    @test norm(problem.f_ext) > 0.0
    # Partition of unity: unit-area body load contributes 1+2 and the
    # unit-length top traction contributes 3 to the sum of vector entries.
    @test sum(problem.f_ext) ≈ 6.0
end

@testset "transactional rollback" begin
    config = small_config(scalar_residual = 0.0)
    problem = build_rlm_problem(config; grid = small_grid())
    update!(problem.ch_u, 1.0)
    state = RLMState(zeros(ndofs(problem.dh_u)), zeros(ndofs(problem.dh_d)), 1.0, 0.0)
    apply!(state.u, problem.ch_u)
    state.P = rlm_nonlinear_energy(problem, state.u, state.d)
    snapshot = copy(state)
    @test_throws PffSAV.RLMStepFailure compute_rlm_bdf1_trial(problem, state)
    @test state.u == snapshot.u
    @test state.d == snapshot.d
    @test state.q == snapshot.q
    @test state.P == snapshot.P
end

@testset "load-relax and diagnostics" begin
    problem = build_rlm_problem(small_config(); grid = small_grid())
    result = solve_rlm_bdf1(problem)
    @test result.success
    @test result.converged
    @test count(d -> d.status == "load_reset", result.diagnostics) == 2
    accepted = filter(d -> d.status == "accepted", result.diagnostics)
    @test length(accepted) == 2
    for diagnostic in accepted
        @test isfinite(diagnostic.raw_energy)
        @test isfinite(diagnostic.proxy_energy)
        @test isfinite(diagnostic.q_minus_one)
        @test diagnostic.c1 <= 1.0e-12
        @test isfinite(diagnostic.discriminant)
        @test diagnostic.scalar_residual <= problem.config.tolerances.scalar_residual
        @test diagnostic.phase_increment >= 0.0
        @test diagnostic.healing >= 0.0
    end

    mktempdir() do directory
        path = joinpath(directory, "diagnostics.csv")
        write_rlm_diagnostics(path, result.diagnostics)
        text = read(path, String)
        @test occursin("raw_energy", text)
        @test occursin("proxy_energy", text)
        @test occursin("q_minus_one", text)
        @test occursin("healing", text)
    end
end
