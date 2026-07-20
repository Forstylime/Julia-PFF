using Test
using LinearAlgebra
using Ferrite
using PffSAV
using Tensors

const MAT = RLMMaterialConfig(E = 1_000.0, nu = 0.25, G_c = 1.0, ell = 0.2,
    kappa = 1.0e-6, viscosity = 1.0)
const TOL = RLMToleranceConfig()
small_grid() = generate_grid(Quadrilateral, (2, 2), Vec(0.0, 0.0), Vec(1.0, 1.0))

function history(values; times = [0.0, 0.001, 0.002])
    RLMPiecewiseLinearHistory(times, values)
end

function small_config(; displacement = history([0.0, 1.0, 1.0]), body = history([0.0, 0.0, 0.0]),
    traction = history([0.0, 0.0, 0.0]), initial_damage = 0.0, dt = 0.001, final_time = 0.002,
    scalar_residual = 1.0e-10)
    RLMConfig(material = MAT, mesh = RLMMeshConfig(path = "unused", quadrature_order = 2),
        load = RLMLoadConfig(fixed_boundary = "left", loaded_boundary = "right", component = 1,
            displacement_amplitude = 1.0e-4, displacement_history = displacement,
            initial_damage = initial_damage, body_force = (1.0, 2.0), body_force_history = body,
            traction_boundary = "top", traction = (0.0, 3.0), traction_history = traction),
        time = RLMTimeConfig(final_time = final_time, dt = dt, alpha = 10.0),
        tolerances = RLMToleranceConfig(scalar_residual = scalar_residual),
        output = RLMOutputConfig(write_csv = false, write_vtk = false, verbose = false))
end

@testset "piecewise-linear physical-time histories" begin
    h = RLMPiecewiseLinearHistory([0, 1, 3], [0, 2, 1])
    @test history_value(h, 0.5) == 1.0
    @test history_value(h, 2.0) == 1.5
    @test_throws ArgumentError history_value(h, -0.1)
    @test_throws ArgumentError RLMPiecewiseLinearHistory([0, 1], [0])
    @test_throws ArgumentError RLMPiecewiseLinearHistory([0, 0], [0, 1])
end

@testset "Miehe split and scalar root" begin
    lambda = MAT.E * MAT.nu / ((1 + MAT.nu) * (1 - 2MAT.nu))
    mu = MAT.E / (2 * (1 + MAT.nu))
    C0 = PffSAV.plane_strain_elasticity_tensor(MAT)
    tension = SymmetricTensor{2, 2, Float64}((0.01, 0.0, 0.0))
    pp, pm, sp, sm, split = miehe_response_2d(tension, MAT, TOL)
    @test pp ≈ (0.5lambda + mu) * 0.01^2
    @test pm ≈ 0.0 atol = 1.0e-14
    @test split.epsilon_plus + split.epsilon_minus ≈ tension
    @test sp + sm ≈ C0 ⊡ tension
    @test solve_rlm_quadratic(1.0, 0.0, -1.0, 1.0, TOL).q ≈ 1.0
    @test !solve_rlm_quadratic(1.0, 0.0, 1.0, 1.0, TOL).success
end

@testset "real-time affine BDF1 step" begin
    problem = build_rlm_problem(small_config(); grid = small_grid())
    update!(problem.ch_u, 0.0)
    state = RLMState(zeros(ndofs(problem.dh_u)), zeros(ndofs(problem.dh_d)), 1.0, 0.0)
    apply!(state.u, problem.ch_u)
    state.P = rlm_nonlinear_energy(problem, state.u, state.d)
    snapshot = copy(state)
    trial = compute_rlm_bdf1_trial(problem, state, 0.001)
    @test state.u == snapshot.u && state.d == snapshot.d && state.q == snapshot.q && state.P == snapshot.P
    free = free_dofs(problem.ch_u)
    equilibrium = problem.K_u * trial.u + trial.q * trial.n_u - problem.f_ext
    @test norm(equilibrium[free]) <= 1.0e-9
    phase = MAT.viscosity / problem.config.time.dt * problem.M_d * (trial.d - state.d) +
            problem.K_AT2 * trial.d + trial.q * trial.n_d
    @test norm(phase) <= 1.0e-9
    @test trial.c1 <= 1.0e-12
    @test trial.scalar_residual <= problem.config.tolerances.scalar_residual
    @test trial.energy_balance_residual <= problem.config.tolerances.energy_balance_rel
    @test isfinite(trial.external_work) && isfinite(trial.viscous_dissipation)
end

@testset "one update per physical time step and continuous q/P" begin
    result = solve_rlm_bdf1(build_rlm_problem(small_config(); grid = small_grid()))
    @test result.success && result.completed
    @test length(result.diagnostics) == 3
    @test all(d -> d.accepted, result.diagnostics)
    @test [d.step for d in result.diagnostics] == [0, 1, 2]
    @test [d.time for d in result.diagnostics] ≈ [0.0, 0.001, 0.002]
    @test all(d -> isfinite(d.cumulative_external_work), result.diagnostics)
    @test all(d -> isfinite(d.energy_balance_residual), result.diagnostics[2:end])
    @test !any(d -> occursin("reset", d.status), result.diagnostics)
    mktempdir() do directory
        path = joinpath(directory, "history.csv")
        write_rlm_time_history(path, result.diagnostics)
        @test length(readlines(path)) == 4
    end
end

@testset "time-dependent body and traction factors" begin
    config = small_config(displacement = history([0.0, 0.0, 0.0]), body = history([0.0, 0.5, 1.0]),
        traction = history([0.0, 0.25, 0.5]))
    problem = build_rlm_problem(config; grid = small_grid())
    update_rlm_external_force!(problem, 0.001)
    @test problem.f_ext ≈ 0.5 .* problem.f_body_reference .+ 0.25 .* problem.f_traction_reference
    @test norm(problem.f_ext) > 0.0
end

@testset "uniform viscous relaxation is first-order in time" begin
    function final_error(dt)
        final_time = 0.004
        times = collect(0.0:dt:final_time)
        zero_history = RLMPiecewiseLinearHistory(times, zeros(length(times)))
        config = small_config(displacement = zero_history, body = zero_history, traction = zero_history,
            initial_damage = 0.1, dt = dt, final_time = final_time)
        result = solve_rlm_bdf1(build_rlm_problem(config; grid = small_grid()))
        @test result.success
        exact = 0.1 * exp(-MAT.G_c / (MAT.viscosity * MAT.ell) * final_time)
        return maximum(abs.(result.state.d .- exact))
    end
    e1, e2 = final_error(0.001), final_error(0.0005)
    @test e2 < 0.65 * e1
end

@testset "transactional real-time rollback" begin
    base = small_config()
    config = RLMConfig(material = base.material, mesh = base.mesh, load = base.load, time = base.time,
        tolerances = RLMToleranceConfig(positive_root = 2.0), output = base.output)
    problem = build_rlm_problem(config; grid = small_grid())
    result = solve_rlm_bdf1(problem)
    @test !result.success && !result.completed
    @test !last(result.diagnostics).accepted
    @test last(result.diagnostics).step == 1
end
