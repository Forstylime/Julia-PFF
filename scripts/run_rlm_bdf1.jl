using PffSAV

project_root = normpath(joinpath(@__DIR__, ".."))

config = RLMConfig(
    material = RLMMaterialConfig(
        E = 25_840.0,
        nu = 0.18,
        G_c = 0.65,
        ell = 10.0,
        kappa = 1.0e-6,
        mobility = 20.0,
    ),
    mesh = RLMMeshConfig(
        path = joinpath(project_root, "l_shape.msh"),
        quadrature_order = 2,
    ),
    load = RLMLoadConfig(
        fixed_boundary = "top",
        loaded_boundary = "right",
        component = 2,
        overlap_policy = :loaded,
        final_displacement = -0.01,
        load_steps = 20,
        initial_damage = 0.0,
        body_force = (0.0, 0.0),
        traction_boundary = nothing,
        traction = (0.0, 0.0),
    ),
    time = RLMTimeConfig(
        dt = 1.0e-3,
        alpha = 1.0,
        min_relax_steps = 1,
        max_relax_steps = 100,
    ),
    tolerances = RLMToleranceConfig(
        principal_zero_abs = 1.0e-14,
        principal_zero_rel = 1.0e-12,
        repeated_eigen_abs = 1.0e-14,
        repeated_eigen_rel = 1.0e-12,
        discriminant_abs = 0.0,
        discriminant_rel = 1.0e-12,
        scalar_residual = 1.0e-10,
        phase = 1.0e-6,
        q = 1.0e-6,
    ),
    output = RLMOutputConfig(
        directory = joinpath(project_root, "data", "sims", "rlm_bdf1"),
        write_csv = true,
        write_vtk = true,
        vtk_every_load_step = 1,
        verbose = true,
    ),
)

problem = build_rlm_problem(config)
result = solve_rlm_bdf1(problem)
println(result.message)
result.success || exit(1)
