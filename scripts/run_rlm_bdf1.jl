using PffSAV
using CairoMakie

project_root = normpath(joinpath(@__DIR__, ".."))
ramp_time = 1.0
# Hold after the ramp; `final_time` must be strictly larger than `ramp_time`
# because the piecewise-linear history requires strictly increasing knots.
final_time = 1.3
dt = 0.00001

config = RLMConfig(
    material = RLMMaterialConfig(E = 25_840.0, nu = 0.18, G_c = 0.65, ell = 10.0,
        kappa = 1.0e-8, viscosity = 0.001),
    mesh = RLMMeshConfig(path = joinpath(project_root, "data", "mesh", "l_shape.msh"), quadrature_order = 2),
    load = RLMLoadConfig(fixed_boundary = "top", loaded_boundary = "right", component = 2,
        overlap_policy = :loaded, displacement_amplitude = -0.8,
        displacement_history = RLMPiecewiseLinearHistory([0.0, ramp_time, final_time], [0.0, 1.0, 1.0]),
        body_force_history = RLMPiecewiseLinearHistory([0.0, final_time], [0.0, 0.0]),
        traction_history = RLMPiecewiseLinearHistory([0.0, final_time], [0.0, 0.0])),
    time = RLMTimeConfig(final_time = final_time, dt = dt, alpha = 1_000.0),
    output = RLMOutputConfig(directory = joinpath(project_root, "data", "sims", "rlm_bdf1_realtime"),
        write_csv = true, write_vtk = true, vtk_every_time_step = 1000, verbose = false),
)

result = solve_rlm_bdf1(build_rlm_problem(config))
println(result.message)
result.success || exit(1)
history = filter(diagnostic -> diagnostic.accepted, result.diagnostics)

times = getfield.(history, :time)
displacements = -getfield.(history, :displacement)
reactions = -getfield.(history, :reaction_force)
damages = getfield.(history, :max_d)
proxy = getfield.(history, :proxy_energy)
work = getfield.(history, :cumulative_external_work)
dissipation = getfield.(history, :cumulative_viscous_dissipation) .+ getfield.(history, :cumulative_numerical_dissipation)

figure = Figure(size = (900, 700))
ax1 = Axis(figure[1, 1]; xlabel = "imposed displacement [mm]", ylabel = "reaction force [N]")
lines!(ax1, displacements, reactions; linewidth = 2)
ax2 = Axis(figure[1, 2]; xlabel = "time", ylabel = "reaction force [N]")
lines!(ax2, times, reactions; linewidth = 2)
ax3 = Axis(figure[2, 1]; xlabel = "time", ylabel = "max damage")
lines!(ax3, times, damages; linewidth = 2)
ax4 = Axis(figure[2, 2]; xlabel = "time", ylabel = "energy")
lines!(ax4, times, proxy; label = "proxy internal", linewidth = 2)
lines!(ax4, times, work; label = "external work", linewidth = 2)
lines!(ax4, times, dissipation; label = "cumulative dissipation", linewidth = 2)
axislegend(ax4)
save(joinpath(project_root, "data", "plots", "rlm_bdf1.png"), figure)
