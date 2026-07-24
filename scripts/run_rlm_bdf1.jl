using PffSAV
using CairoMakie
using JLD2

project_root = normpath(joinpath(@__DIR__, ".."))
result_directory = joinpath(project_root, "data", "jld2")
ramp_time = 50.0
# Hold after the ramp; `final_time` must be strictly larger than `ramp_time`
# because the piecewise-linear history requires strictly increasing knots.
final_time = 100.0
dt = 0.001

config = RLMConfig(
    material = RLMMaterialConfig(E = 25_840.0, nu = 0.18, G_c = 0.65, ell = 10.0,
        kappa = 1.0e-8, viscosity = 0.1),
    mesh = RLMMeshConfig(path = joinpath(project_root, "data", "mesh", "l_shape_2.msh"), quadrature_order = 2),
    load = RLMLoadConfig(fixed_boundary = "top", loaded_boundary = "right", component = 2,
        overlap_policy = :loaded, displacement_amplitude = -0.8,
        displacement_history = RLMPiecewiseLinearHistory([0.0, ramp_time, final_time], [0.0, 1.0, 1.0]),
        body_force_history = RLMPiecewiseLinearHistory([0.0, final_time], [0.0, 0.0]),
        traction_history = RLMPiecewiseLinearHistory([0.0, final_time], [0.0, 0.0])),
    time = RLMTimeConfig(final_time = final_time, dt = dt, alpha = 1_000.0),
    output = RLMOutputConfig(directory = joinpath(project_root, "data", "sims", "rlm_bdf1_0_01"),
        write_csv = true, write_vtk = true, vtk_every_time_step = 1000, verbose = false),
)

result = solve_rlm_bdf1(build_rlm_problem(config))
println(result.message)
result.success || exit(1)
history = filter(diagnostic -> diagnostic.accepted, result.diagnostics)

times = getfield.(history, :time)
imposed_displacements = getfield.(history, :displacement)
displacements = -imposed_displacements # positive magnitude for the plotting convention below
reactions_rlm = -getfield.(history, :reaction_rlm)
reactions_phys = -getfield.(history, :reaction_phys)
reaction_difference = reactions_phys .- reactions_rlm
damages = getfield.(history, :max_d)
relaxed_energy = getfield.(history, :relaxed_internal_energy)
work = getfield.(history, :cumulative_external_work)
dissipation = getfield.(history, :cumulative_viscous_dissipation) .+ getfield.(history, :cumulative_numerical_dissipation)
q_minus_one = getfield.(history, :q_minus_one)
energy_plus_dissipation_minus_work = relaxed_energy .+ dissipation .- work
u_final = copy(result.state.u)
d_final = copy(result.state.d)
imposed_displacement_final = last(imposed_displacements)

mkpath(result_directory)
result_path = joinpath(result_directory, "rlm_bdf1_results.jld2")
jldsave(
    result_path;
    displacement = imposed_displacements,
    imposed_displacement_final,
    u_final,
    d_final,
    physical_time = times,
    reaction_rlm = reactions_rlm,
    reaction_phys = reactions_phys,
    diagnostics = result.diagnostics,
    dt,
    total_time = final_time,
    ramp_time,
    eta = config.material.viscosity,
    final_displacement_amplitude = config.load.displacement_amplitude,
)
println("RLM-BDF1 simulation data saved to $result_path")

figure = Figure(size = (1300, 1500))
ax1 = Axis(figure[1, 1]; xlabel = "imposed displacement [mm]", ylabel = "reaction force [N]")
lines!(ax1, displacements, reactions_rlm; label = "Reactions_RLM", linewidth = 2)
lines!(ax1, displacements, reactions_phys; label = "Reactions_phys", linewidth = 2)
axislegend(ax1, position = :lt)
ax2 = Axis(figure[1, 2]; xlabel = "time", ylabel = "reaction force [N]")
lines!(ax2, times, reactions_rlm; label = "Reactions_RLM", linewidth = 2, color = :green)
lines!(ax2, times, reactions_phys; label = "Reactions_phys", linestyle = :dash, linewidth = 2, color = :red)
lines!(ax2, times, reaction_difference; label = "phys - RLM", linestyle = :dash, linewidth = 2, color = :blue)
axislegend(ax2, position = :lt)
ax3 = Axis(figure[2, 1]; xlabel = "time", ylabel = "max damage")
lines!(ax3, times, damages; linewidth = 2)
ax4 = Axis(figure[2, 2]; xlabel = "time", ylabel = "energy")
lines!(ax4, times, relaxed_energy; label = "relaxed internal", linewidth = 2)
lines!(ax4, times, work; label = "external work", linewidth = 2)
lines!(ax4, times, dissipation; label = "cumulative dissipation", linewidth = 2)
axislegend(ax4, position = :lt)
ax5 = Axis(figure[3, 1]; xlabel = "time", ylabel = "q - 1")
lines!(ax5, times, q_minus_one; linewidth = 2)
hlines!(ax5, [0.0]; color = :gray50, linestyle = :dash)
ax6 = Axis(figure[3, 2]; xlabel = "time", ylabel = "energy")
lines!(ax6, times, energy_plus_dissipation_minus_work;
    label = "relaxed + dissipation - work", linewidth = 2)
hlines!(ax6, [first(relaxed_energy)]; color = :gray50, linestyle = :dash,
    label = "initial relaxed energy")
axislegend(ax6, position = :lb)
save(joinpath(project_root, "data", "plots", "rlm_bdf1_0_01.png"), figure)

figure_2 = Figure(size = (900, 700))
ax7 = Axis(figure_2[1, 1]; xlabel = "time", ylabel = "reaction force difference [N]")
lines!(ax7, times, reaction_difference; label = "phys - RLM", linestyle = :solid, linewidth = 2, color = :blue)
axislegend(ax7, position = :rt)
save(joinpath(project_root, "data", "plots", "rlm_bdf1_0_01_reaction_difference.png"), figure_2)
