using CairoMakie
using JLD2
using PffSAV

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const PLOT_DIRECTORY = joinpath(PROJECT_ROOT, "data", "plots")
const RESULT_DIRECTORY = joinpath(PROJECT_ROOT, "data", "jld2")

# Set ETA > 0 to run the BDF1 viscous phase-field model.  Keeping ETA == 0
# reproduces the original quasi-static staggered calculation.
const ETA = 0.1
const DT = 0.005 # physical BDF1 time step; set TIME_POINTS = nothing to use it
const TOTAL_TIME = 100.0 # total physical time for the BDF1 simulation
const N_STEPS = Int(TOTAL_TIME / DT) # number of time steps for the staggered simulation
const RAMP_TIME = 50.0 # linear ramp ends here; hold the final displacement afterwards
const TIME_POINTS = nothing # or e.g. [0.0, 0.02, 0.1, 1.0]
const MODE_NAME = ETA > 0 ? "staggered_bdf1" : "staggered"
const VTK_DIRECTORY = joinpath(PROJECT_ROOT, "data", "sims", MODE_NAME)
const STAGGERED_CONFIG = StaggeredConfig(
    time = StaggeredTimeConfig(
        n_steps = N_STEPS,
        dt = isnothing(TIME_POINTS) ? DT : nothing,
        total_time = TOTAL_TIME,
        time_points = TIME_POINTS,
    ),
    load = StaggeredLoadConfig(ramp_time = RAMP_TIME),
    solver = StaggeredSolverConfig(
        tolerance = 1.0e-5,
        max_staggered_iterations = 1_000,
        max_newton_iterations = 20,
        enforce_irreversibility = false,
        viscosity = ETA,
    ),
    output = StaggeredOutputConfig(
        directory = VTK_DIRECTORY,
        write_vtk = true,
        vtk_interval = Int(N_STEPS ÷ 50),
        verbose = true,
    ),
)

function plot_results(setup, displacement, force, _elastic_energy, surface_energy, times)
    direction_sign = sign(setup.final_displacement)
    plotted_displacement = direction_sign .* displacement
    plotted_force = direction_sign .* force

    peak_index = argmax(plotted_force)
    peak_force = plotted_force[peak_index]
    peak_displacement = plotted_displacement[peak_index]
    println(
        "Peak load: F_max = $(round(peak_force; digits = 4)) N ",
        "@ ū = $(round(peak_displacement; digits = 4)) mm",
    )

    load_figure = Figure(size = (900, 400))
    load_axis_u = Axis(
        load_figure[1, 1];
        xlabel = L"\bar{u}~\mathrm{[mm]}",
        ylabel = L"F_{\mathrm{reaction}}~\mathrm{[N]}",
        title = "Displacement – Reaction Force",
        xgridvisible = true,
        ygridvisible = true,
        xgridcolor = :lightgray,
        ygridcolor = :lightgray,
    )
    lines!(load_axis_u, plotted_displacement, plotted_force; color = :red, linewidth = 2)

    load_axis_t = Axis(
        load_figure[1, 2];
        xlabel = L"Time~\mathrm{[s]}",
        ylabel = L"F_{\mathrm{reaction}}~\mathrm{[N]}",
        title = "Time – Reaction Force",
        xgridvisible = true,
        ygridvisible = true,
        xgridcolor = :lightgray,
        ygridcolor = :lightgray,
    )
    lines!(load_axis_t, times, plotted_force; color = :blue, linewidth = 2)

    load_plot_path = joinpath(PLOT_DIRECTORY, "load_displacement_staggered_2.png")
    save(load_plot_path, load_figure)

    energy_figure = Figure(size = (600, 400))
    energy_axis = Axis(
        energy_figure[1, 1];
        xlabel = L"\bar{u}~\mathrm{[mm]}",
        ylabel = L"\mathcal{G}_f~\mathrm{[N\,mm]}",
        title = L"\mathcal{G}_f\ \mathrm{(surface)\ Evolution}",
        xgridvisible = true,
        ygridvisible = true,
        xgridcolor = :lightgray,
        ygridcolor = :lightgray,
    )
    lines!(
        energy_axis,
        plotted_displacement,
        surface_energy;
        color = :orange,
        label = L"\mathcal{G}_f\ \mathrm{(surface)}",
        linewidth = 2,
    )
    energy_plot_path = joinpath(PLOT_DIRECTORY, "energy_evolution_staggered_2.png")
    save(energy_plot_path, energy_figure)

    println("Load curve saved to $load_plot_path")
    println("Energy evolution curve saved to $energy_plot_path")
    return nothing
end

function main()
    material = MaterialParams(;
        dim = 2,
        E = 25_840.0,
        ν = 0.18,
        gc = 0.65,
        l = 10.0,
        k = 1.0e-8,
    )
    setup = setup_tension(;
        msh_file = joinpath(PROJECT_ROOT, "data", "mesh", "l_shape_2.msh"),
        final_displacement = -0.8,
        fixed_face = "top",
        dir = 2,
    )

    result = solve_staggered(setup, material, STAGGERED_CONFIG)
    println(result.message)
    result.success || error("Staggered simulation did not complete; no comparison data was written.")
    displacement = result.displacements
    force = result.reaction_forces
    elastic_energy = result.elastic_energies
    surface_energy = result.surface_energies
    u_final = result.u_final
    d_final = result.d_final
    imposed_displacement_final = last(displacement)

    mkpath(PLOT_DIRECTORY)
    mkpath(RESULT_DIRECTORY)
    times = result.times
    plot_results(setup, displacement, force, elastic_energy, surface_energy, times)

    result_path = joinpath(RESULT_DIRECTORY, "$(MODE_NAME)_results_2.jld2")
    jldsave(
        result_path;
        displacement,
        imposed_displacement_final,
        u_final,
        d_final,
        force,
        elastic_energy,
        surface_energy,
        physical_time = times,
        diagnostics = result.diagnostics,
        cumulative_viscous_dissipation = result.cumulative_viscous_dissipation,
        eta = STAGGERED_CONFIG.solver.viscosity,
        dt = STAGGERED_CONFIG.time.dt,
        total_time = STAGGERED_CONFIG.time.total_time,
        ramp_time = STAGGERED_CONFIG.load.ramp_time,
        final_displacement_amplitude = setup.final_displacement,
        enforce_irreversibility = STAGGERED_CONFIG.solver.enforce_irreversibility,
    )
    println("Staggered simulation data saved to $result_path")
    return nothing
end

main()
