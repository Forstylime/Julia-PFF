using CairoMakie
using JLD2
using PffSAV

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const PLOT_DIRECTORY = joinpath(PROJECT_ROOT, "data", "plots")
const RESULT_DIRECTORY = joinpath(PROJECT_ROOT, "data", "jld2")

# Set ETA > 0 to run the BDF1 viscous phase-field model.  Keeping ETA == 0
# reproduces the original quasi-static staggered calculation.
const ETA = 0.0001
const N_STEPS = 1000
const DT = 0.001 # physical BDF1 time step; set TIME_POINTS = nothing to use it
const TOTAL_TIME = N_STEPS * DT
const TIME_POINTS = nothing # or e.g. [0.0, 0.02, 0.1, 1.0]
const MODE_NAME = ETA > 0 ? "staggered_bdf1" : "staggered"
const VTK_DIRECTORY = joinpath(PROJECT_ROOT, "data", "sims", MODE_NAME)

"""Load factor evaluated at physical time t; replace for a different history."""
load_history(t) = t / TOTAL_TIME

function plot_results(setup, displacement, force, _elastic_energy, surface_energy)
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

    load_figure = Figure(size = (600, 400))
    load_axis = Axis(
        load_figure[1, 1];
        xlabel = L"\bar{u}~\mathrm{[mm]}",
        ylabel = L"F_{\mathrm{reaction}}~\mathrm{[N]}",
        title = "Displacement – Reaction Force",
        xgridvisible = true,
        ygridvisible = true,
        xgridcolor = :lightgray,
        ygridcolor = :lightgray,
    )
    lines!(load_axis, plotted_displacement, plotted_force; color = :red, linewidth = 2)
    load_plot_path = joinpath(PLOT_DIRECTORY, "load_displacement_staggered.png")
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
    energy_plot_path = joinpath(PLOT_DIRECTORY, "energy_evolution_staggered.png")
    save(energy_plot_path, energy_figure)

    println("Load-displacement curve saved to $load_plot_path")
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
        msh_file = joinpath(PROJECT_ROOT, "data", "mesh", "l_shape.msh"),
        final_displacement = -0.8,
        fixed_face = "top",
        dir = 2,
    )

    solver_dt = isnothing(TIME_POINTS) ? DT : nothing
    result = solve_staggered(
        setup,
        material;
        n_steps = N_STEPS,
        max_iter = 2_000,
        enforce_irreversibility = false,
        eta = ETA,
        dt = solver_dt,
        total_time = TOTAL_TIME,
        time_points = TIME_POINTS,
        load_history = load_history,
        return_diagnostics = true,
        output_directory = VTK_DIRECTORY,
    )
    displacement = result.displacements
    force = result.reaction_forces
    elastic_energy = result.elastic_energies
    surface_energy = result.surface_energies

    mkpath(PLOT_DIRECTORY)
    mkpath(RESULT_DIRECTORY)
    plot_results(setup, displacement, force, elastic_energy, surface_energy)

    result_path = joinpath(RESULT_DIRECTORY, "$(MODE_NAME)_results.jld2")
    jldsave(
        result_path;
        displacement,
        force,
        elastic_energy,
        surface_energy,
        physical_time = result.times,
        diagnostics = result.diagnostics,
        cumulative_viscous_dissipation = result.cumulative_viscous_dissipation,
        eta = ETA,
        dt = solver_dt,
        total_time = TOTAL_TIME,
        enforce_irreversibility = true,
    )
    println("Staggered simulation data saved to $result_path")
    return nothing
end

main()
