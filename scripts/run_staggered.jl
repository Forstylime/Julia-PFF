using CairoMakie
using JLD2
using PffSAV

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const PLOT_DIRECTORY = joinpath(PROJECT_ROOT, "data", "plots")
const RESULT_DIRECTORY = joinpath(PROJECT_ROOT, "data", "jld2")
const VTK_DIRECTORY = joinpath(PROJECT_ROOT, "data", "sims", "staggered")

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

    displacement, force, elastic_energy, surface_energy = solve_staggered(
        setup,
        material;
        n_steps = 100,
        max_iter = 2_000,
        enforce_irreversibility = true,
        output_directory = VTK_DIRECTORY,
    )

    mkpath(PLOT_DIRECTORY)
    mkpath(RESULT_DIRECTORY)
    plot_results(setup, displacement, force, elastic_energy, surface_energy)

    result_path = joinpath(RESULT_DIRECTORY, "staggered_results.jld2")
    jldsave(
        result_path;
        displacement,
        force,
        elastic_energy,
        surface_energy,
    )
    println("Staggered simulation data saved to $result_path")
    return nothing
end

main()
