using CairoMakie
using JLD2
using LinearAlgebra

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const RESULT_DIRECTORY = joinpath(PROJECT_ROOT, "data", "jld2")
const PLOT_DIRECTORY = joinpath(PROJECT_ROOT, "data", "plots")
const STAGGERED_RESULT = joinpath(RESULT_DIRECTORY, "staggered_bdf1_results_2.jld2")
const RLM_RESULT = joinpath(RESULT_DIRECTORY, "rlm_bdf1_results.jld2")

function require_dataset(data, name::AbstractString, path::AbstractString)
    haskey(data, name) || throw(ArgumentError("$path does not contain dataset '$name'"))
    return data[name]
end

function main()
    isfile(STAGGERED_RESULT) || throw(ArgumentError("Staggered result not found: $STAGGERED_RESULT"))
    isfile(RLM_RESULT) || throw(ArgumentError("RLM-BDF1 result not found: $RLM_RESULT"))

    staggered = load(STAGGERED_RESULT)
    rlm = load(RLM_RESULT)
    u_staggered = require_dataset(staggered, "u_final", STAGGERED_RESULT)
    u_rlm = require_dataset(rlm, "u_final", RLM_RESULT)
    length(u_staggered) == length(u_rlm) || throw(ArgumentError(
        "The two final displacement vectors have different lengths " *
        "($(length(u_staggered)) and $(length(u_rlm)); use the same mesh and DOF layout.",
    ))

    imposed_staggered = require_dataset(staggered, "imposed_displacement_final", STAGGERED_RESULT)
    imposed_rlm = require_dataset(rlm, "imposed_displacement_final", RLM_RESULT)
    staggered_dt = require_dataset(staggered, "dt", STAGGERED_RESULT)
    rlm_dt = require_dataset(rlm, "dt", RLM_RESULT)
    if staggered_dt isa Real && rlm_dt isa Real
        isapprox(staggered_dt, rlm_dt) || @warn(
            "The time steps differ; the final fields are not a time-step convergence comparison.",
            staggered_dt,
            rlm_dt,
        )
    else
        @warn "At least one result does not use a constant time step; compare the physical-time grids separately."
    end
    imposed_difference = imposed_rlm - imposed_staggered
    u_difference = u_rlm - u_staggered
    u_difference_norm = norm(u_difference)
    u_relative_difference = u_difference_norm / max(norm(u_staggered), eps(Float64))
    dof = eachindex(u_difference)

    mkpath(PLOT_DIRECTORY)
    figure = Figure(size = (1100, 800))
    values_axis = Axis(
        figure[1, 1];
        xlabel = "displacement DOF index",
        ylabel = "u [mm]",
        title = "Final displacement vector",
    )
    lines!(values_axis, dof, u_staggered; label = "Staggered", linewidth = 2)
    lines!(values_axis, dof, u_rlm; label = "RLM-BDF1", linewidth = 2, linestyle = :dash)
    axislegend(values_axis, position = :rb)

    difference_axis = Axis(
        figure[2, 1];
        xlabel = "displacement DOF index",
        ylabel = "u_RLM - u_Staggered [mm]",
        title = "Final displacement-vector difference",
    )
    lines!(difference_axis, dof, u_difference; color = :firebrick, linewidth = 1.5)
    hlines!(difference_axis, [0.0]; color = :gray50, linestyle = :dash)
    Label(
        figure[3, 1],
        "Imposed displacement difference: $(imposed_difference) mm    " *
        "||u_RLM - u_Staggered||₂: $(u_difference_norm)    " *
        "relative difference: $(u_relative_difference)",
        tellwidth = false,
    )

    plot_path = joinpath(PLOT_DIRECTORY, "staggered_rlm_final_displacement_difference.png")
    save(plot_path, figure)
    comparison_path = joinpath(RESULT_DIRECTORY, "staggered_rlm_final_displacement_difference.jld2")
    jldsave(
        comparison_path;
        imposed_staggered,
        imposed_rlm,
        imposed_difference,
        u_staggered,
        u_rlm,
        u_difference,
        u_difference_norm,
        u_relative_difference,
    )

    println("Imposed displacement difference: $imposed_difference mm")
    println("Final displacement-vector L2 difference: $u_difference_norm")
    println("Final displacement-vector relative difference: $u_relative_difference")
    println("Comparison plot saved to $plot_path")
    println("Comparison data saved to $comparison_path")
    return nothing
end

main()
