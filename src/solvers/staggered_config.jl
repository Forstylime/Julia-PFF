"""Physical-time discretization used by the staggered phase-field solver."""
Base.@kwdef struct StaggeredTimeConfig
    n_steps::Int = 100
    dt::Union{Nothing,Float64} = nothing
    total_time::Float64 = 1.0
    time_points::Union{Nothing,Vector{Float64}} = nothing
end

"""Displacement-control history for the staggered solver."""
Base.@kwdef struct StaggeredLoadConfig
    load_history::Union{Nothing,Function} = nothing
    ramp_time::Union{Nothing,Float64} = nothing
end

"""Nonlinear iteration and phase-field evolution options."""
Base.@kwdef struct StaggeredSolverConfig
    tolerance::Float64 = 1.0e-5
    max_staggered_iterations::Int = 20
    max_newton_iterations::Int = 10
    enforce_irreversibility::Bool = true
    viscosity::Float64 = 0.0
end

"""VTK and console-output options for the staggered solver."""
Base.@kwdef struct StaggeredOutputConfig
    directory::Union{Nothing,String} = nothing
    write_vtk::Bool = true
    vtk_interval::Int = 5
    verbose::Bool = true
end

"""Complete, explicit configuration for a staggered phase-field simulation."""
Base.@kwdef struct StaggeredConfig
    time::StaggeredTimeConfig = StaggeredTimeConfig()
    load::StaggeredLoadConfig = StaggeredLoadConfig()
    solver::StaggeredSolverConfig = StaggeredSolverConfig()
    output::StaggeredOutputConfig = StaggeredOutputConfig()
end

"""Accepted or rejected staggered physical-time step diagnostic."""
Base.@kwdef struct StaggeredDiagnostic
    step::Int
    time::Float64
    dt::Float64
    load_factor::Float64
    imposed_displacement::Float64
    staggered_iterations::Int
    newton_iterations::Int
    displacement_residual::Float64
    converged::Bool
    damage_increment::Float64
    damage_min::Float64
    damage_max::Float64
    viscous_dissipation::Float64
    cumulative_viscous_dissipation::Float64
end

"""Complete result of a staggered simulation, including final raw DOF vectors."""
struct StaggeredResult
    success::Bool
    completed::Bool
    message::String
    config::StaggeredConfig
    times::Vector{Float64}
    displacements::Vector{Float64}
    reaction_forces::Vector{Float64}
    elastic_energies::Vector{Float64}
    surface_energies::Vector{Float64}
    diagnostics::Vector{StaggeredDiagnostic}
    u_final::Vector{Float64}
    d_final::Vector{Float64}
    cumulative_viscous_dissipation::Float64
    total_newton_iterations::Int
end
