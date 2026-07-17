module PffSAV

using LinearAlgebra
using Serialization
using SparseArrays
using Ferrite
using FerriteGmsh
using Tensors

# physics
include("physics/types.jl")
include("physics/constitutive.jl")
include("physics/energies.jl")

# fem setup
include("fem/setup.jl")
include("fem/assembly.jl")

# Miehe--RLM-PC--BDF1 (separate from the legacy staggered implementation)
include("rlm/config.jl")
include("rlm/miehe2d.jl")
include("rlm/assembly.jl")
include("rlm/scalar_solver.jl")
include("rlm/solver_bdf1.jl")

# solvers
include("solvers/staggered.jl")

# utilities
include("utils/utils_fun.jl")

export
    # --- struct ---
    MaterialParams,
    TensionSetup,
    RLMMaterialConfig,
    RLMMeshConfig,
    RLMLoadConfig,
    RLMTimeConfig,
    RLMToleranceConfig,
    RLMOutputConfig,
    RLMConfig,
    RLMState,
    RLMDiagnostic,
    RLMTrial,
    RLMResult,
    RLMProblem,
    # --- setups ---
    setup_tension,
    build_rlm_problem,
    # --- solvers ---
    solve_staggered,
    solve_rlm_bdf1,
    compute_rlm_bdf1_trial,
    solve_rlm_quadratic,
    scalar_equation_residual,
    # --- Miehe/RLM physics ---
    degradation,
    degradation_derivative,
    miehe_split_2d,
    miehe_response_2d,
    assemble_rlm_nonlinear_forces!,
    rlm_nonlinear_energy,
    rlm_raw_energy,
    rlm_proxy_energy,
    phase_field_metrics,
    write_rlm_diagnostics,
    # --- utilities ---
    compute_driving_force!,
    compute_reaction_forces,
    get_right_dofs
end
