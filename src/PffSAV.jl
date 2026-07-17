module PffSAV

using LinearAlgebra
using Serialization
using SparseArrays
using Ferrite
using FerriteGmsh
using Tensors

# physics
include("physics/types.jl")
include("physics/ch_types.jl")
include("physics/constitutive.jl")
include("physics/energies.jl")
include("physics/ch_energies.jl")

# fem setup
include("fem/setup.jl")
include("fem/ch_setup.jl")
include("fem/assembly.jl")
include("fem/ch_assembly.jl")

# solvers
include("solvers/sav.jl")
include("solvers/sav_quasistatic.jl")
include("solvers/staggered.jl")
include("solvers/rlm_ch.jl")

# utilities
include("utils/utils_fun.jl")

export
    # --- struct ---
    MaterialParams,
    NumericalParams,
    SimulationState,
    TensionSetup,
    TensionSetupSAV,
    CHParams,
    CHState,
    CHSetup,
    # --- setups ---
    setup_tension,
    setup_tension_sav,
    setup_ch,
    # --- solvers ---
    solve_sav,
    solve_sav_quasistatic,
    solve_staggered,
    solve_rlm_ch,
    # --- state ---
    update_states!,
    # --- utilities ---
    compute_driving_force!,
    compute_sav_scalars,
    compute_reaction_forces,
    get_right_dofs,
    # --- constitutive ---
    evaluate_damaged_stress,
    tensile_energy_density,
    elastic_energy_density_tensile,
    ch_F,
    ch_f
end
