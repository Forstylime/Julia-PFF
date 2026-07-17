# scripts/run_ch_rlm.jl

using PffSAV

function main()
    # 1. Parameter setup
    # Typical Cahn-Hilliard parameters for spinodal decomposition
    M = 1.0          # Mobility
    ϵ = 0.1         # Interface width
    α = 1e3          # RLM relaxation parameter (sufficiently large)
    Δt = 1e-4        # Time step size
    
    params = CHParams(M, ϵ, α, Δt)
    
    # 2. Setup grid and FEM structures
    # Using L-shape mesh from Gmsh
    println("Setting up Cahn-Hilliard L-shape case...")
    setup = setup_ch(msh_file = "data/mesh/l_shape_uni.msh")
    
    # 3. Solve using RLM-PC
    println("Starting RLM Cahn-Hilliard simulation...")
    # T_final=0.05 to see initial spinodal decomposition
    energies = solve_rlm_ch(setup, params; T_final=0.1, outdir="data/sims/ch_rlm")
    
    println("Simulation complete.")
end

main()
