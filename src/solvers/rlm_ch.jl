# src/solvers/rlm_ch.jl

using Ferrite
using LinearAlgebra
using SparseArrays

"""
    solve_rlm_ch(setup::CHSetup, params::CHParams; T_final=1.0)

Executes the RLM-PC algorithm for the Cahn-Hilliard equation.
"""
function solve_rlm_ch(
    setup::CHSetup, params::CHParams;
    T_final = 1.0,
    outdir = "data/sims/ch_rlm"
)
    mkpath(outdir)
    
    grid = setup.grid
    dh = setup.dh
    ch = setup.ch
    
    ndofs_total = ndofs(dh)
    
    # 1. Initialize CellValues
    qr = QuadratureRule{RefQuadrilateral}(2)
    cv = CellValues(qr, Lagrange{RefQuadrilateral, 1}())
    
    # 2. Initialize State
    state = CHState(ndofs_total)
    
    # Initial Condition: Random noise for spinodal decomposition
    # Only initialize ϕ, μ starts as 0
    phi_dofs = Ferrite.dof_range(dh, :phi)
    for cell in CellIterator(dh)
        loc_dofs = celldofs(cell)
        for i in phi_dofs
            global_dof = loc_dofs[i]
            # Assign random noise between -0.05 and 0.05
            state.ϕ_n[global_dof] = 0.1 * rand() - 0.05
            state.ϕ_nm1[global_dof] = state.ϕ_n[global_dof] # Initially phi^{n-1} = phi^n
        end
    end
    
    # Initialize historical states
    state.ϕ_star_n .= state.ϕ_n
    
    # 3. Assemble and Factorize constant LHS
    K = allocate_matrix(dh)
    assemble_CH_LHS!(K, dh, cv, params)
    
    apply!(K, ch)
    K_fact = lu(K)
    
    # RHS arrays
    b1 = zeros(ndofs_total)
    b2 = zeros(ndofs_total)
    
    # 4. Time stepping loop
    n_steps = Int(round(T_final / params.Δt))
    println("LHS pre-factorization complete. Starting RLM-PC time stepping loop. Total steps: $n_steps")
    
    # Relaxed energy tracking
    energies = Float64[]
    
    # Save initial state
    VTKGridFile(joinpath(outdir, "ch_step_0"), dh) do vtk
        write_solution(vtk, dh, state.ϕ_n)
    end
    
    for step in 1:n_steps
        # Step 1: Assemble RHS and solve for U1, U2
        assemble_CH_RHS!(b1, b2, dh, cv, params, state)
        
        apply!(b1, ch)
        apply_zero!(b2, ch)
        
        U1 = K_fact \ b1
        U2 = K_fact \ b2
        
        # Step 2: Prediction state and RLM coefficients
        U_star_np1 = U1 .+ U2
        
        B2, C2 = compute_rlm_integrals(dh, cv, params, state, U_star_np1)
        A2 = params.α
        
        # Solve quadratic equation A2 * q^2 + B2 * q + C2 = 0
        discriminant = B2^2 - 4 * A2 * C2
        if discriminant < 0
            @warn "Discriminant is negative ($discriminant). Setting to 0."
            discriminant = 0.0
        end
        
        q_p = (-B2 + sqrt(discriminant)) / (2 * A2)
        q_m = (-B2 - sqrt(discriminant)) / (2 * A2)
        
        # Root Selection Criterion
        if abs(q_p - state.q_n) < abs(q_m - state.q_n)
            q_np1 = q_p
        else
            q_np1 = q_m
        end
        
        # Step 3: Update physical variables
        U_np1 = U1 .+ q_np1 .* U2
        
        # Extract phi and mu from U_np1 to update state
        # Ferrite's dof structure interleaves phi and mu, so we use U_np1 directly as it has the same layout
        state.ϕ_nm1 .= state.ϕ_n
        state.ϕ_n .= U_np1
        state.μ_n .= U_np1 # Both are technically using the full U_np1, but our loops only read their specific dofs
        
        # Historical tracking for next step
        state.ϕ_star_n .= U_star_np1
        state.q_n = q_np1
        
        # Calculate relaxed energy (simplified calculation for tracking)
        # Note: True tracking requires assembling the physical energy. Here we use C2 for validation.
        push!(energies, q_np1-1)
        
        if step % 10 == 0 || step == n_steps
            println("Step $step / $n_steps | q-1: $(round(q_np1-1, digits=9))")
            VTKGridFile(joinpath(outdir, "ch_step_$step"), dh) do vtk
                write_solution(vtk, dh, U_np1)
            end
        end
    end
    
    println("Simulation finished! Results saved to $outdir")
    return energies
end
