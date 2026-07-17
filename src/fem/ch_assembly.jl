# src/fem/ch_assembly.jl

using Ferrite
using SparseArrays
using LinearAlgebra

"""
    assemble_CH_LHS!(K, dh, cv, params)

Assemble the constant left-hand side matrix for the Cahn-Hilliard equation.
"""
function assemble_CH_LHS!(
    K::SparseMatrixCSC,
    dh::DofHandler,
    cv::CellValues,
    params::CHParams
)
    assembler = start_assemble(K)
    ndofs_cell = ndofs_per_cell(dh)
    Ke = zeros(ndofs_cell, ndofs_cell)
    
    # Pre-extract DOF indices for :phi and :mu (e.g., 1:4 and 5:8)
    phi_dofs = Ferrite.dof_range(dh, :phi)
    mu_dofs = Ferrite.dof_range(dh, :mu)
    n_basefuncs = getnbasefunctions(cv)
    
    for cell in CellIterator(dh)
        reinit!(cv, cell)
        fill!(Ke, 0.0)
        
        for q_point in 1:getnquadpoints(cv)
            dΩ = getdetJdV(cv, q_point)
            
            # 1. Test function w (from phi_dofs)
            for i in 1:n_basefuncs
                w = shape_value(cv, q_point, i)
                ∇w = shape_gradient(cv, q_point, i)
                
                # Trial function ϕ (from phi_dofs)
                for j in 1:n_basefuncs
                    ϕ = shape_value(cv, q_point, j)
                    Ke[phi_dofs[i], phi_dofs[j]] += (1.0 / params.Δt) * ϕ * w * dΩ
                end
                
                # Trial function μ (from mu_dofs)
                for j in 1:n_basefuncs
                    ∇μ = shape_gradient(cv, q_point, j)
                    Ke[phi_dofs[i], mu_dofs[j]] += (params.M / 2.0) * (∇μ ⋅ ∇w) * dΩ
                end
            end
            
            # 2. Test function v (from mu_dofs)
            for i in 1:n_basefuncs
                v = shape_value(cv, q_point, i)
                ∇v = shape_gradient(cv, q_point, i)
                
                # Trial function ϕ (from phi_dofs)
                for j in 1:n_basefuncs
                    ∇ϕ = shape_gradient(cv, q_point, j)
                    Ke[mu_dofs[i], phi_dofs[j]] -= 0.5 * params.ϵ^2 * (∇ϕ ⋅ ∇v) * dΩ
                end
                
                # Trial function μ (from mu_dofs)
                for j in 1:n_basefuncs
                    μ = shape_value(cv, q_point, j)
                    Ke[mu_dofs[i], mu_dofs[j]] += 0.5 * μ * v * dΩ
                end
            end
        end
        assemble!(assembler, celldofs(cell), Ke)
    end
end

"""
    assemble_CH_RHS!(b1, b2, dh, cv, params, state)

Assemble the RHS vectors for the Cahn-Hilliard equation.
"""
function assemble_CH_RHS!(
    b1::Vector{Float64},
    b2::Vector{Float64},
    dh::DofHandler,
    cv::CellValues,
    params::CHParams,
    state::CHState
)
    ndofs_cell = ndofs_per_cell(dh)
    be1 = zeros(ndofs_cell)
    be2 = zeros(ndofs_cell)
    
    phi_dofs = Ferrite.dof_range(dh, :phi)
    mu_dofs = Ferrite.dof_range(dh, :mu)
    n_basefuncs = getnbasefunctions(cv)
    
    fill!(b1, 0.0)
    fill!(b2, 0.0)
    
    for cell in CellIterator(dh)
        reinit!(cv, cell)
        fill!(be1, 0.0)
        fill!(be2, 0.0)
        
        loc_dofs = celldofs(cell)
        ϕ_n_loc = state.ϕ_n[loc_dofs[phi_dofs]]
        ϕ_nm1_loc = state.ϕ_nm1[loc_dofs[phi_dofs]]
        μ_n_loc = state.μ_n[loc_dofs[mu_dofs]]
        
        for q_point in 1:getnquadpoints(cv)
            dΩ = getdetJdV(cv, q_point)
            
            ϕ_n_q = function_value(cv, q_point, ϕ_n_loc)
            ϕ_nm1_q = function_value(cv, q_point, ϕ_nm1_loc)
            μ_n_q = function_value(cv, q_point, μ_n_loc)
            
            ∇ϕ_n_q = function_gradient(cv, q_point, ϕ_n_loc)
            ∇μ_n_q = function_gradient(cv, q_point, μ_n_loc)
            
            # Extrapolate
            ϕ_bar_q = 1.5 * ϕ_n_q - 0.5 * ϕ_nm1_q
            f_bar_q = ch_f(ϕ_bar_q)
            
            # 1. Test function w (from phi_dofs)
            for i in 1:n_basefuncs
                w = shape_value(cv, q_point, i)
                ∇w = shape_gradient(cv, q_point, i)
                
                # b1 term for phi test function
                val1 = (ϕ_n_q / params.Δt) * w - (params.M / 2.0) * (∇μ_n_q ⋅ ∇w)
                be1[phi_dofs[i]] += val1 * dΩ
                
                # b2 term for phi test function is 0
            end
            
            # 2. Test function v (from mu_dofs)
            for i in 1:n_basefuncs
                v = shape_value(cv, q_point, i)
                ∇v = shape_gradient(cv, q_point, i)
                
                # b1 term for mu test function
                val1 = -0.5 * μ_n_q * v - 0.5 * params.ϵ^2 * (∇ϕ_n_q ⋅ ∇v) + 0.5 * state.q_n * f_bar_q * v
                be1[mu_dofs[i]] += val1 * dΩ
                
                # b2 term for mu test function
                val2 = 0.5 * f_bar_q * v
                be2[mu_dofs[i]] += val2 * dΩ
            end
        end
        assemble!(b1, loc_dofs, be1)
        assemble!(b2, loc_dofs, be2)
    end
end

"""
    compute_rlm_integrals(dh, cv, params, state, U_star_np1)

Computes the integrals B2 and C2 for the RLM-PC algorithm.
"""
function compute_rlm_integrals(
    dh::DofHandler,
    cv::CellValues,
    params::CHParams,
    state::CHState,
    U_star_np1::Vector{Float64}
)
    phi_dofs = Ferrite.dof_range(dh, :phi)
    
    B2 = 0.0
    C2_int = 0.0
    
    for cell in CellIterator(dh)
        reinit!(cv, cell)
        
        loc_dofs = celldofs(cell)
        ϕ_n_loc = state.ϕ_n[loc_dofs[phi_dofs]]
        ϕ_nm1_loc = state.ϕ_nm1[loc_dofs[phi_dofs]]
        ϕ_star_n_loc = state.ϕ_star_n[loc_dofs[phi_dofs]]
        ϕ_star_np1_loc = U_star_np1[loc_dofs[phi_dofs]]
        
        for q_point in 1:getnquadpoints(cv)
            dΩ = getdetJdV(cv, q_point)
            
            ϕ_n_q = function_value(cv, q_point, ϕ_n_loc)
            ϕ_nm1_q = function_value(cv, q_point, ϕ_nm1_loc)
            ϕ_star_n_q = function_value(cv, q_point, ϕ_star_n_loc)
            ϕ_star_np1_q = function_value(cv, q_point, ϕ_star_np1_loc)
            
            ϕ_bar_q = 1.5 * ϕ_n_q - 0.5 * ϕ_nm1_q
            f_bar_q = ch_f(ϕ_bar_q)
            
            # B2 integral term
            B2 += -0.5 * f_bar_q * (ϕ_star_np1_q - ϕ_n_q) * dΩ
            
            # C2 integral term: ∫ ( F(ϕ^{n+1,*}) - F(ϕ^{n,*}) ) dx
            C2_int += (ch_F(ϕ_star_np1_q) - ch_F(ϕ_star_n_q)) * dΩ
        end
    end
    
    C2 = C2_int + state.q_n * B2 - params.α * (state.q_n)^2
    
    return B2, C2
end
