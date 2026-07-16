# src/fem/assembly_rlm.jl

using Ferrite
using LinearAlgebra
using SparseArrays

"""
    assemble_RLM_LHS!(K_u, K_d, dh_u, dh_d, cv_u, cv_d, mat, rlm_params)

装配 RLM 算法中四个分支共享的全域常数 LHS 矩阵。
"""
function assemble_RLM_LHS!(
    K_u::SparseMatrixCSC, K_d::SparseMatrixCSC,
    dh_u::DofHandler, dh_d::DofHandler,
    cv_u::CellValues, cv_d::CellValues,
    mat::MaterialParams, rlm_params::RLMParams
)
    n_base_u = getnbasefunctions(cv_u)
    n_base_d = getnbasefunctions(cv_d)

    ke_u = zeros(n_base_u, n_base_u)
    ke_d = zeros(n_base_d, n_base_d)

    assembler_u = start_assemble(K_u)
    assembler_d = start_assemble(K_d)

    C0 = mat.C0
    coef_d = 1.5 / rlm_params.Δt + rlm_params.M * mat.gc / mat.l
    coef_grad = rlm_params.M * mat.gc * mat.l

    for (cell_u, cell_d) in zip(CellIterator(dh_u), CellIterator(dh_d))
        reinit!(cv_u, cell_u)
        reinit!(cv_d, cell_d)

        fill!(ke_u, 0.0)
        fill!(ke_d, 0.0)

        for qp in 1:getnquadpoints(cv_u)
            dΩ = getdetJdV(cv_u, qp)

            for i in 1:n_base_u
                δε_i = shape_symmetric_gradient(cv_u, qp, i)
                for j in 1:n_base_u
                    δε_j = shape_symmetric_gradient(cv_u, qp, j)
                    ke_u[i, j] += (δε_i ⊡ (C0 ⊡ δε_j)) * dΩ
                end
            end

            for i in 1:n_base_d
                δd_i = shape_value(cv_d, qp, i)
                ∇δd_i = shape_gradient(cv_d, qp, i)
                for j in 1:n_base_d
                    δd_j = shape_value(cv_d, qp, j)
                    ∇δd_j = shape_gradient(cv_d, qp, j)
                    ke_d[i, j] += (coef_d * δd_i * δd_j + coef_grad * (∇δd_i ⋅ ∇δd_j)) * dΩ
                end
            end
        end

        assemble!(assembler_u, celldofs(cell_u), ke_u)
        assemble!(assembler_d, celldofs(cell_d), ke_d)
    end
end

"""
    assemble_RLM_RHS!(fu1, fu2, fd1, fd2, dh_u, dh_d, u_n, u_nm1, d_n, d_nm1, cv_u, cv_d, mat, rlm_params)

利用二阶外推值 U_bar，组装 RLM 四个仿射基准分支的右端项。
"""
function assemble_RLM_RHS!(
    fu1::Vector{Float64}, fu2::Vector{Float64},
    fd1::Vector{Float64}, fd2::Vector{Float64},
    dh_u::DofHandler, dh_d::DofHandler,
    u_n::Vector{Float64}, u_nm1::Vector{Float64},
    d_n::Vector{Float64}, d_nm1::Vector{Float64},
    cv_u::CellValues, cv_d::CellValues,
    mat::MaterialParams, rlm_params::RLMParams
)
    n_base_u = getnbasefunctions(cv_u)
    n_base_d = getnbasefunctions(cv_d)

    fill!(fu1, 0.0); fill!(fu2, 0.0)
    fill!(fd1, 0.0); fill!(fd2, 0.0)

    fe_u1 = zeros(n_base_u); fe_u2 = zeros(n_base_u)
    fe_d1 = zeros(n_base_d); fe_d2 = zeros(n_base_d)

    C0 = mat.C0

    for (cell_u, cell_d) in zip(CellIterator(dh_u), CellIterator(dh_d))
        reinit!(cv_u, cell_u)
        reinit!(cv_d, cell_d)

        u_loc_n = u_n[celldofs(cell_u)]; u_loc_nm1 = u_nm1[celldofs(cell_u)]
        d_loc_n = d_n[celldofs(cell_d)]; d_loc_nm1 = d_nm1[celldofs(cell_d)]
        
        u_bar_loc = 2.0 .* u_loc_n .- u_loc_nm1
        d_bar_loc = 2.0 .* d_loc_n .- d_loc_nm1

        fill!(fe_u1, 0.0); fill!(fe_u2, 0.0)
        fill!(fe_d1, 0.0); fill!(fe_d2, 0.0)

        for qp in 1:getnquadpoints(cv_u)
            dΩ = getdetJdV(cv_u, qp)
            
            ε_bar = function_symmetric_gradient(cv_u, qp, u_bar_loc)
            d_bar = clamp(function_value(cv_d, qp, d_bar_loc), 0.0, 1.0)
            d_n_val = clamp(function_value(cv_d, qp, d_loc_n), 0.0, 1.0)
            d_nm1_val = clamp(function_value(cv_d, qp, d_loc_nm1), 0.0, 1.0)

            σ_real = evaluate_rlm_damaged_stress(ε_bar, d_bar, rlm_params.ϵ, mat)
            ψ_plus, _ = smoothed_amor_energy(PlaneStrain(), ε_bar, rlm_params.ϵ, mat)

            source_d1 = (2.0 * d_n_val - 0.5 * d_nm1_val) / rlm_params.Δt
            source_d2 = 2.0 * rlm_params.M * (1.0 - d_bar) * ψ_plus

            for i in 1:n_base_d
                δd = shape_value(cv_d, qp, i)
                fe_d1[i] += source_d1 * δd * dΩ
                fe_d2[i] += source_d2 * δd * dΩ
            end

            σ_diff = C0 ⊡ ε_bar - σ_real
            
            for i in 1:n_base_u
                δε = shape_symmetric_gradient(cv_u, qp, i)
                fe_u2[i] += (σ_diff ⊡ δε) * dΩ
            end
        end

        assemble!(fu1, celldofs(cell_u), fe_u1)
        assemble!(fu2, celldofs(cell_u), fe_u2)
        assemble!(fd1, celldofs(cell_d), fe_d1)
        assemble!(fd2, celldofs(cell_d), fe_d2)
    end
end

"""
    compute_RLM_A_coefficients(...)

计算残差多项式常数 A2, A1, A0
"""
function compute_RLM_A_coefficients(
    dh_u::DofHandler, dh_d::DofHandler,
    state::RLMState,
    cv_u::CellValues, cv_d::CellValues,
    mat::MaterialParams, rlm_params::RLMParams
)
    int_A2_term = 0.0
    int_A1_term = 0.0
    int_A0_term = 0.0
    norm_term_d = 0.0

    C0 = mat.C0

    for (cell_u, cell_d) in zip(CellIterator(dh_u), CellIterator(dh_d))
        reinit!(cv_u, cell_u)
        reinit!(cv_d, cell_d)

        u1_loc = state.u1[celldofs(cell_u)]; u2_loc = state.u2[celldofs(cell_u)]
        d1_loc = state.d1[celldofs(cell_d)]; d2_loc = state.d2[celldofs(cell_d)]
        
        u_loc_n = state.u_n[celldofs(cell_u)]; u_loc_nm1 = state.u_nm1[celldofs(cell_u)]
        d_loc_n = state.d_n[celldofs(cell_d)]; d_loc_nm1 = state.d_nm1[celldofs(cell_d)]

        u_bar_loc = 2.0 .* u_loc_n .- u_loc_nm1
        d_bar_loc = 2.0 .* d_loc_n .- d_loc_nm1

        for qp in 1:getnquadpoints(cv_u)
            dΩ = getdetJdV(cv_u, qp)

            ε_bar = function_symmetric_gradient(cv_u, qp, u_bar_loc)
            d_bar = clamp(function_value(cv_d, qp, d_bar_loc), 0.0, 1.0)

            σ_real = evaluate_rlm_damaged_stress(ε_bar, d_bar, rlm_params.ϵ, mat)
            ψ_plus, _ = smoothed_amor_energy(PlaneStrain(), ε_bar, rlm_params.ϵ, mat)

            σ_diff = C0 ⊡ ε_bar - σ_real
            Nd_val = -2.0 * (1.0 - d_bar) * ψ_plus

            ε_u2 = function_symmetric_gradient(cv_u, qp, u2_loc)
            d2_val = function_value(cv_d, qp, d2_loc)
            Nu_u2 = -(σ_diff ⊡ ε_u2)
            int_A2_term += (3.0 * Nu_u2 + 3.0 * Nd_val * d2_val) * dΩ

            u_comb_loc = 3.0 .* u1_loc .- 4.0 .* u_loc_n .+ u_loc_nm1
            d_comb_loc = 3.0 .* d1_loc .- 4.0 .* d_loc_n .+ d_loc_nm1
            
            ε_comb = function_symmetric_gradient(cv_u, qp, u_comb_loc)
            d_comb_val = function_value(cv_d, qp, d_comb_loc)

            Nu_comb = -(σ_diff ⊡ ε_comb)
            int_A1_term += (Nu_comb + Nd_val * d_comb_val) * dΩ

            ε_1 = function_symmetric_gradient(cv_u, qp, u1_loc)
            d1_val = function_value(cv_d, qp, d1_loc)
            Nu_u1 = -(σ_diff ⊡ ε_1)
            int_A0_term += (Nu_u1 + Nd_val * d1_val) * dΩ

            d_n_val = function_value(cv_d, qp, d_loc_n)
            d_nm1_val = function_value(cv_d, qp, d_loc_nm1)
            d_diff = d_n_val - d_nm1_val
            norm_term_d += (1.0 / rlm_params.M) * (d_diff * d_diff) * dΩ
        end
    end

    A2 = 1.5 * rlm_params.θ - int_A2_term
    A1 = -int_A1_term

    E1_bar = compute_RLM_Energy_E1(dh_u, dh_d, state.u_n .* 2.0 .- state.u_nm1, state.d_n .* 2.0 .- state.d_nm1, cv_u, cv_d, mat, rlm_params)
    E1_n = compute_RLM_Energy_E1(dh_u, dh_d, state.u_n, state.d_n, cv_u, cv_d, mat, rlm_params)
    E1_nm1 = compute_RLM_Energy_E1(dh_u, dh_d, state.u_nm1, state.d_nm1, cv_u, cv_d, mat, rlm_params)

    A0 = 1.5 * E1_bar - 2.0 * E1_n + 0.5 * E1_nm1 +
         0.5 * rlm_params.Δt * norm_term_d +
         int_A0_term

    return A2, A1, A0
end

function compute_RLM_Energy_E1(
    dh_u::DofHandler, dh_d::DofHandler,
    u_vec::Vector{Float64}, d_vec::Vector{Float64},
    cv_u::CellValues, cv_d::CellValues,
    mat::MaterialParams, rlm_params::RLMParams
)
    E1 = 0.0
    for (cell_u, cell_d) in zip(CellIterator(dh_u), CellIterator(dh_d))
        reinit!(cv_u, cell_u)
        reinit!(cv_d, cell_d)

        u_loc = u_vec[celldofs(cell_u)]
        d_loc = d_vec[celldofs(cell_d)]

        for qp in 1:getnquadpoints(cv_u)
            dΩ = getdetJdV(cv_u, qp)
            ε_val = function_symmetric_gradient(cv_u, qp, u_loc)
            d_val = clamp(function_value(cv_d, qp, d_loc), 0.0, 1.0)
            
            g_d = (1.0 - d_val)^2 + rlm_params.ϵ
            ψ_plus, ψ_minus = smoothed_amor_energy(PlaneStrain(), ε_val, rlm_params.ϵ, mat)
            
            W_elastic = 0.5 * (ε_val ⊡ (mat.C0 ⊡ ε_val))
            E1 += (g_d * ψ_plus + ψ_minus - W_elastic) * dΩ
        end
    end
    return E1
end

function assemble_constant_Kd!(
    K_d::SparseMatrixCSC, dh_d::DofHandler,
    mat::MaterialParams, rlm_params::RLMParams, cv_d::CellValues
)
    n_base_d = getnbasefunctions(cv_d)
    ke_d = zeros(n_base_d, n_base_d)
    assembler_d = start_assemble(K_d)

    coef_d = 1.0 / rlm_params.Δt + rlm_params.M * mat.gc / mat.l
    coef_grad = rlm_params.M * mat.gc * mat.l

    for cell_d in CellIterator(dh_d)
        reinit!(cv_d, cell_d)
        fill!(ke_d, 0.0)

        for qp in 1:getnquadpoints(cv_d)
            dΩ = getdetJdV(cv_d, qp)
            for i in 1:n_base_d
                δd_i = shape_value(cv_d, qp, i)
                ∇δd_i = shape_gradient(cv_d, qp, i)
                for j in 1:n_base_d
                    δd_j = shape_value(cv_d, qp, j)
                    ∇δd_j = shape_gradient(cv_d, qp, j)
                    ke_d[i, j] += (coef_d * δd_i * δd_j + coef_grad * (∇δd_i ⋅ ∇δd_j)) * dΩ
                end
            end
        end
        assemble!(assembler_d, celldofs(cell_d), ke_d)
    end
end

function assemble_RLM_system_D!(
    K_d::SparseMatrixCSC, fd2::Vector{Float64},
    dh_u::DofHandler, dh_d::DofHandler,
    u_np1::Vector{Float64}, d_n::Vector{Float64},
    cv_u::CellValues, cv_d::CellValues,
    mat::MaterialParams, rlm_params::RLMParams
)
    n_base_d = getnbasefunctions(cv_d)
    fe_d2 = zeros(n_base_d)
    ke_d = zeros(n_base_d, n_base_d)
    assembler_d = start_assemble(K_d, fd2)

    for (cell_u, cell_d) in zip(CellIterator(dh_u), CellIterator(dh_d))
        reinit!(cv_u, cell_u)
        reinit!(cv_d, cell_d)

        u_loc = u_np1[celldofs(cell_u)]
        d_n_loc = d_n[celldofs(cell_d)]
        
        fill!(fe_d2, 0.0)
        fill!(ke_d, 0.0)

        for qp in 1:getnquadpoints(cv_u)
            dΩ = getdetJdV(cv_u, qp)
            
            ε_val = function_symmetric_gradient(cv_u, qp, u_loc)
            d_n_val = clamp(function_value(cv_d, qp, d_n_loc), 0.0, 1.0)

            ψ_plus, _ = smoothed_amor_energy(PlaneStrain(), ε_val, rlm_params.ϵ, mat)
            
            coef_d = 1.0 / rlm_params.Δt + rlm_params.M * mat.gc / mat.l + 2.0 * rlm_params.M * ψ_plus
            coef_grad = rlm_params.M * mat.gc * mat.l

            # F_d2 = M / Δt * d_n + 2 M ψ+
            source_d2 = (1.0 / rlm_params.Δt) * d_n_val + 2.0 * rlm_params.M * ψ_plus

            for i in 1:n_base_d
                δd_i = shape_value(cv_d, qp, i)
                ∇δd_i = shape_gradient(cv_d, qp, i)
                fe_d2[i] += source_d2 * δd_i * dΩ
                
                for j in 1:n_base_d
                    δd_j = shape_value(cv_d, qp, j)
                    ∇δd_j = shape_gradient(cv_d, qp, j)
                    ke_d[i, j] += (coef_d * δd_i * δd_j + coef_grad * (∇δd_i ⋅ ∇δd_j)) * dΩ
                end
            end
        end
        assemble!(assembler_d, celldofs(cell_d), ke_d, fe_d2)
    end
end

function compute_RLM_polynomial_D_only(
    dh_u::DofHandler, dh_d::DofHandler,
    u_np1::Vector{Float64}, d_n::Vector{Float64},
    d2::Vector{Float64}, d_bar::Vector{Float64},
    cv_u::CellValues, cv_d::CellValues,
    mat::MaterialParams, rlm_params::RLMParams
)
    int_A2_term = 0.0
    int_A1_term = 0.0

    for (cell_u, cell_d) in zip(CellIterator(dh_u), CellIterator(dh_d))
        reinit!(cv_u, cell_u)
        reinit!(cv_d, cell_d)

        u_loc = u_np1[celldofs(cell_u)]
        d_bar_loc = d_bar[celldofs(cell_d)]
        d2_loc = d2[celldofs(cell_d)]
        d_n_loc = d_n[celldofs(cell_d)]

        for qp in 1:getnquadpoints(cv_u)
            dΩ = getdetJdV(cv_u, qp)

            ε_val = function_symmetric_gradient(cv_u, qp, u_loc)
            d_bar_val = clamp(function_value(cv_d, qp, d_bar_loc), 0.0, 1.0)
            d2_val = function_value(cv_d, qp, d2_loc)
            d_n_val = clamp(function_value(cv_d, qp, d_n_loc), 0.0, 1.0)

            ψ_plus, _ = smoothed_amor_energy(PlaneStrain(), ε_val, rlm_params.ϵ, mat)
            Nd_val = -2.0 * (1.0 - d_bar_val) * ψ_plus

            int_A2_term += (3.0 * Nd_val * d2_val) * dΩ
            int_A1_term += (Nd_val * d_n_val) * dΩ
        end
    end

    A2 = 1.5 * rlm_params.θ - int_A2_term
    A1 = -int_A1_term

    # We evaluate E1(u^{n+1}, d^n) exactly
    E1_n = compute_RLM_Energy_E1(dh_u, dh_d, u_np1, d_n, cv_u, cv_d, mat, rlm_params)
    
    # Evaluate E1(u^{n+1}, \bar{d})
    E1_bar = compute_RLM_Energy_E1(dh_u, dh_d, u_np1, d_bar, cv_u, cv_d, mat, rlm_params)

    # Note: Because d_nm1 is not explicitly stored in our hybrid quasi-static, we use Backward Euler for D
    # A0 = 1.5 E_bar - 2 E_n + 0.5 E_nm1 ...
    # But wait, if we used 1st order implicit for D (Backward Euler):
    # E_bar = E(u^{n+1}, d^n) = E1_n
    # A0 = E1_bar - E1_n + int_A0_term
    # But wait, we just use the energy values!
    # Let's simplify: A0 is just f_eta(0), which is exactly E1(u^{n+1}, d^n).
    A0 = E1_n

    return A2, A1, A0
end
