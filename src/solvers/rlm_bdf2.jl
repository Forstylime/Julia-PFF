# src/solvers/rlm_bdf2.jl

using Ferrite
using LinearAlgebra
using SparseArrays

"""
一维求根函数：Brent 方法 (无导数)
"""
function brent_root(f, x0::Float64, x1::Float64; tol::Float64=1e-9, max_iter::Int=100)
    a = x0
    b = x1
    fa = f(a)
    fb = f(b)

    if fa * fb > 0
        # 如果初始区间没有变号，为了鲁棒性，强制拓宽区间寻找变号点
        # (在物理问题中，eta 约等于 1.0)
        for i in 1:10
            a -= 0.5
            b += 0.5
            fa = f(a)
            fb = f(b)
            if fa * fb <= 0
                break
            end
        end
        if fa * fb > 0
            # 仍然没找到，退化返回 1.0
            return 1.0
        end
    end

    if abs(fa) < abs(fb)
        a, b = b, a
        fa, fb = fb, fa
    end

    c = a
    fc = fa
    s = 0.0
    fs = 0.0
    d = 0.0
    mflag = true

    for iter in 1:max_iter
        if fa != fc && fb != fc
            # 逆二次插值
            s = a * fb * fc / ((fa - fb) * (fa - fc)) +
                b * fa * fc / ((fb - fa) * (fb - fc)) +
                c * fa * fb / ((fc - fa) * (fc - fb))
        else
            # 割线法
            s = b - fb * (b - a) / (fb - fa)
        end

        # 判断是否需要二分
        condition1 = !((s > (3 * a + b) / 4) && (s < b)) && !((s < (3 * a + b) / 4) && (s > b))
        condition2 = mflag && abs(s - b) >= abs(b - c) / 2
        condition3 = !mflag && abs(s - b) >= abs(c - d) / 2
        condition4 = mflag && abs(b - c) < tol
        condition5 = !mflag && abs(c - d) < tol

        if condition1 || condition2 || condition3 || condition4 || condition5
            s = (a + b) / 2
            mflag = true
        else
            mflag = false
        end

        fs = f(s)
        d = c
        c = b
        fc = fb

        if fa * fs < 0
            b = s
            fb = fs
        else
            a = s
            fa = fs
        end

        if abs(fa) < abs(fb)
            a, b = b, a
            fa, fb = fb, fa
        end

        if abs(b - a) < tol || fb == 0.0
            return b
        end
    end
    return b
end

function compute_E1(dh_u::DofHandler, dh_d::DofHandler, u_vec::Vector{Float64}, d_vec::Vector{Float64}, mat::MaterialParams, rlm_params::RLMParams, cv_u::CellValues, cv_d::CellValues)
    E1 = 0.0
    for (cell_u, cell_d) in zip(CellIterator(dh_u), CellIterator(dh_d))
        reinit!(cv_u, cell_u)
        reinit!(cv_d, cell_d)
        
        u_loc = u_vec[celldofs(cell_u)]
        d_loc = d_vec[celldofs(cell_d)]
        
        for qp in 1:getnquadpoints(cv_u)
            dΩ = getdetJdV(cv_u, qp)
            
            ε_val = function_symmetric_gradient(cv_u, qp, u_loc)
            d_val = function_value(cv_d, qp, d_loc)
            # clamp damage to [0,1] for energy evaluation to prevent unphysical values during line search
            d_val = clamp(d_val, 0.0, 1.0)
            
            g_d = (1.0 - d_val)^2 + rlm_params.ϵ
            ψ_plus, ψ_minus = smoothed_amor_energy(PlaneStrain(), ε_val, rlm_params.ϵ, mat)
            
            W_elastic = 0.5 * (ε_val ⊡ (mat.C0 ⊡ ε_val))
            E1 += (g_d * ψ_plus + ψ_minus - W_elastic) * dΩ
        end
    end
    return E1
end

function assemble_RLM_BDF2_LHS!(
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

    # K_u = L_u = C0
    # K_d = (3 / (2 * M * Δt) + Gc / l) * M_d + Gc * l * Lap_d
    coef_d_mass = 1.5 / (rlm_params.M * rlm_params.Δt) + mat.gc / mat.l
    coef_d_lap  = mat.gc * mat.l

    for (cell_u, cell_d) in zip(CellIterator(dh_u), CellIterator(dh_d))
        reinit!(cv_u, cell_u)
        reinit!(cv_d, cell_d)

        fill!(ke_u, 0.0)
        fill!(ke_d, 0.0)

        for qp in 1:getnquadpoints(cv_u)
            dΩ = getdetJdV(cv_u, qp)

            # u LHS
            for i in 1:n_base_u
                δε_i = shape_symmetric_gradient(cv_u, qp, i)
                for j in 1:n_base_u
                    δε_j = shape_symmetric_gradient(cv_u, qp, j)
                    ke_u[i, j] += (δε_i ⊡ mat.C0 ⊡ δε_j) * dΩ
                end
            end

            # d LHS
            for i in 1:n_base_d
                δd_i = shape_value(cv_d, qp, i)
                ∇δd_i = shape_gradient(cv_d, qp, i)
                for j in 1:n_base_d
                    δd_j = shape_value(cv_d, qp, j)
                    ∇δd_j = shape_gradient(cv_d, qp, j)
                    ke_d[i, j] += (coef_d_mass * δd_i * δd_j + coef_d_lap * (∇δd_i ⋅ ∇δd_j)) * dΩ
                end
            end
        end

        assemble!(assembler_u, celldofs(cell_u), ke_u)
        assemble!(assembler_d, celldofs(cell_d), ke_d)
    end
end

function assemble_RLM_BDF2_RHS!(
    F_u2::Vector{Float64}, F_d1::Vector{Float64}, F_d2::Vector{Float64},
    dh_u::DofHandler, dh_d::DofHandler,
    u_bar::Vector{Float64}, d_bar::Vector{Float64},
    d_n::Vector{Float64}, d_nm1::Vector{Float64},
    cv_u::CellValues, cv_d::CellValues,
    mat::MaterialParams, rlm_params::RLMParams
)
    n_base_u = getnbasefunctions(cv_u)
    n_base_d = getnbasefunctions(cv_d)
    
    fe_u2 = zeros(n_base_u)
    fe_d1 = zeros(n_base_d)
    fe_d2 = zeros(n_base_d)

    for (cell_u, cell_d) in zip(CellIterator(dh_u), CellIterator(dh_d))
        reinit!(cv_u, cell_u)
        reinit!(cv_d, cell_d)

        u_bar_loc = u_bar[celldofs(cell_u)]
        d_bar_loc = d_bar[celldofs(cell_d)]
        d_n_loc   = d_n[celldofs(cell_d)]
        d_nm1_loc = d_nm1[celldofs(cell_d)]

        fill!(fe_u2, 0.0)
        fill!(fe_d1, 0.0)
        fill!(fe_d2, 0.0)

        for qp in 1:getnquadpoints(cv_u)
            dΩ = getdetJdV(cv_u, qp)
            
            ε_bar = function_symmetric_gradient(cv_u, qp, u_bar_loc)
            d_bar_val = clamp(function_value(cv_d, qp, d_bar_loc), 0.0, 1.0)
            d_n_val   = clamp(function_value(cv_d, qp, d_n_loc), 0.0, 1.0)
            d_nm1_val = clamp(function_value(cv_d, qp, d_nm1_loc), 0.0, 1.0)

            # --- u2 RHS (N_u) ---
            g_d = (1.0 - d_bar_val)^2 + rlm_params.ϵ
            σ_plus, σ_minus = smoothed_amor_stress(PlaneStrain(), ε_bar, rlm_params.ϵ, mat)
            σ_real = g_d * σ_plus + σ_minus
            σ_linear = mat.C0 ⊡ ε_bar

            # --- d2 RHS (N_d) ---
            ψ_plus, _ = smoothed_amor_energy(PlaneStrain(), ε_bar, rlm_params.ϵ, mat)
            # N_d = -2(1-d_bar)*ψ_plus,  F_d2 = - M * N_d = 2 M (1-d_bar) ψ_plus
            source_d2 = 2.0 * rlm_params.M * (1.0 - d_bar_val) * ψ_plus

            # --- d1 RHS ---
            source_d1 = (4.0 * d_n_val - d_nm1_val) / (2.0 * rlm_params.M * rlm_params.Δt)

            for i in 1:n_base_u
                δε_i = shape_symmetric_gradient(cv_u, qp, i)
                fe_u2[i] += (δε_i ⊡ (σ_linear - σ_real)) * dΩ
            end

            for i in 1:n_base_d
                δd_i = shape_value(cv_d, qp, i)
                fe_d1[i] += source_d1 * δd_i * dΩ
                fe_d2[i] += source_d2 * δd_i * dΩ
            end
        end

        assemble!(F_u2, celldofs(cell_u), fe_u2)
        assemble!(F_d1, celldofs(cell_d), fe_d1)
        assemble!(F_d2, celldofs(cell_d), fe_d2)
    end
end

function solve_rlm_bdf2(setup::TensionSetup, mat::MaterialParams, rlm_params::RLMParams; n_steps=500)
    grid = setup.grid
    dh_u = setup.dh_u
    dh_d = setup.dh_d
    ch_u = setup.ch_u
    ch_d = setup.ch_d

    ndofs_u = ndofs(dh_u)
    ndofs_d = ndofs(dh_d)
    right_x_dofs = get_right_dofs(grid, dh_u, setup.dir)

    qr = QuadratureRule{RefQuadrilateral}(2)
    cv_u = CellValues(qr, Lagrange{RefQuadrilateral, 1}()^2)
    cv_d = CellValues(qr, Lagrange{RefQuadrilateral, 1}())

    K_u = allocate_matrix(dh_u)
    K_d = allocate_matrix(dh_d)

    println("初始化装配全域 LHS 常数矩阵...")
    assemble_RLM_BDF2_LHS!(K_u, K_d, dh_u, dh_d, cv_u, cv_d, mat, rlm_params)

    # 抽取用于 apply_rhs! 的信息
    rhsdata_u = get_rhs_data(ch_u, K_u)
    rhsdata_d = get_rhs_data(ch_d, K_d)

    apply!(K_u, ch_u)
    apply!(K_d, ch_d)

    println("LU 分解 LHS 矩阵...")
    factor_Ku = cholesky(Symmetric(K_u))
    factor_Kd = cholesky(Symmetric(K_d))

    u_n = zeros(ndofs_u)
    u_nm1 = zeros(ndofs_u)
    d_n = zeros(ndofs_d)
    d_nm1 = zeros(ndofs_d)
    
    eta_n = 1.0
    eta_nm1 = 1.0

    F_u1 = zeros(ndofs_u)
    F_u2 = zeros(ndofs_u)
    F_d1 = zeros(ndofs_d)
    F_d2 = zeros(ndofs_d)

    reaction_forces = Float64[0.0]
    displacements = Float64[0.0]

    mkpath("data/sims/rlm_bdf2")

    t_start = time()

    for step in 1:n_steps
        current_disp = setup.final_displacement * step / n_steps
        update!(ch_u, current_disp)

        u_bar = 2 .* u_n .- u_nm1
        d_bar = 2 .* d_n .- d_nm1

        fill!(F_u2, 0.0)
        fill!(F_d1, 0.0)
        fill!(F_d2, 0.0)

        assemble_RLM_BDF2_RHS!(F_u2, F_d1, F_d2, dh_u, dh_d, u_bar, d_bar, d_n, d_nm1, cv_u, cv_d, mat, rlm_params)

        F_u2_raw = copy(F_u2)
        F_d2_raw = copy(F_d2)

        fill!(F_u1, 0.0)

        apply_rhs!(rhsdata_u, F_u1, ch_u, false)
        apply_rhs!(rhsdata_u, F_u2, ch_u, true)
        apply_rhs!(rhsdata_d, F_d1, ch_d, false)
        apply_rhs!(rhsdata_d, F_d2, ch_d, true)

        u1 = factor_Ku \ F_u1
        u2 = factor_Ku \ F_u2
        d1 = factor_Kd \ F_d1
        d2 = factor_Kd \ F_d2

        apply!(u1, ch_u)
        apply_zero!(u2, ch_u)
        apply!(d1, ch_d)
        apply_zero!(d2, ch_d)

        A2 = 1.5 * rlm_params.θ + 3 * dot(u2, F_u2_raw) + (3.0 / rlm_params.M) * dot(d2, F_d2_raw)
        A1 = dot(3 .* u1 .- 4 .* u_n .+ u_nm1, F_u2_raw) + (1.0 / rlm_params.M) * dot(3 .* d1 .- 4 .* d_n .+ d_nm1, F_d2_raw)
        
        
        E1_n = compute_E1(dh_u, dh_d, u_n, d_n, mat, rlm_params, cv_u, cv_d)
        E1_nm1 = compute_E1(dh_u, dh_d, u_nm1, d_nm1, mat, rlm_params, cv_u, cv_d)
        A0 = -2 * E1_n + 0.5 * E1_nm1 - rlm_params.θ * (2 * eta_n^2 - 0.5 * eta_nm1^2)

        function f_eta(η)
            u_eval = u1 .+ η .* u2
            d_eval = d1 .+ η .* d2
            return 1.5 * compute_E1(dh_u, dh_d, u_eval, d_eval, mat, rlm_params, cv_u, cv_d) + A2 * η^2 + A1 * η + A0
        end

        eta_opt = brent_root(f_eta, 0.8, 1.2)
        
        u_np1 = u1 .+ eta_opt .* u2
        d_np1 = d1 .+ eta_opt .* d2
        
        # 强制单调递增损伤
        d_np1 .= max.(d_np1, d_n)
        d_np1 .= clamp.(d_np1, 0.0, 1.0)
        
        apply!(u_np1, ch_u)
        apply!(d_np1, ch_d)

        # 更新历史
        u_nm1 .= u_n
        u_n .= u_np1
        d_nm1 .= d_n
        d_n .= d_np1
        eta_nm1 = eta_n
        eta_n = eta_opt

        # 计算反力
        K_u_nl = allocate_matrix(dh_u)
        R_u = zeros(ndofs_u)
        assemble_u!(K_u_nl, R_u, dh_u, dh_d, u_np1, d_np1, mat, cv_u, cv_d)
        f_total = sum(R_u[dof] for dof in right_x_dofs)
        
        push!(displacements, current_disp)
        push!(reaction_forces, f_total)

        println("=== 载荷步 $step / $n_steps | 位移: $(round(current_disp, digits=4)) | η = $(round(eta_opt, digits=7)) ===")

        if step % 5 == 0 || step == n_steps
            VTKGridFile(joinpath("data/sims/rlm_bdf2", "fracture_step_$step"), dh_u) do vtk
                write_solution(vtk, dh_u, u_np1)
                write_solution(vtk, dh_d, d_np1)
            end
        end
    end

    println("仿真结束！VTK 文件保存在 data/sims/rlm_bdf2 目录下。")
    println("计算耗时: $(round(time() - t_start, digits=2)) 秒")
    return displacements, reaction_forces
end
