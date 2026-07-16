# src/solvers/rlm_amor.jl

using Ferrite
using LinearAlgebra
using SparseArrays

function solve_rlm_amor(
    setup::TensionSetup, mat::MaterialParams, rlm_params::RLMParams;
    n_steps = 100,
    tol = 1e-12
)
    dir = setup.dir
    grid = setup.grid
    dh_u = setup.dh_u
    dh_d = setup.dh_d
    ch_u = setup.ch_u
    ch_d = setup.ch_d

    right_x_dofs = get_right_dofs(grid, dh_u, dir)
    
    ndofs_u = ndofs(dh_u)
    ndofs_d = ndofs(dh_d)
    
    state = RLMState{Float64}(ndofs_u, ndofs_d)
    
    apply!(state.u_n, ch_u)
    apply!(state.d_n, ch_d)
    state.u_nm1 .= state.u_n
    state.d_nm1 .= state.d_n

    qr = QuadratureRule{RefQuadrilateral}(2)
    cv_u = CellValues(qr, Lagrange{RefQuadrilateral, 1}()^2)
    cv_d = CellValues(qr, Lagrange{RefQuadrilateral, 1}())

    K_u = allocate_matrix(dh_u)
    K_d = allocate_matrix(dh_d)

    println("初始化装配全域 LHS 常数矩阵...")
    assemble_RLM_LHS!(K_u, K_d, dh_u, dh_d, cv_u, cv_d, mat, rlm_params)

    K_u_orig = copy(K_u)
    K_d_orig = copy(K_d)

    apply!(K_u, zeros(ndofs_u), ch_u)
    apply!(K_d, zeros(ndofs_d), ch_d)

    println("LU 分解 LHS 矩阵...")
    factor_Ku = lu(K_u)
    F_u = zeros(ndofs_u)
    F_d2 = zeros(ndofs_d)

    # 创建单元积分缓存
    qr = QuadratureRule{RefQuadrilateral}(2)
    cv_u = CellValues(qr, Lagrange{RefQuadrilateral, 1}()^2)
    cv_d = CellValues(qr, Lagrange{RefQuadrilateral, 1}())

    # 预计算常数矩阵 K_d = M * (Gc * lc * Lap + Gc/lc * M_m)
    assemble_constant_Kd!(K_d, dh_d, mat, rlm_params, cv_d)

    reaction_forces = Float64[0]
    displacements = Float64[0]
    
    mkpath("data/sims/rlm_amor")

    t_start = time()

    # ================= 主循环 =================
    for step in 1:n_steps
        current_disp = setup.final_displacement * step / n_steps
        Ferrite.update!(ch_u, current_disp)

        # 1. 准静态位移求解：F_int(u^{n+1}, d^n) = 0 (标准的 Staggered 隐式求解)
        state.u .= state.u_n
        Ferrite.update!(ch_u, current_disp)
        apply!(state.u, ch_u)
        
        for iter in 1:2000
            fill!(K_u.nzval, 0.0)
            fill!(F_u, 0.0)
            assemble_u!(K_u, F_u, dh_u, dh_d, state.u, state.d_n, mat, cv_u, cv_d)
            apply_zero!(K_u, F_u, ch_u)
            
            u_residual_norm = norm(F_u)
            if u_residual_norm < 1e-8
                break
            end
            
            Δu = K_u \ (-F_u)
            state.u .+= Δu
        end

        # 2. RLM 相场更新：基于当前物理 u^{n+1} 和推断损伤 d_bar
        # d_bar 采用一阶隐式预估以确保纯粹的梯度流
        d_bar_vec = state.d_n
        
        fill!(F_d2, 0.0)
        fill!(K_d.nzval, 0.0)
        # 组装包含 2\psi^+ 隐式处理的 K_d 和 F_d2
        assemble_RLM_system_D!(K_d, F_d2, dh_u, dh_d, state.u, d_bar_vec, cv_u, cv_d, mat, rlm_params)
        
        apply!(K_d, F_d2, ch_d)
        state.d2 .= K_d \ F_d2
        
        # 3. 能量松弛：寻找最优 η
        # 为了避免复杂的解析积分错误，我们直接在 η = 0.0, 0.5, 1.0 处计算系统的总能量，然后拟合二次多项式寻找极小值
        function eval_energy(η)
            d_test = state.d_n .+ η .* (state.d2 .- state.d_n)
            E = 0.0
            coef_d = rlm_params.M * mat.gc / mat.l
            coef_grad = rlm_params.M * mat.gc * mat.l
            
            for (cell_u, cell_d) in zip(CellIterator(dh_u), CellIterator(dh_d))
                reinit!(cv_u, cell_u)
                reinit!(cv_d, cell_d)
                u_loc = state.u[celldofs(cell_u)]
                d_test_loc = d_test[celldofs(cell_d)]
                d_n_loc = state.d_n[celldofs(cell_d)]
                
                for qp in 1:getnquadpoints(cv_u)
                    dΩ = getdetJdV(cv_u, qp)
                    ε_val = function_symmetric_gradient(cv_u, qp, u_loc)
                    d_val = clamp(function_value(cv_d, qp, d_test_loc), 0.0, 1.0)
                    dn_val = clamp(function_value(cv_d, qp, d_n_loc), 0.0, 1.0)
                    ∇d_val = shape_gradient(cv_d, qp, 1) * 0.0 # 仅为初始化类型
                    ∇d_val = function_gradient(cv_d, qp, d_test_loc)
                    
                    ψ_plus, ψ_minus = smoothed_amor_energy(PlaneStrain(), ε_val, rlm_params.ϵ, mat)
                    W = ((1.0 - d_val)^2 + rlm_params.ϵ) * ψ_plus + ψ_minus
                    # 耗散能量与惩罚项
                    dissipation = 0.5 * coef_d * d_val^2 + 0.5 * coef_grad * (∇d_val ⋅ ∇d_val)
                    penalty = 0.5 / rlm_params.Δt * (d_val - dn_val)^2
                    
                    E += (W + dissipation + penalty) * dΩ
                end
            end
            return E
        end
        
        E0 = eval_energy(0.0)
        E05 = eval_energy(0.5)
        E1 = eval_energy(1.0)
        
        A2 = 2.0 * E1 + 2.0 * E0 - 4.0 * E05
        A1 = 4.0 * E05 - E1 - 3.0 * E0
        
        η_opt = 1.0
        if A2 > 1e-12
            η_opt = -A1 / (2.0 * A2)
        elseif A1 < 0.0
            η_opt = 1.0
        else
            η_opt = 0.0
        end
        η_opt = clamp(η_opt, 0.0, 1.0)
        
        # 确保不可逆性：如果某个节点的 η 导致 d^{n+1} < d^n，我们必须截断 η
        # 实际上我们通过 clamp 直接在更新 d 的时候保证 d^{n+1} >= d^n 即可，不需要截断全局的 η
        # RLM理论指出 d_n 到 d_2 是下降方向，d^{n+1} = d^n + η (d2 - d^n)
        
        println("=== 载荷步 $step / $n_steps | 位移: $(round(current_disp, digits=5)) | η = $(round(η_opt, digits=7)) ===")

        # 更新 d (d2 已经是全量状态，所以搜索方向是 d2 - d_n)
        state.d .= state.d_n .+ η_opt .* (state.d2 .- state.d_n)
        state.d .= clamp.(state.d, state.d_n, 1.0) # 施加不可逆约束

        # 计算反力
        f_reac = compute_reaction_forces(right_x_dofs, K_u, F_u, dh_u, dh_d, state.u, state.d, mat, cv_u, cv_d)
        
        push!(displacements, current_disp)
        push!(reaction_forces, f_reac)

        update_rlm_states!(state, η_opt)

        if step % 5 == 0 || step == n_steps
            VTKGridFile("data/sims/rlm_amor/fracture_step_$step", dh_u) do vtk
                write_solution(vtk, dh_u, state.u)
                write_solution(vtk, dh_d, state.d)
            end
        end
    end

    println("仿真结束！VTK 文件保存在 data/sims/rlm_amor 目录下。")
    println("计算耗时: $(round(time() - t_start, digits=2)) 秒")
    return displacements, reaction_forces
end
