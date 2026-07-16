# src/solvers/sav.jl

using Ferrite
using LinearAlgebra
using SparseArrays

"""
    solve_sav(setup, mat, nat; T_final=0.01, enforce_irreversibility=true)

执行 SAV 完全线性化交替求解。
"""
function solve_sav(
    setup::TensionSetup, mat::MaterialParams, nat::NumericalParams;
    T_final = 0.01,
    enforce_irreversibility::Bool = true
)
    # --- 1. 提取网格与自由度 ---
    grid  = setup.grid
    dh_u  = setup.dh_u
    dh_d  = setup.dh_d
    ch_u = setup.ch_u
    ch_d  = setup.ch_d

    ndofs_u = ndofs(dh_u)
    ndofs_d = ndofs(dh_d)

    right_x_dofs = get_right_dofs(grid, dh_u, setup.dir)

    # --- 2. 准备积分法则和 CellValues ---
    qr   = QuadratureRule{RefQuadrilateral}(2)
    cv_u = CellValues(qr, Lagrange{RefQuadrilateral,1}()^2)
    cv_d = CellValues(qr, Lagrange{RefQuadrilateral,1}())

    # --- 3. 初始化仿真状态 ---
    state = SimulationState{Float64}(ndofs_u, ndofs_d, 0.0)
    # 计算 t=0 时的初始非线性能量
    E_nl_0 = nonlinear_energy(dh_u, dh_d, state.u_n, state.d_n, mat, cv_u, cv_d)
    r_0 = sqrt(E_nl_0 + nat.S0)
    state = SimulationState{Float64}(ndofs_u, ndofs_d, r_0)
    state.r = r_0

    apply!(state.u_n, ch_u)
    apply!(state.d_n, ch_d)

    n_qpoints = getncells(grid) * getnquadpoints(cv_u)
    driving_force = zeros(n_qpoints)

    # --- 4. 分配全局稀疏矩阵 ---
    K_u = allocate_matrix(dh_u, ch_u)
    K_d = allocate_matrix(dh_d, ch_d)
    K_u_nl = allocate_matrix(dh_u)   # 用于反力计算（不与 SAV 矩阵共享）
    R_u    = zeros(ndofs_u)

    # RHS 工作向量
    fd1 = zeros(ndofs_d)
    fd2 = zeros(ndofs_d)
    fu1 = zeros(ndofs_u)
    fu2 = zeros(ndofs_u)

    # --- 5. 记录用数组 ---
    reaction_forces = Float64[0.0]
    displacements   = Float64[0.0]
    elastic_energies = Float64[0.0]
    surface_energies = Float64[0.0]

    # VTK 输出目录
    outdir = enforce_irreversibility ? "data/sims/sav" : "data/sims/sav2"
    mkpath(outdir)

    t_start = time()

    # --- 6. 组装并因子分解恒定的 SAV 左端矩阵 ---
    assemble_SAV_LHS!(K_d, K_u, dh_u, dh_d, cv_u, cv_d, mat, nat)

    rhsdata_d  = get_rhs_data(ch_d, K_d)
    rhsdata_u1 = get_rhs_data(ch_u, K_u)

    apply!(K_d, ch_d)
    apply!(K_u, ch_u)

    # 因子分解
    F_Kd = cholesky(Symmetric(K_d))
    F_Ku = cholesky(Symmetric(K_u))

    n_steps = Int(round(T_final / nat.Δt))
    println("矩阵分解完成 ($(typeof(F_Kd).name.name), $(typeof(F_Ku).name.name))，开始 SAV 格式求解，总步数: $n_steps")

    # --- 7. 时间步循环 ---
    for step in 1:n_steps
        u_n = state.u_n
        d_n = state.d_n
        v_n = state.v_n
        r_n = state.r_n

        current_t = step * nat.Δt
        current_disp = current_t * setup.final_displacement
        update!(ch_u, current_t)
        println("=== 载荷步 $step / $n_steps | 位移: $(round(current_disp, digits=5)) ===")

        # --- 计算非线性能量和驱动力 ---
        E_nl_n = nonlinear_energy(dh_u, dh_d, u_n, d_n, mat, cv_u, cv_d)
        compute_driving_force!(driving_force, dh_u, u_n, mat, cv_u, enforce_irreversibility)
        ξ_n = driving_force ./ sqrt(E_nl_n + nat.S0)

        # --- 组装未处理边界的原始 RHS ---
        fill!(fd1, 0.0); fill!(fd2, 0.0)
        fill!(fu1, 0.0); fill!(fu2, 0.0)
        assemble_SAV_RHS!(fd1, fd2, fu1, fu2, dh_u, dh_d,
            u_n, v_n, d_n, ξ_n, E_nl_n, cv_u, cv_d, mat, nat)

        # --- 施加边界条件到 RHS （核心修正：利用 RHSData 进行列消元） ---
        apply_rhs!(rhsdata_d, fd1, ch_d, false)  # d1 施加真实边界条件
        apply_rhs!(rhsdata_d, fd2, ch_d, true)   # d2 施加齐次（零）边界条件

        apply_rhs!(rhsdata_u1, fu1, ch_u, false) # u1 施加真实位移加载
        apply_rhs!(rhsdata_u1, fu2, ch_u, true)  # u2 施加齐次（零）位移约束

        # --- 步骤 1: 回代求解 d1, d2 并后处理约束 ---
        d1 = F_Kd \ fd1
        d2 = F_Kd \ fd2
        apply!(d1, ch_d)
        apply_zero!(d2, ch_d)

        # --- 步骤 2: 回代求解 u1, u2 并后处理约束 ---
        u1 = F_Ku \ fu1
        u2 = F_Ku \ fu2
        apply!(u1, ch_u)
        apply_zero!(u2, ch_u)

        # --- 步骤 3: 计算全局标量 r^{n+1} = B / A ---
        A, B = compute_sav_scalars(dh_u, dh_d, u_n, d_n,
            u1, u2, d1, d2, r_n, E_nl_n, nat, mat, cv_u, cv_d)
        r_np1 = B / A

        # --- 步骤 4: 线性叠加重建物理场 ---
        d_np1 = d1 + r_np1 .* d2
        u_np1 = u1 + r_np1 .* u2
        v_np1 = (u_np1 - u_n) / nat.Δt

        # 确保重建后的物理场严格满足边界
        apply!(u_np1, ch_u)
        apply!(d_np1, ch_d)

        # --- 存入 state 并滚动历史 ---
        copyto!(state.u, u_np1)
        copyto!(state.d, d_np1)
        copyto!(state.v, v_np1)
        state.r = r_np1
        update_states!(state)

        # --- 计算反力 ---
        fill!(R_u, 0.0)
        assemble_u!(K_u_nl, R_u, dh_u, dh_d, u_np1, d_np1, mat, cv_u, cv_d)
        f_total = sum(R_u[dof] for dof in right_x_dofs)

        push!(displacements, current_disp)
        push!(reaction_forces, f_total)

        # --- 能量 ---
        push!(elastic_energies, elastic_energy(dh_u, dh_d, u_np1, d_np1, mat, cv_u, cv_d))
        push!(surface_energies, surface_energy(dh_d, d_np1, mat, cv_d))

        # --- VTK 输出 ---
        if step % 5 == 0 || step == n_steps
            VTKGridFile(joinpath(outdir, "fracture_step_$step"), dh_u) do vtk
                write_solution(vtk, dh_u, u_np1)
                write_solution(vtk, dh_d, d_np1)
            end
        end
    end

    println("仿真结束！VTK 文件保存在 $outdir 目录下。")
    println("计算耗时: $(round(time() - t_start, digits=2)) 秒")
    return displacements, reaction_forces, elastic_energies, surface_energies
end
