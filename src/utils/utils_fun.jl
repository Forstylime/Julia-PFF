"""
计算相场演化的驱动力，即历史场变量 H 或者 ψ+。
如果 enforce_irreversibility = true，则取历史最大值 (L型试件)；
如果 enforce_irreversibility = false，则直接使用当前应变能 (CT型试件)。
"""
function compute_driving_force!(
    driving_force::Vector{Float64}, 
    dh_u::DofHandler, u_global::Vector{Float64}, 
    mat::MaterialParams, cv_u::CellValues,
    enforce_irreversibility::Bool
)
    qp_count = 1
    for cell in CellIterator(dh_u)
        reinit!(cv_u, cell)
        u_loc = u_global[celldofs(cell)]
        
        for qp in 1:getnquadpoints(cv_u)
            # 计算当前积分点的应变
            ε_q = function_symmetric_gradient(cv_u, qp, u_loc)
            
            # 调用你写的拉伸应变能密度函数
            Ψ_plus = tensile_energy_density(ε_q, mat)
            
            # 核心分支判断：
            if enforce_irreversibility
                # L型试件：取历史最大值
                driving_force[qp_count] = max(driving_force[qp_count], Ψ_plus)
            else
                # CT型试件：当前即为驱动力，允许减小
                driving_force[qp_count] = Ψ_plus
            end
            
            qp_count += 1
        end
    end
end

"""
    get_right_dofs(grid, dh_u, dir; tol=1e-12)
提取位于右边界的节点对应的特定方向位移自由度编号，用于计算反力。
dir = 1 代表水平方向(x), dir = 2 代表竖向(y)。
"""
function get_right_dofs(grid, dh_u, dir::Int; tol=1e-12)
    @assert dir == 1 || dir == 2 "dir 必须是 1 (x方向) 或 2 (y方向)"
    
    node_dofs_u = zeros(Int, 2, getnnodes(grid))
    for cell_id in 1:getncells(grid)
        cell = getcells(grid, cell_id)
        dofs = celldofs(dh_u, cell_id)
        for (local_node, node_id) in pairs(cell.nodes)
            node_dofs_u[1, node_id] = dofs[(local_node - 1) * 2 + 1]
            node_dofs_u[2, node_id] = dofs[(local_node - 1) * 2 + 2]
        end
    end

    # 找到最右侧边界的 x 坐标
    coords_x = [node.x[1] for node in grid.nodes]
    right_x = maximum(coords_x)
    
    # 筛选出位于右边界的节点
    right_nodes = findall(x -> isapprox(x, right_x; atol=tol), coords_x)
    
    # 根据传入的 dir 参数提取对应的自由度
    right_dofs = [node_dofs_u[dir, node_id] for node_id in right_nodes]
    
    return unique(right_dofs)
end

"""
    compute_reaction_forces()

计算F_{reaction}, 用于提取反力。
"""
function compute_reaction_forces(f_reac_dof, K_u, R_u, dh_u, dh_d, u_n, d_n, mat, cv_u, cv_d)
    # 计算整张网格的内力前必须清空矩阵和向量，否则会累加之前的虚假力或边界条件残留
    fill!(K_u.nzval, 0.0)
    fill!(R_u, 0.0)
    assemble_u!(K_u, R_u, dh_u, dh_d, u_n, d_n, mat, cv_u, cv_d)
        
    # 提取右边竖向位移自由度的反力并求和，得到总反力
    f_reac = sum(R_u[dof] for dof in f_reac_dof)
    return f_reac
end


"""
    compute_sav_scalars(dh_u, dh_d, u_n, d_n, u1, u2, d1, d2, r_n, E_nl_n, mat, cv_u, cv_d)

计算全局标量积分 A 和 B (公式 54, 55)，用于决定该时间步的比例因子 r^(n+1) = B/A。
"""
function compute_sav_scalars(
    dh_u::DofHandler, dh_d::DofHandler,
    u_n::Vector{Float64}, d_n::Vector{Float64},
    u1::Vector{Float64}, u2::Vector{Float64},
    d1::Vector{Float64}, d2::Vector{Float64},
    r_n::Float64, E_nl_n::Float64, nat::NumericalParams,
    mat::MaterialParams, cv_u::CellValues, cv_d::CellValues
)
    integral_A = 0.0
    integral_B = 0.0
    denom = 2.0 * sqrt(E_nl_n + nat.S0)

    for (cell_u, cell_d) in zip(CellIterator(dh_u), CellIterator(dh_d))
        reinit!(cv_u, cell_u)
        reinit!(cv_d, cell_d)
        
        # 提取相关场
        u_loc_n = u_n[celldofs(cell_u)]; d_loc_n = d_n[celldofs(cell_d)]
        u_loc_1 = u1[celldofs(cell_u)];  d_loc_1 = d1[celldofs(cell_d)]
        u_loc_2 = u2[celldofs(cell_u)];  d_loc_2 = d2[celldofs(cell_d)]
        
        for qp in 1:getnquadpoints(cv_u)
            dΩ = getdetJdV(cv_u, qp)
            
            # 计算应变与损伤值
            ε_n = function_symmetric_gradient(cv_u, qp, u_loc_n)
            d_q_n = function_value(cv_d, qp, d_loc_n)
            
            ε_2 = function_symmetric_gradient(cv_u, qp, u_loc_2)
            d_q_2 = function_value(cv_d, qp, d_loc_2)
            
            ε_1 = function_symmetric_gradient(cv_u, qp, u_loc_1)
            d_q_1 = function_value(cv_d, qp, d_loc_1)
            
            # 物理量评估
            σ_n = evaluate_damaged_stress(ε_n, d_q_n, mat)
            σ_tilde_n = σ_n - mat.C0 ⊡ ε_n     # σ̃ = σ - C₀:ε (公式 60)
            δ_d_ψ = -2.0 * (1.0 - d_q_n) * elastic_energy_density_tensile(ε_n, mat)

            # 积分子项 (公式 88, 89，使用 σ̃_n 而不是 σ_n)
            integral_A += (σ_tilde_n ⊡ ε_2 + δ_d_ψ * d_q_2) * dΩ
            integral_B += (σ_tilde_n ⊡ (ε_1 - ε_n) + δ_d_ψ * (d_q_1 - d_q_n)) * dΩ
        end
    end
    
    A = 1.0 - integral_A / denom
    B = r_n + integral_B / denom
    return A, B
end

# ==============================================================================
# 自定义求根函数 (Brent's Method)
# 接口仿照 Roots.jl，以便未来无缝切换
# ==============================================================================

struct Brent end

"""
    find_zero(f, bracket::Tuple{T, T}, ::Brent; tol=1e-12, max_iter=100) where {T<:Real}

使用 Brent 法寻找一维连续函数 f 的根。如果初始 bracket [a, b] 对应的函数值同号，
会自动向外扩展区间进行搜寻（步进法扩展策略）。
"""
function find_zero(f, bracket::Tuple{T, T}, ::Brent; tol=1e-12, max_iter=100) where {T<:Real}
    a, b = bracket
    fa = f(a)
    fb = f(b)

    # 自动扩展区间 Fallback 策略
    expand_iter = 0
    while fa * fb > 0 && expand_iter < 20
        # 如果同号，向外扩展 50% 的区间长度
        width = (b - a)
        if abs(fa) < abs(fb)
            a = a - 0.5 * width
            fa = f(a)
        else
            b = b + 0.5 * width
            fb = f(b)
        end
        expand_iter += 1
    end
    
    if fa * fb > 0
        error("Brent method failed to find a bracketing interval. Last bracket: [$a, $b], fa=$fa, fb=$fb")
    end

    if abs(fa) < abs(fb)
        a, b = b, a
        fa, fb = fb, fa
    end

    c = a
    fc = fa
    d = 0.0
    mflag = true

    for _ in 1:max_iter
        if fa != fc && fb != fc
            # 逆二次插值 (Inverse quadratic interpolation)
            s = a * fb * fc / ((fa - fb) * (fa - fc)) +
                b * fa * fc / ((fb - fa) * (fb - fc)) +
                c * fa * fb / ((fc - fa) * (fc - fb))
        else
            # 割线法 (Secant method)
            s = b - fb * (b - a) / (fb - fa)
        end

        # 判断是否需要回退到二分法 (Bisection method)
        # 这里条件 1 的原意是检查 s 是否在 (3a+b)/4 和 b 之间，但要注意 a 和 b 的大小关系
        if a < b
            cond1 = (s < (3 * a + b) / 4 || s > b)
        else
            cond1 = (s > (3 * a + b) / 4 || s < b)
        end
        cond2 = mflag && (abs(s - b) >= abs(b - c) / 2)
        cond3 = !mflag && (abs(s - b) >= abs(c - d) / 2)
        cond4 = mflag && (abs(b - c) < tol)
        cond5 = !mflag && (abs(c - d) < tol)

        if cond1 || cond2 || cond3 || cond4 || cond5
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

        if fb == 0.0 || abs(b - a) < tol
            return b
        end
    end

    @warn "Brent method reached max_iter without fully converging."
    return b
end