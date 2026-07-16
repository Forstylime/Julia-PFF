# src/fem/assembly.jl

using Ferrite
using Tensors
using LinearAlgebra
using SparseArrays

"""
    assemble_u!(K, R, dh_u, dh_d, u, d, mat, cv_u, cv_d)

组装位移场 u 的切线刚度矩阵 K 和残差向量 R (内力)。
这是非线性牛顿迭代的核心步骤。
"""
function assemble_u!(
    K::AbstractMatrix{T}, R::AbstractVector{T}, # 允许 T 类型 (Float64 或 Dual)
    dh_u::DofHandler, dh_d::DofHandler,
    u_global::AbstractVector{T}, d_global::AbstractVector, # 位移设为 T
    mat::MaterialParams, 
    cv_u::CellValues, cv_d::CellValues
) where T <: Real # 使用参数化类型 T
    # 初始化汇编器，自动把单元矩阵加到全局稀疏矩阵的对应位置
    assembler = start_assemble(K, R)
    
    n_basefuncs_u = getnbasefunctions(cv_u)
    
    # 获取每个单元的局部矩阵和向量缓存
    Ke = zeros(T, n_basefuncs_u, n_basefuncs_u)
    Re = zeros(T, n_basefuncs_u)
    
    # 遍历每一个网格单元 (同时遍历 u 和 d 的自由度布局)
    for (cell_u, cell_d) in zip(CellIterator(dh_u), CellIterator(dh_d))
        # 重新初始化形函数和雅可比矩阵
        reinit!(cv_u, cell_u)
        reinit!(cv_d, cell_d)
        
        # 提取当前单元上节点的全局自由度值
        u_loc = u_global[celldofs(cell_u)]
        d_loc = d_global[celldofs(cell_d)]
        
        fill!(Ke, 0.0)
        fill!(Re, 0.0)
        
        # 遍历单元内的积分点 (Gauss Quadrature points)
        for q_point in 1:getnquadpoints(cv_u)
            dΩ = getdetJdV(cv_u, q_point) # 雅可比行列式乘积分权重
            
            # 计算该积分点上的应变和相场值
            ε_q = function_symmetric_gradient(cv_u, q_point, u_loc)
            d_q = function_value(cv_d, q_point, d_loc)
            
            # 能量谱分解，在特征值重复时做二阶自动微分容易产生 NaN 的问题，
            # 已在 constitutive.jl 中添加微小扰动 ε_pert 以避免这个问题。
            ψ(ε) = elastic_energy_density(ε, d_q, mat)
            
            # 直接使用 Tensors.jl 求一阶导得到应力 σ，求二阶导得到四阶材料刚度张量 ℂ
            σ = Tensors.gradient(ψ, ε_q)
            ℂ = Tensors.hessian(ψ, ε_q)
            
            # 组装单元残差和刚度矩阵
            for i in 1:n_basefuncs_u
                δε = shape_symmetric_gradient(cv_u, q_point, i)
                # 残差: 内力向量 (σ:δε)
                Re[i] += (σ ⊡ δε) * dΩ 
                
                for j in 1:n_basefuncs_u
                    Δε = shape_symmetric_gradient(cv_u, q_point, j)
                    # 切线刚度: (δε:ℂ:Δε)
                    Ke[i, j] += (δε ⊡ ℂ ⊡ Δε) * dΩ
                end
            end
        end
        # 将单元矩阵推入全局矩阵
        assemble!(assembler, celldofs(cell_u), Ke, Re)
    end
end

"""
    assemble_d!(K, F, dh_d, H, mat, cv_d)

组装相场 d 的线性方程组 K*d = F
"""
function assemble_d!(
    K::SparseMatrixCSC, F::Vector{Float64}, 
    dh_d::DofHandler, H::Vector{Float64}, 
    mat::MaterialParams, cv_d::CellValues
)
    assembler = start_assemble(K, F)
    n_basefuncs_d = getnbasefunctions(cv_d)
    
    Ke = zeros(n_basefuncs_d, n_basefuncs_d)
    Fe = zeros(n_basefuncs_d)
    
    qp_count = 1
    for cell in CellIterator(dh_d)
        reinit!(cv_d, cell)
        
        fill!(Ke, 0.0)
        fill!(Fe, 0.0)
        
        for q_point in 1:getnquadpoints(cv_d)
            dΩ = getdetJdV(cv_d, q_point)
            H_q = H[qp_count] # 取出该积分点最新的历史变量
            
            # 将弱形式 (17) 重排为 A*d = B 的形式
            coef_d    = mat.gc / mat.l + 2.0 * H_q
            coef_grad = mat.gc * mat.l
            
            for i in 1:n_basefuncs_d
                δd  = shape_value(cv_d, q_point, i)
                ∇δd = shape_gradient(cv_d, q_point, i)
                
                # 右端项载荷向量
                Fe[i] += (2.0 * H_q * δd) * dΩ
                
                for j in 1:n_basefuncs_d
                    Δd  = shape_value(cv_d, q_point, j)
                    ∇Δd = shape_gradient(cv_d, q_point, j)
                    
                    # 刚度矩阵
                    Ke[i, j] += (coef_d * δd * Δd + coef_grad * (∇δd ⋅ ∇Δd)) * dΩ
                end
            end
            qp_count += 1
        end
        assemble!(assembler, celldofs(cell), Ke, Fe)
    end
end


"""
    assemble_SAV_LHS!(Kd, Ku, dh_u, dh_d, mat, cv_u, cv_d, dt, ρ, η)

在物理模拟开始前（t=0）仅执行一次。
组装恒定的相场左端矩阵 Kd 和位移场左端矩阵 Ku。
"""
function assemble_SAV_LHS!(
    Kd::SparseMatrixCSC, Ku::SparseMatrixCSC,
    dh_u::DofHandler, dh_d::DofHandler, 
    cv_u::CellValues, cv_d::CellValues,
    mat::MaterialParams,
    nat::NumericalParams
)
    # 1. 组装相场 LHS (Kd)
    assembler_d = start_assemble(Kd)
    n_base_d = getnbasefunctions(cv_d)
    Ke_d = zeros(n_base_d, n_base_d)
    
    # 这里的系数完全不含 H，因此是恒定常数
    coef_d = mat.η / nat.Δt + mat.gc / mat.l
    coef_grad = mat.gc * mat.l
    
    for cell in CellIterator(dh_d)
        reinit!(cv_d, cell)
        fill!(Ke_d, 0.0)
        for qp in 1:getnquadpoints(cv_d)
            dΩ = getdetJdV(cv_d, qp)
            for i in 1:n_base_d
                δd = shape_value(cv_d, qp, i)
                ∇δd = shape_gradient(cv_d, qp, i)
                for j in 1:n_base_d
                    Δd = shape_value(cv_d, qp, j)
                    ∇Δd = shape_gradient(cv_d, qp, j)
                    Ke_d[i, j] += (coef_d * δd * Δd + coef_grad * (∇δd ⋅ ∇Δd)) * dΩ
                end
            end
        end
        assemble!(assembler_d, celldofs(cell), Ke_d)
    end

    # 2. 组装位移场 LHS (Ku)
    assembler_u = start_assemble(Ku)
    n_base_u = getnbasefunctions(cv_u)
    Ke_u = zeros(n_base_u, n_base_u)
    
    # C0 是无损材料的弹性刚度张量 (第四阶张量)
    C0 = mat.C0 
    coef_mass = mat.ρ / (nat.Δt^2)
    
    for cell in CellIterator(dh_u)
        reinit!(cv_u, cell)
        fill!(Ke_u, 0.0)
        for qp in 1:getnquadpoints(cv_u)
            dΩ = getdetJdV(cv_u, qp)
            for i in 1:n_base_u
                δu = shape_value(cv_u, qp, i) # 用于质量矩阵
                δε = shape_symmetric_gradient(cv_u, qp, i) # 用于刚度矩阵
                for j in 1:n_base_u
                    Δu = shape_value(cv_u, qp, j)
                    Δε = shape_symmetric_gradient(cv_u, qp, j)
                    
                    # 动力学质量项 + 无损线性弹性刚度项
                    Ke_u[i, j] += (coef_mass * (δu ⋅ Δu) + (δε ⊡ C0 ⊡ Δε)) * dΩ
                end
            end
        end
        assemble!(assembler_u, celldofs(cell), Ke_u)
    end
end


"""
    assemble_SAV_RHS!(fd1, fd2, fu1, fu2, dh_u, dh_d, u_n, v_n, d_n, ξ_n, E_nl_n, mat, cv_u, cv_d, dt, ρ, η)

在每个时间步迭代中调用。高效组装四个基解的右端项向量 fd1, fd2, fu1, fu2。
"""
function assemble_SAV_RHS!(
    fd1::Vector{Float64}, fd2::Vector{Float64},
    fu1::Vector{Float64}, fu2::Vector{Float64},
    dh_u::DofHandler, dh_d::DofHandler,
    u_n::Vector{Float64}, v_n::Vector{Float64}, d_n::Vector{Float64},
    ξ_n::Vector{Float64}, E_nl_n::Float64,
    cv_u::CellValues, cv_d::CellValues,
    mat::MaterialParams,
    nat::NumericalParams
)
    n_base_d = getnbasefunctions(cv_d)
    n_base_u = getnbasefunctions(cv_u)

    fe_d1 = zeros(n_base_d); fe_d2 = zeros(n_base_d)
    fe_u1 = zeros(n_base_u); fe_u2 = zeros(n_base_u)

    # 标量分母
    denom = sqrt(E_nl_n + nat.S0)
    C0 = mat.C0

    qp_count = 1
    for (cell_u, cell_d) in zip(CellIterator(dh_u), CellIterator(dh_d))
        reinit!(cv_u, cell_u)
        reinit!(cv_d, cell_d)

        u_loc = u_n[celldofs(cell_u)]
        v_loc = v_n[celldofs(cell_u)]
        d_loc = d_n[celldofs(cell_d)]

        fill!(fe_d1, 0.0); fill!(fe_d2, 0.0)
        fill!(fe_u1, 0.0); fill!(fe_u2, 0.0)

        for qp in 1:getnquadpoints(cv_u)
            dΩ = getdetJdV(cv_u, qp)

            # 基础物理量计算
            d_q = function_value(cv_d, qp, d_loc)
            ε_q = function_symmetric_gradient(cv_u, qp, u_loc)
            u_q = function_value(cv_u, qp, u_loc)
            v_q = function_value(cv_u, qp, v_loc)

            # 计算当前步的物理受损应力 σ_n
            σ_n = evaluate_damaged_stress(ε_q, d_q, mat)
            ξ_q = ξ_n[qp_count]

            # 1. 组装相场基解的 RHS
            for i in 1:n_base_d
                δd = shape_value(cv_d, qp, i)
                fe_d1[i] += (mat.η / nat.Δt * d_q * δd) * dΩ
                fe_d2[i] += (2.0 * (1.0 - d_q) * ξ_q * δd) * dΩ
            end

            # 2. 组装位移场基解的 RHS
            for i in 1:n_base_u
                δu = shape_value(cv_u, qp, i)
                δε = shape_symmetric_gradient(cv_u, qp, i)

                # u1 动力学惯性项
                inertia_q = mat.ρ / (nat.Δt^2) * (u_q + nat.Δt * v_q)
                fe_u1[i] += (inertia_q ⋅ δu) * dΩ

                # u2 偏置应力弱形式散度项
                σ_diff = σ_n - C0 ⊡ ε_q
                fe_u2[i] += (- (δε ⊡ σ_diff) / denom) * dΩ
            end
            qp_count += 1
        end
        assemble!(fd1, celldofs(cell_d), fe_d1)
        assemble!(fd2, celldofs(cell_d), fe_d2)
        assemble!(fu1, celldofs(cell_u), fe_u1)
        assemble!(fu2, celldofs(cell_u), fe_u2)
    end
end

"""
    assemble_SAV_QS_LHS!(Kd, Ku, dh_u, dh_d, mat, cv_u, cv_d, dt, ρ, η)

在物理模拟开始前（t=0）仅执行一次。
组装恒定的相场左端矩阵 Kd 和位移场左端矩阵 Ku。这是准静态版本（无惯性项）。
"""
function assemble_SAV_QS_LHS!(
    Kd::SparseMatrixCSC, Ku::SparseMatrixCSC,
    dh_u::DofHandler, dh_d::DofHandler, 
    cv_u::CellValues, cv_d::CellValues,
    mat::MaterialParams,
    nat::NumericalParams
)
    # 1. 组装相场 LHS (Kd)
    assembler_d = start_assemble(Kd)
    n_base_d = getnbasefunctions(cv_d)
    Ke_d = zeros(n_base_d, n_base_d)
    
    coef_d = mat.η / nat.Δt + mat.gc / mat.l
    coef_grad = mat.gc * mat.l
    
    for cell in CellIterator(dh_d)
        reinit!(cv_d, cell)
        fill!(Ke_d, 0.0)
        for qp in 1:getnquadpoints(cv_d)
            dΩ = getdetJdV(cv_d, qp)
            for i in 1:n_base_d
                δd = shape_value(cv_d, qp, i)
                ∇δd = shape_gradient(cv_d, qp, i)
                for j in 1:n_base_d
                    Δd = shape_value(cv_d, qp, j)
                    ∇Δd = shape_gradient(cv_d, qp, j)
                    Ke_d[i, j] += (coef_d * δd * Δd + coef_grad * (∇δd ⋅ ∇Δd)) * dΩ
                end
            end
        end
        assemble!(assembler_d, celldofs(cell), Ke_d)
    end

    # 2. 组装位移场 LHS (Ku)
    assembler_u = start_assemble(Ku)
    n_base_u = getnbasefunctions(cv_u)
    Ke_u = zeros(n_base_u, n_base_u)
    
    # C0 是无损材料的弹性刚度张量 (第四阶张量)
    C0 = mat.C0 
    
    for cell in CellIterator(dh_u)
        reinit!(cv_u, cell)
        fill!(Ke_u, 0.0)
        for qp in 1:getnquadpoints(cv_u)
            dΩ = getdetJdV(cv_u, qp)
            for i in 1:n_base_u
                δε = shape_symmetric_gradient(cv_u, qp, i) # 用于刚度矩阵
                for j in 1:n_base_u
                    Δε = shape_symmetric_gradient(cv_u, qp, j)
                    
                    # 仅有无损线性弹性刚度项（去除动力学质量项）
                    Ke_u[i, j] += (δε ⊡ C0 ⊡ Δε) * dΩ
                end
            end
        end
        assemble!(assembler_u, celldofs(cell), Ke_u)
    end
end


"""
    assemble_SAV_QS_RHS!(fd1, fd2, fu1, fu2, dh_u, dh_d, u_n, d_n, ξ_n, E_nl_n, cv_u, cv_d, mat, nat)

在每个时间步迭代中调用。准静态版本：高效组装四个基解的右端项向量 fd1, fd2, fu1, fu2（无速度和惯性项）。
"""
function assemble_SAV_QS_RHS!(
    fd1::Vector{Float64}, fd2::Vector{Float64},
    fu1::Vector{Float64}, fu2::Vector{Float64},
    dh_u::DofHandler, dh_d::DofHandler,
    u_n::Vector{Float64}, d_n::Vector{Float64},
    ξ_n::Vector{Float64}, E_nl_n::Float64,
    cv_u::CellValues, cv_d::CellValues,
    mat::MaterialParams,
    nat::NumericalParams
)
    n_base_d = getnbasefunctions(cv_d)
    n_base_u = getnbasefunctions(cv_u)

    fe_d1 = zeros(n_base_d); fe_d2 = zeros(n_base_d)
    fe_u1 = zeros(n_base_u); fe_u2 = zeros(n_base_u)

    # 标量分母
    denom = sqrt(E_nl_n + nat.S0)
    C0 = mat.C0

    qp_count = 1
    for (cell_u, cell_d) in zip(CellIterator(dh_u), CellIterator(dh_d))
        reinit!(cv_u, cell_u)
        reinit!(cv_d, cell_d)

        u_loc = u_n[celldofs(cell_u)]
        d_loc = d_n[celldofs(cell_d)]

        fill!(fe_d1, 0.0); fill!(fe_d2, 0.0)
        fill!(fe_u1, 0.0); fill!(fe_u2, 0.0)

        for qp in 1:getnquadpoints(cv_u)
            dΩ = getdetJdV(cv_u, qp)

            # 基础物理量计算
            d_q = function_value(cv_d, qp, d_loc)
            ε_q = function_symmetric_gradient(cv_u, qp, u_loc)

            # 计算当前步的物理受损应力 σ_n
            σ_n = evaluate_damaged_stress(ε_q, d_q, mat)
            ξ_q = ξ_n[qp_count]

            # 1. 组装相场基解的 RHS
            for i in 1:n_base_d
                δd = shape_value(cv_d, qp, i)
                fe_d1[i] += (mat.η / nat.Δt * d_q * δd) * dΩ
                fe_d2[i] += (2.0 * (1.0 - d_q) * ξ_q * δd) * dΩ
            end

            # 2. 组装位移场基解的 RHS
            for i in 1:n_base_u
                δε = shape_symmetric_gradient(cv_u, qp, i)

                # u1 动力学惯性项 (准静态下已被移除，假定无体力，这里是0。如果有体力则应添加。)
                # fe_u1[i] += 0.0 
                
                # u2 偏置应力弱形式散度项
                σ_diff = σ_n - C0 ⊡ ε_q
                fe_u2[i] += (- (δε ⊡ σ_diff) / denom) * dΩ
            end
            qp_count += 1
        end
        assemble!(fd1, celldofs(cell_d), fe_d1)
        assemble!(fd2, celldofs(cell_d), fe_d2)
        assemble!(fu1, celldofs(cell_u), fe_u1)
        assemble!(fu2, celldofs(cell_u), fe_u2)
    end
end