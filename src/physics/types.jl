# ===================================================
# src/types.jl
# 存储相场断裂 SAV 数值格式的物理/数值参数以及状态容器
# ===================================================

using LinearAlgebra

"""
    MaterialParams{dim, T<:AbstractFloat}

存储固体材料的物理性质及裂纹传播的本构参数。
dim 为空间维度，T 为浮点精度。

# 用户可设关键字参数
- `E`: 杨氏模量
- `ν`: 泊松比
- `gc`: 临界能量释放率（断裂韧性）
- `l`: 正则化相场尺度参数 ℓ_c
- `η`: 相场演化的动力学粘性参数
- `ρ`: 材料质量密度
- `k`: 极小残留刚度参数，防止完全损坏时刚度矩阵奇异

# 自动计算字段
- `λ`, `μ`: Lamé 常数（由 E, ν 自动计算）
- `C0`: 无损材料的四阶弹性张量（由 λ, μ, dim 自动计算）
"""
struct MaterialParams{dim, T<:AbstractFloat}
    E::T
    ν::T
    λ::T
    μ::T
    gc::T
    l::T
    η::T
    ρ::T
    C0::SymmetricTensor{4, dim, T}
    k::T
end

function MaterialParams(;
    dim::Int = 2,
    E::T = 25840.0,
    ν::T = 0.18,
    gc::T = 0.65,
    l::T = 10.0,
    η::T = 0.05,
    ρ::T = 1.0,
    k::T = 1e-6,
) where {T<:AbstractFloat}
    λ_val = (E * ν) / ((1 + ν) * (1 - 2ν))  # plane strain
    μ_val = E / (2 * (1 + ν))

    @assert λ_val > 0 "Lamé 常数 λ 必须大于 0"
    @assert μ_val > 0 "剪切模量 μ 必须大于 0"
    @assert gc > 0 "断裂韧性 gc 必须大于 0"
    @assert l > 0 "长度尺度 l 必须大于 0"
    @assert η >= 0 "粘性参数 η 必须非负"
    @assert ρ > 0 "密度 ρ 必须大于 0"
    @assert 0 < k < 1e-2 "残留刚度 k 应为极小的正数"

    I2 = one(SymmetricTensor{2, dim, T})
    I4_sym = one(SymmetricTensor{4, dim, T})
    C0_val = λ_val * (I2 ⊗ I2) + 2μ_val * I4_sym

    return MaterialParams{dim, T}(E, ν, λ_val, μ_val, gc, l, η, ρ, C0_val, k)
end


"""
    NumericalParams{T<:AbstractFloat}

存储时间离散、显式不可逆惩罚以及 SAV 方法相关的算法控制参数。

# 成员变量
- `Δt::T`: 物理时间步长 \$\\Delta t\$
- `S0::T`: SAV 方法中的正定平移常数 \$S_0\$，确保非凸能量根号内的值恒大于 0
"""
Base.@kwdef struct NumericalParams{T<:AbstractFloat}
    Δt::T
    S0::T

    function NumericalParams(Δt::T, S0::T) where {T<:AbstractFloat}
        @assert Δt > 0 "时间步长 Δt 必须大于 0"
        @assert S0 >= 0 "SAV 平移常数 S0 必须 >= 0"
        return new{T}(Δt, S0)
    end
end


"""
    SimulationState{T<:AbstractFloat}

仿真状态容器。存储当前步以及前两步的历史物理量。

# 成员变量
- `u::Vector{T}`: 当前步位移场矢量 \$\\mathbf{u}^{n+1}\$
- `u_n::Vector{T}`: \$t^n\$ 步位移场矢量 \$\\mathbf{u}^n\$
- `u_nm1::Vector{T}`: \$t^{n-1}\$ 步位移场矢量 \$\\mathbf{u}^{n-1}\$ (用于位移外推)
- `v::Vector{T}`: 当前步速度场矢量 \$\\mathbf{v}^{n+1}\$ (用于动态 Newmark 步进)
- `v_n::Vector{T}`: \$t^n\$ 步速度场矢量 \$\\mathbf{v}^n\$
- `a::Vector{T}`: 当前步加速度场矢量 \$\\mathbf{a}^{n+1}\$
- `a_n::Vector{T}`: \$t^n\$ 步加速度场矢量 \$\\mathbf{a}^n\$

- `d::Vector{T}`: 当前步相场矢量 \$\\mathbf{d}^{n+1}\$
- `d_n::Vector{T}`: \$t^n\$ 步相场矢量 \$\\mathbf{d}^n\$
- `d_nm1::Vector{T}`: \$t^{n-1}\$ 步相场矢量 \$\\mathbf{d}^{n-1}\$ (用于相场外推)

- `r::T`: 当前步 SAV 标量辅助变量 \$r^{n+1}\$
- `r_n::T`: \$t^n\$ 步 SAV 标量辅助变量 \$r^n\$
"""
mutable struct SimulationState{T<:AbstractFloat}
    # 位移相关的时序向量 (全自由度向量)
    u::Vector{T}
    u_n::Vector{T}
    u_nm1::Vector{T}
    v::Vector{T}
    v_n::Vector{T}
    a::Vector{T}
    a_n::Vector{T}

    # 相场相关的时序向量 (全自由度向量)
    d::Vector{T}
    d_n::Vector{T}
    d_nm1::Vector{T}

    # SAV 全局时间相关标量
    r::T
    r_n::T

    # 默认构造函数：根据位移和相场的全局自由度数量进行初始化
    function SimulationState{T}(ndofs_u::Int, ndofs_d::Int, initial_r::T) where {T<:AbstractFloat}
        return new{T}(
            zeros(T, ndofs_u), # u
            zeros(T, ndofs_u), # u_n
            zeros(T, ndofs_u), # u_nm1
            zeros(T, ndofs_u), # v
            zeros(T, ndofs_u), # v_n
            zeros(T, ndofs_u), # a
            zeros(T, ndofs_u), # a_n
            
            zeros(T, ndofs_d), # d
            zeros(T, ndofs_d), # d_n
            zeros(T, ndofs_d), # d_nm1
            
            initial_r,         # r
            initial_r          # r_n
        )
    end
end


"""
    update_states!(state::SimulationState)

在一个时间步求解完全结束、准备向下一个物理时间步推进时，滚动更新仿真状态历史。
使 \$t^{n-1} \\leftarrow t^n\$, \$t^n \\leftarrow t^{n+1}\$。
"""
function update_states!(state::SimulationState{T}) where {T}
    # 滚动位移、速度与加速度历史
    copyto!(state.u_nm1, state.u_n)
    copyto!(state.u_n, state.u)
    copyto!(state.v_n, state.v)
    copyto!(state.a_n, state.a)

    # 滚动相场历史
    copyto!(state.d_nm1, state.d_n)
    copyto!(state.d_n, state.d)

    # 滚动 SAV 标量
    state.r_n = state.r
    return nothing
end