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