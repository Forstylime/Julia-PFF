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
    C0::SymmetricTensor{4, dim, T}
    k::T
end

function MaterialParams(;
    dim::Integer = 2,
    E::Real = 25_840.0,
    ν::Real = 0.18,
    gc::Real = 0.65,
    l::Real = 10.0,
    k::Real = 1.0e-6,
)
    dim > 0 || throw(ArgumentError("dim must be positive"))

    T = promote_type(
        typeof(float(E)),
        typeof(float(ν)),
        typeof(float(gc)),
        typeof(float(l)),
        typeof(float(k)),
    )
    E, ν, gc, l, k = T.((E, ν, gc, l, k))

    E > 0 || throw(ArgumentError("Young's modulus E must be positive"))
    0 < ν < 0.5 || throw(ArgumentError("Poisson's ratio ν must satisfy 0 < ν < 0.5"))
    gc > 0 || throw(ArgumentError("fracture toughness gc must be positive"))
    l > 0 || throw(ArgumentError("length scale l must be positive"))
    0 < k < 1.0e-2 || throw(ArgumentError("residual stiffness k must satisfy 0 < k < 1e-2"))

    λ_val = (E * ν) / ((1 + ν) * (1 - 2ν)) # Plane strain
    μ_val = E / (2 * (1 + ν))

    I2 = one(SymmetricTensor{2, dim, T})
    I4_sym = one(SymmetricTensor{4, dim, T})
    C0_val = λ_val * (I2 ⊗ I2) + 2μ_val * I4_sym

    return MaterialParams{Int(dim),T}(E, ν, λ_val, μ_val, gc, l, C0_val, k)
end
