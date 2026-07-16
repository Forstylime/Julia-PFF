# ==============================================================================
# src/constitutive.jl
# 积分点层面的局部物理计算：谱分解、正负能量分裂、退化 Cauchy 应力
# ==============================================================================

using Tensors
using LinearAlgebra

# ------------------------------------------------------------------------------
# 1. 物理模式特征定义 (Traits for Elastic Mode)
# ------------------------------------------------------------------------------
abstract type ElasticMode end

"""平面应变模式（常用于 2D 厚板、裂纹垂直于平面的标准断裂基准题）"""
struct PlaneStrain <: ElasticMode end

"""平面应力模式（常用于 2D 薄板）"""
struct PlaneStress <: ElasticMode end

"""三维实体模式"""
struct ThreeDimensional <: ElasticMode end


# ------------------------------------------------------------------------------
# 2. 适应不同物理模式的 Lamé 常数调节器
# ------------------------------------------------------------------------------
"""
    effective_lambda(::ElasticMode, params::MaterialParams{T}) where {T}

根据物理分析模式返回等效的第一 Lamé 常数 \$\\lambda^*\$。
- 对于 `PlaneStrain` 和 `ThreeDimensional`，\$\\lambda^* = \\lambda\$
- 对于 `PlaneStress`，\$\\lambda^* = \\frac{2\\lambda\\mu}{\\lambda + 2\\mu}\$
"""
@inline effective_lambda(::PlaneStrain, params::MaterialParams{dim}) where {dim} = params.λ
@inline effective_lambda(::ThreeDimensional, params::MaterialParams{dim}) where {dim} = params.λ
@inline effective_lambda(::PlaneStress, params::MaterialParams{dim}) where {dim} =
    (2 * params.λ * params.μ) / (params.λ + 2 * params.μ)


# ------------------------------------------------------------------------------
# 3. 核心应变张量谱分解 (Spectral Decomposition)
# ------------------------------------------------------------------------------
"""
    spectral_decomposition(ε::SymmetricTensor{2, dim, T}) where {dim, T}

对对称应变张量 \$\\boldsymbol{\\varepsilon}\$ 进行谱分解，返回拉伸和压缩对应的正负应变部分 \$\\boldsymbol{\\varepsilon}_+\$ 和 \$\\boldsymbol{\\varepsilon}_-\$。
该实现基于 Miehe 谱分解公式 (Eq. 11)。
"""
function spectral_decomposition(ε::SymmetricTensor{2, dim, T}) where {dim, T}
    # 利用 Tensors.jl 的高性能特征值与特征向量求解器
    E = eigen(ε)
    vals = E.values     # 升序排列的特征值向量 Vec{dim, T}
    vecs = E.vectors    # 特征向量张量 Tensor{2, dim, T}，列向量为特征向量

    ε_plus = zero(SymmetricTensor{2, dim, T})
    ε_minus = zero(SymmetricTensor{2, dim, T})

    # 遍历特征空间组装拉伸与压缩张量
    for i in 1:dim
        val = vals[i]
        # 提取第 i 个特征向量 (Tensors.jl 的列提取语法)
        v = vecs[:, i] 
        
        # 麦考利夹（Macaulay brackets）操作
        val_plus = max(val, zero(T))
        val_minus = val - val_plus # 精确互补，防止浮点微小间隙

        # 计算特征方向的对称外积
        v_otimes_v = symmetric(v ⊗ v)

        ε_plus += val_plus * v_otimes_v
        ε_minus += val_minus * v_otimes_v
    end

    return ε_plus, ε_minus
end


# ------------------------------------------------------------------------------
# 4. 正负应变能密度计算 (Elastic Energy Densities)
# ------------------------------------------------------------------------------
"""
    elastic_energy_densities(mode::ElasticMode, ε::SymmetricTensor{2, dim, T}, params::MaterialParams{T}) where {dim, T}

计算无损伤状态下的拉伸应变能密度 \$\\psi_0^+\$ 和压缩应变能密度 \$\\psi_0^-\$ (Eq. 12)。
- `mode`: 选择 `PlaneStrain()`、`PlaneStress()` 或 `ThreeDimensional()`
"""
function elastic_energy_densities(mode::ElasticMode, ε::SymmetricTensor{2, dim, T}, params::MaterialParams{dim}) where {dim, T}
    # 1. 对应变张量进行拉压分解
    ε_plus, ε_minus = spectral_decomposition(ε)

    # 2. 计算迹（Trace）及其正负分裂
    tr_ε = tr(ε)
    tr_ε_plus = max(tr_ε, zero(T))
    tr_ε_minus = tr_ε - tr_ε_plus

    # 3. 提取特征 Lamé 参数
    λ_eff = effective_lambda(mode, params)
    μ = params.μ

    # 4. 计算应变能密度分量 (Eq. 12)
    # 注：Tensors.jl 中双收缩运算符 ⊡ (输入 \ddot) 代表内积，ε_plus ⊡ ε_plus 等价于 tr(ε_plus^2)
    ψ_plus = 0.5 * λ_eff * tr_ε_plus^2 + μ * (ε_plus ⊡ ε_plus)
    ψ_minus = 0.5 * λ_eff * tr_ε_minus^2 + μ * (ε_minus ⊡ ε_minus)

    return ψ_plus, ψ_minus
end


# ------------------------------------------------------------------------------
# 5. 应力张量计算 (Cauchy Stress)
# ------------------------------------------------------------------------------
"""
    elastic_stresses_0(mode::ElasticMode, ε::SymmetricTensor{2, dim, T}, params::MaterialParams{T}) where {dim, T}

计算无损伤状态下的拉伸应力张量 \$\\boldsymbol{\\sigma}_0^+\$ 和压缩应力张量 \$\\boldsymbol{\\sigma}_0^-\$。
"""
function elastic_stresses_0(mode::ElasticMode, ε::SymmetricTensor{2, dim, T}, params::MaterialParams{dim}) where {dim, T}
    ε_plus, ε_minus = spectral_decomposition(ε)

    tr_ε = tr(ε)
    tr_ε_plus = max(tr_ε, zero(T))
    tr_ε_minus = tr_ε - tr_ε_plus

    λ_eff = effective_lambda(mode, params)
    μ = params.μ
    I2 = one(SymmetricTensor{2, dim, T}) # 2阶单位对称张量

    σ0_plus = λ_eff * tr_ε_plus * I2 + 2 * μ * ε_plus
    σ0_minus = λ_eff * tr_ε_minus * I2 + 2 * μ * ε_minus

    return σ0_plus, σ0_minus
end


"""
    cauchy_stress(mode::ElasticMode, ε::SymmetricTensor{2, dim, T}, d::T, params::MaterialParams{T}) where {dim, T}

计算经典的（非 SAV）退化 Cauchy 应力张量 \$\\boldsymbol{\\sigma}\$ (Eq. 17)：
\$\\boldsymbol{\\sigma} = g(d)\\boldsymbol{\\sigma}_0^+ + \\boldsymbol{\\sigma}_0^-\$
"""
function cauchy_stress(mode::ElasticMode, ε::SymmetricTensor{2, dim, T}, d::T, params::MaterialParams{dim}) where {dim, T}
    σ0_plus, σ0_minus = elastic_stresses_0(mode, ε, params)
    g_d = (one(T) - d)^2 + params.k # 退化函数 g(d) = (1-d)^2 + k
    return g_d * σ0_plus + σ0_minus
end


"""
    elastic_energy_density(ε::SymmetricTensor{2,2}, d::Real, mat::MaterialParams{T})

计算单个积分点上的纯弹性能密度 Ψ(ε, d)。
结合了微裂纹闭合效应 (MCR-effect)。
"""
function elastic_energy_density(ε::SymmetricTensor{2,dim,T}, d::Real, mat::MaterialParams{dim}) where {dim,T}
    ε_pert = SymmetricTensor{2, dim, T}((T(1e-14), T(0.0), T(-1e-14)))
    ε = ε + ε_pert

    Ψ_plus, Ψ_minus = elastic_energy_densities(PlaneStrain(), ε, mat)

    g_d = (1.0 - d)^2 + mat.k

    return g_d * Ψ_plus + Ψ_minus
end


# ------------------------------------------------------------------------------
# 6. SAV 辅助函数
# ------------------------------------------------------------------------------
"""
    evaluate_damaged_stress(ε, d, mat)

计算受损 Cauchy 应力: σ = g(d)σ₀⁺ + σ₀⁻
`cauchy_stress` 的便捷封装，默认使用 PlaneStrain 模式。
"""
function evaluate_damaged_stress(
    ε::SymmetricTensor{2, dim, T}, d::T, mat::MaterialParams{dim}
) where {dim, T}
    return cauchy_stress(PlaneStrain(), ε, d, mat)
end

"""
    tensile_energy_density(ε, mat)

返回无损拉伸应变能密度 ψ₀⁺(ε)。
"""
function tensile_energy_density(
    ε::SymmetricTensor{2, dim, T}, mat::MaterialParams{dim}
) where {dim, T}
    ψ_plus, _ = elastic_energy_densities(PlaneStrain(), ε, mat)
    return ψ_plus
end

"""
    elastic_energy_density_tensile(ε, mat)

`tensile_energy_density` 的别名，用于 SAV 标量计算。
"""
function elastic_energy_density_tensile(
    ε::SymmetricTensor{2, dim, T}, mat::MaterialParams{dim}
) where {dim, T}
    return tensile_energy_density(ε, mat)
end

# ------------------------------------------------------------------------------
# 7. RLM Smoothed Amor Split 辅助函数
# ------------------------------------------------------------------------------

"""
    smoothed_amor_energy(mode::ElasticMode, ε::SymmetricTensor{2, dim, T}, ϵ::T, mat::MaterialParams{dim})

计算基于 ϵ-光滑化正则化的 Amor 劈裂应变能密度 ψ₀⁺ 和 ψ₀⁻。
"""
function smoothed_amor_energy(
    mode::ElasticMode, ε::SymmetricTensor{2, dim, T}, ϵ::T, mat::MaterialParams{dim}
) where {dim, T}
    λ_eff = effective_lambda(mode, mat)
    μ = mat.μ
    # 按照 Amor 劈裂，K = λ + 2μ/dim
    K = λ_eff + 2 * μ / dim

    tr_ε = tr(ε)
    # 光滑化 Macaulay 括号
    tr_plus_ϵ  = 0.5 * (tr_ε + sqrt(tr_ε^2 + ϵ^2))
    tr_minus_ϵ = 0.5 * (tr_ε - sqrt(tr_ε^2 + ϵ^2))

    # 偏应变张量
    I2 = one(SymmetricTensor{2, dim, T})
    ε_dev = ε - (tr_ε / dim) * I2

    ψ_plus  = 0.5 * K * tr_plus_ϵ^2 + μ * (ε_dev ⊡ ε_dev)
    ψ_minus = 0.5 * K * tr_minus_ϵ^2

    return ψ_plus, ψ_minus
end

"""
    smoothed_amor_stress(mode::ElasticMode, ε::SymmetricTensor{2, dim, T}, ϵ::T, mat::MaterialParams{dim})

计算基于 ϵ-光滑化正则化的 Amor 劈裂应力张量 σ₀⁺ 和 σ₀⁻。
"""
function smoothed_amor_stress(
    mode::ElasticMode, ε::SymmetricTensor{2, dim, T}, ϵ::T, mat::MaterialParams{dim}
) where {dim, T}
    λ_eff = effective_lambda(mode, mat)
    μ = mat.μ
    K = λ_eff + 2 * μ / dim

    tr_ε = tr(ε)
    denom = sqrt(tr_ε^2 + ϵ^2)
    
    # 避免除以 0
    ratio = denom > 1e-15 ? (tr_ε / denom) : zero(T)

    tr_plus_ϵ  = 0.5 * (tr_ε + denom)
    tr_minus_ϵ = 0.5 * (tr_ε - denom)

    # 导数项： 0.5 * (1 ± tr_ε / sqrt(tr_ε^2 + ϵ^2))
    diff_plus  = 0.5 * (1.0 + ratio)
    diff_minus = 0.5 * (1.0 - ratio)

    I2 = one(SymmetricTensor{2, dim, T})
    ε_dev = ε - (tr_ε / dim) * I2

    σ0_plus  = K * tr_plus_ϵ * diff_plus * I2 + 2 * μ * ε_dev
    σ0_minus = K * tr_minus_ϵ * diff_minus * I2

    return σ0_plus, σ0_minus
end

"""
    evaluate_rlm_damaged_stress(ε::SymmetricTensor{2, dim, T}, d::T, ϵ::T, mat::MaterialParams{dim})

计算 RLM 框架下受损的 Cauchy 应力：σ = g(d)σ₀⁺(ϵ) + σ₀⁻(ϵ)
"""
function evaluate_rlm_damaged_stress(
    ε::SymmetricTensor{2, dim, T}, d::T, ϵ::T, mat::MaterialParams{dim}
) where {dim, T}
    σ0_plus, σ0_minus = smoothed_amor_stress(PlaneStrain(), ε, ϵ, mat)
    g_d = (one(T) - d)^2 + mat.k
    return g_d * σ0_plus + σ0_minus
end