"""Numerical Miehe split of a two-dimensional symmetric strain tensor."""
struct MieheSplit2D
    epsilon_plus::SymmetricTensor{2, 2, Float64, 3}
    epsilon_minus::SymmetricTensor{2, 2, Float64, 3}
    principal_min::Float64
    principal_max::Float64
    principal_tolerance::Float64
    repeated_tolerance::Float64
end

@inline degradation(d::Real, kappa::Real) = (1.0 - kappa) * (1.0 - d)^2 + kappa
@inline degradation_derivative(d::Real, kappa::Real) = -2.0 * (1.0 - kappa) * (1.0 - d)

"""
    miehe_split_2d(epsilon, tolerances)

Analytic 2x2 spectral split. `hypot` and an invariant projector formula avoid
storing or differentiating individual eigenvectors. Principal strains inside
the configured zero band are treated as numerical zero.
"""
function miehe_split_2d(
    epsilon::SymmetricTensor{2, 2, T},
    tolerances::RLMToleranceConfig,
) where {T<:Real}
    a = Float64(epsilon[1, 1])
    b = Float64(epsilon[1, 2])
    c = Float64(epsilon[2, 2])
    mean_strain = 0.5 * (a + c)
    dev_normal = 0.5 * (a - c)
    radius = hypot(dev_normal, b)

    lambda_min_raw = mean_strain - radius
    lambda_max_raw = mean_strain + radius
    principal_scale = max(abs(lambda_min_raw), abs(lambda_max_raw))
    zero_tol = tolerances.principal_zero_abs +
               tolerances.principal_zero_rel * principal_scale
    repeated_tol = tolerances.repeated_eigen_abs +
                   tolerances.repeated_eigen_rel * principal_scale

    # Values inside the zero band are inactive in the positive branch. They are
    # retained in the complementary branch below, so epsilon_plus +
    # epsilon_minus and psi_plus + psi_minus still match the unsplit tensor and
    # energy to roundoff.
    p_min = lambda_min_raw > zero_tol ? lambda_min_raw : 0.0
    p_max = lambda_max_raw > zero_tol ? lambda_max_raw : 0.0

    epsilon_float = SymmetricTensor{2, 2, Float64}((a, b, c))
    if radius <= repeated_tol && lambda_min_raw > zero_tol
        epsilon_plus = epsilon_float
    elseif radius <= repeated_tol && lambda_max_raw <= zero_tol
        epsilon_plus = zero(epsilon_float)
    elseif radius == 0.0
        epsilon_plus = mean_strain > zero_tol ? epsilon_float : zero(epsilon_float)
    else
        p_average = 0.5 * (p_max + p_min)
        p_beta = (p_max - p_min) / (2.0 * radius)
        epsilon_plus = SymmetricTensor{2, 2, Float64}((
            p_average + p_beta * dev_normal,
            p_beta * b,
            p_average - p_beta * dev_normal,
        ))
    end
    epsilon_minus = epsilon_float - epsilon_plus

    return MieheSplit2D(
        epsilon_plus,
        epsilon_minus,
        lambda_min_raw,
        lambda_max_raw,
        zero_tol,
        repeated_tol,
    )
end

@inline function _trace_parts(epsilon, zero_tolerance)
    trace_value = Float64(tr(epsilon))
    trace_plus = trace_value > zero_tolerance ? trace_value : 0.0
    return trace_plus, trace_value - trace_plus
end

"""Return `(psi_plus, psi_minus, sigma_plus, sigma_minus, split)`."""
function miehe_response_2d(
    epsilon::SymmetricTensor{2, 2, T},
    material::RLMMaterialConfig,
    tolerances::RLMToleranceConfig,
) where {T<:Real}
    split = miehe_split_2d(epsilon, tolerances)
    trace_plus, trace_minus = _trace_parts(epsilon, split.principal_tolerance)
    lambda = lame_lambda(material)
    mu = lame_mu(material)
    identity2 = one(SymmetricTensor{2, 2, Float64})

    psi_plus = 0.5 * lambda * trace_plus^2 +
               mu * (split.epsilon_plus ⊡ split.epsilon_plus)
    psi_minus = 0.5 * lambda * trace_minus^2 +
                mu * (split.epsilon_minus ⊡ split.epsilon_minus)
    sigma_plus = lambda * trace_plus * identity2 + 2.0 * mu * split.epsilon_plus
    sigma_minus = lambda * trace_minus * identity2 + 2.0 * mu * split.epsilon_minus
    return psi_plus, psi_minus, sigma_plus, sigma_minus, split
end

function plane_strain_elasticity_tensor(material::RLMMaterialConfig)
    identity2 = one(SymmetricTensor{2, 2, Float64})
    identity4 = one(SymmetricTensor{4, 2, Float64})
    return lame_lambda(material) * (identity2 ⊗ identity2) +
           2.0 * lame_mu(material) * identity4
end
