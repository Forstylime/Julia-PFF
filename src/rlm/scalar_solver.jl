struct RLMRootResult
    success::Bool
    code::Symbol
    message::String
    q::Float64
    discriminant::Float64
    discriminant_used::Float64
    residual::Float64
    candidates::Vector{Float64}
end

@inline function scalar_equation_residual(A, B, C, q, epsilon)
    numerator = abs(muladd(A, q^2, muladd(B, q, C)))
    denominator = abs(A) * q^2 + abs(B) * abs(q) + abs(C) + epsilon
    return numerator / denominator
end

function _coefficient_tolerance(A, B, C, tolerances)
    scale = max(abs(A), abs(B), abs(C))
    return tolerances.coefficient_abs + tolerances.coefficient_rel * scale
end

function _unique_candidates(candidates, relative_tolerance)
    unique_roots = Float64[]
    for root in candidates
        any(existing -> abs(root - existing) <=
            relative_tolerance * max(1.0, abs(root), abs(existing)), unique_roots) ||
            push!(unique_roots, root)
    end
    return unique_roots
end

"""
    solve_rlm_quadratic(A, B, C, q_previous, tolerances)

Solve the document's unmodified scalar equation with its cancellation-resistant
`qhat` formula. The result never silently accepts a negative discriminant,
nonpositive root, nonfinite root, or excessive original-equation residual.
"""
function solve_rlm_quadratic(
    A::Real,
    B::Real,
    C::Real,
    q_previous::Real,
    tolerances::RLMToleranceConfig,
)
    A64, B64, C64 = Float64(A), Float64(B), Float64(C)
    qn = Float64(q_previous)
    all(isfinite, (A64, B64, C64, qn)) || return RLMRootResult(
        false, :nonfinite_coefficients, "scalar coefficients and q^n must be finite",
        NaN, NaN, NaN, Inf, Float64[],
    )

    discriminant = muladd(-4.0 * A64, C64, B64^2)
    discriminant_scale = B64^2 + abs(4.0 * A64 * C64)
    discriminant_tolerance = tolerances.discriminant_abs +
                             tolerances.discriminant_rel * discriminant_scale
    if discriminant < -discriminant_tolerance
        return RLMRootResult(
            false,
            :negative_discriminant,
            "discriminant $discriminant is below -tolerance $(-discriminant_tolerance)",
            NaN,
            discriminant,
            discriminant,
            Inf,
            Float64[],
        )
    end
    discriminant_used = discriminant < 0.0 ? 0.0 : discriminant
    coefficient_tolerance = _coefficient_tolerance(A64, B64, C64, tolerances)
    candidates = Float64[]

    if abs(A64) <= coefficient_tolerance
        if abs(B64) <= coefficient_tolerance
            if abs(C64) <= coefficient_tolerance
                if isfinite(qn) && qn > tolerances.positive_root
                    push!(candidates, qn)
                else
                    return RLMRootResult(
                        false,
                        :indeterminate_without_positive_reference,
                        "the scalar equation is numerically 0=0 but q^n is not an admissible positive root",
                        NaN,
                        discriminant,
                        discriminant_used,
                        Inf,
                        candidates,
                    )
                end
            else
                return RLMRootResult(
                    false,
                    :inconsistent_degenerate_equation,
                    "A and B are numerically zero while C is not",
                    NaN,
                    discriminant,
                    discriminant_used,
                    Inf,
                    candidates,
                )
            end
        else
            push!(candidates, -C64 / B64)
        end
    else
        sqrt_discriminant = sqrt(discriminant_used)
        # Treat both +0.0 and -0.0 as the document's mathematical B=0 case.
        signed_sqrt = B64 == 0.0 ? sqrt_discriminant :
                      copysign(sqrt_discriminant, B64)
        qhat = -0.5 * (B64 + signed_sqrt)
        if qhat != 0.0
            push!(candidates, qhat / A64)
            push!(candidates, C64 / qhat)
        else
            # This is the repeated zero root (B=C=D=0). It is retained only so
            # the standard candidate filters can explicitly reject it.
            push!(candidates, -B64 / (2.0 * A64))
        end
    end

    candidates = _unique_candidates(candidates, tolerances.duplicate_root_rel)
    eligible = Tuple{Float64, Float64}[]
    best_positive_residual = Inf
    for root in candidates
        isfinite(root) || continue
        root > tolerances.positive_root || continue
        residual = scalar_equation_residual(
            A64,
            B64,
            C64,
            root,
            tolerances.scalar_denominator_epsilon,
        )
        best_positive_residual = min(best_positive_residual, residual)
        residual <= tolerances.scalar_residual || continue
        push!(eligible, (root, residual))
    end

    if isempty(eligible)
        return RLMRootResult(
            false,
            :no_admissible_positive_root,
            "no finite positive root satisfies the original scalar residual tolerance",
            NaN,
            discriminant,
            discriminant_used,
            best_positive_residual,
            candidates,
        )
    end

    selected = length(eligible) == 1 ? eligible[1] :
               eligible[argmin(abs(item[1] - qn) for item in eligible)]
    return RLMRootResult(
        true,
        :ok,
        "accepted positive root",
        selected[1],
        discriminant,
        discriminant_used,
        selected[2],
        candidates,
    )
end
