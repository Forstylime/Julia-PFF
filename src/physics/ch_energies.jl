# src/physics/ch_energies.jl

"""
    ch_F(ϕ)

Double-well potential free energy density.
"""
function ch_F(ϕ::Float64)
    return 0.25 * (ϕ^2 - 1.0)^2
end

"""
    ch_f(ϕ)

Derivative of the double-well potential free energy density.
"""
function ch_f(ϕ::Float64)
    return ϕ^3 - ϕ
end
