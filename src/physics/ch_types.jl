# src/physics/ch_types.jl

"""
    CHParams

Cahn-Hilliard equation parameters.
"""
struct CHParams
    M::Float64       # Mobility
    ϵ::Float64       # Interface width
    α::Float64       # RLM relaxation parameter
    Δt::Float64      # Time step size
end

"""
    CHState

Stores the historical states needed for the RLM-PC algorithm.
"""
mutable struct CHState
    ϕ_n::Vector{Float64}
    ϕ_nm1::Vector{Float64}
    μ_n::Vector{Float64}
    q_n::Float64
    ϕ_star_n::Vector{Float64}
end

function CHState(ndofs_phi::Int)
    return CHState(
        zeros(ndofs_phi), # ϕ_n
        zeros(ndofs_phi), # ϕ_nm1
        zeros(ndofs_phi), # μ_n
        1.0,              # q_n (initially 1.0)
        zeros(ndofs_phi)  # ϕ_star_n
    )
end
