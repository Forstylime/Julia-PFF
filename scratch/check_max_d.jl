using JLD2
@load "data/jld2/rlm_results.jld2" disp force
println("Max force: ", maximum(abs.(force)))

# Read VTK to check d
using WriteVTK
println("Check data/sims/rlm_amor directly")
