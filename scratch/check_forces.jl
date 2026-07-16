using JLD2
@load "data/jld2/rlm_results.jld2" disp force
println("Forces sample: ", force[1:50:end])
println("Max abs force: ", maximum(abs.(force)))
