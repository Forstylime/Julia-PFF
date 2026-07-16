using JLD2

@load "data/jld2/rlm_results.jld2" disp force
peak_idx = argmax(abs.(force))
println("Max force: ", force[peak_idx])
println("Max force at disp: ", disp[peak_idx])
