using JLD2
f = jldopen("data/jld2/rlm_results.jld2")
println(keys(f))
close(f)
