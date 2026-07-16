using JLD2
@load "data/jld2/staggered_results.jld2" disp force psi_energy gf_energy
println("Staggered Initial Force: ", force[1], " at ", disp[1])
println("Staggered First 10 forces: ", force[1:min(10, length(force))])
println("Staggered First 10 disps: ", disp[1:min(10, length(disp))])
println("Staggered Peak Force (Min force for compression): ", minimum(force), " at disp: ", disp[argmin(force)])
