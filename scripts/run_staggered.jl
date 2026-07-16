using PffSAV
using Ferrite
using FerriteGmsh

# 快速示例：从项目根目录运行
#   julia --project=. scripts/run_staggered.jl

mat = MaterialParams(;
    dim = 2,
    E  = 25840.0,
    ν  = 0.18,
    gc = 0.65,
    l  = 10.0,
    η  = 1e-5,
    ρ  = 2.4e-9,
    k  = 1e-8
)

setup = setup_tension(
    msh_file = "data/mesh/l_shape.msh",
    final_displacement = -0.8, # 施加位移
    fixed_face = "top", # 固定边界名称
    dir = 2, # 位移施加方向，x方向为1，y方向为2
)


disp, force, psi_energy, gf_energy = solve_staggered(setup, mat;
    n_steps = 100,
    max_iter = 2000, 
    enforce_irreversibility = true,
) # 实际最终位移为 setup.final_displacement

# 可视化载荷-位移曲线 and 能量演变
mkpath("data/plots")
mkpath("data/jld2")
using CairoMakie

disp_plot = sign(setup.final_displacement)*disp
force_plot = sign(setup.final_displacement)*force

# 找到峰值载荷及其对应位移
peak_idx = argmax(force_plot)
peak_force = force_plot[peak_idx]
peak_disp = disp_plot[peak_idx]
println("峰值载荷: F_max = $(round(peak_force, digits=4)) N @ ū = $(round(peak_disp, digits=4)) mm")

# 图一：载荷-位移曲线
fig_load = Figure(size = (600, 400))
ax_load = Axis(fig_load[1, 1],
    xlabel = L"\bar{u}~\mathrm{[mm]}",
    ylabel = L"F_{\mathrm{reaction}}~\mathrm{[N]}",
    title = "Displacement - Reaction Force Curve",
    #limits = ((0, maximum(disp_plot)), (0, nothing)),
    xgridvisible = true,
    ygridvisible = true,
    xgridcolor = :lightgray,
    ygridcolor = :lightgray,
)
lines!(ax_load, disp_plot, force_plot; linewidth = 2, color = :red, linestyle = :solid) # 红色虚线
save("data/plots/load_displacement_staggered.png", fig_load)

# 图二：能量演变
fig_energy = Figure(size = (600, 400))
ax_energy = Axis(fig_energy[1, 1],
    xlabel = L"\bar{u}~\mathrm{[mm]}",
    ylabel = L"\mathcal{G}_f\ \mathrm{(surface)} [N·mm]",
    title = L"\mathcal{G}_f\ \mathrm{(surface)} Evolution",
    #limits = ((0, maximum(disp_plot)), (0, 1.5)),
    xgridvisible = true,
    ygridvisible = true,
    xgridcolor = :lightgray,
    ygridcolor = :lightgray,
)
#lines!(ax_energy, disp_plot, psi_energy; linewidth = 2, color = :steelblue, label = L"\Psi\ \mathrm{(elastic)}")
lines!(ax_energy, disp_plot, gf_energy; linewidth = 2, color = :darkorange, label = L"\mathcal{G}_f\ \mathrm{(surface)}")
save("data/plots/energy_evolution_staggered.png", fig_energy)

println("载荷-位移曲线已保存至 data/plots/load_displacement_staggered.png。")
println("能量演变曲线已保存至 data/plots/energy_evolution_staggered.png。")

# 保存数据以供后续分析
using JLD2
@save "data/jld2/staggered_results.jld2" disp force psi_energy gf_energy
println("Staggered 仿真数据已保存至 data/jld2/staggered_results.jld2。")