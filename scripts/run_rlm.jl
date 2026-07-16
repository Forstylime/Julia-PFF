using PffSAV
using Ferrite
using FerriteGmsh

# 快速示例：从项目根目录运行
#   julia --project=. scripts/run_rlm.jl

mat = MaterialParams(;
    dim=2,
    E=25840.0,
    ν=0.18,
    gc=2.7,
    l=0.015,
    k=1e-8
)

rlm_params = RLMParams(;
    Δt=1e6, # 使其退化为类似Staggered的无粘性行为
    M=1.0,  
    ϵ=1e-5, 
    θ=1e-10   
)

setup = setup_tension(
    msh_file="data/mesh/l_shape.msh",
    final_displacement=-0.8, # 施加位移
    fixed_face="top", # 固定边界名称
    dir=2, # 位移施加方向，x方向为1，y方向为2
)


disp, force = solve_rlm_amor(setup, mat, rlm_params;
    n_steps=500,
    tol=1e-12,
)

# 可视化载荷-位移曲线
mkpath("data/plots")
mkpath("data/jld2")
using CairoMakie

disp_plot = sign(setup.final_displacement) * disp
force_plot = sign(setup.final_displacement) * force

# 找到峰值载荷及其对应位移
peak_idx = argmax(force_plot)
peak_force = force_plot[peak_idx]
peak_disp = disp_plot[peak_idx]
println("峰值载荷: F_max = $(round(peak_force, digits=4)) N @ ū = $(round(peak_disp, digits=4)) mm")

# 图一：载荷-位移曲线
fig_load = Figure(size=(600, 400))
ax_load = Axis(fig_load[1, 1],
    xlabel=L"\bar{u}~\mathrm{[mm]}",
    ylabel=L"F_{\mathrm{reaction}}~\mathrm{[N]}",
    title="Displacement - Reaction Force Curve (RLM Amor)",
    xgridvisible=true,
    ygridvisible=true,
    xgridcolor=:lightgray,
    ygridcolor=:lightgray,
)
lines!(ax_load, disp_plot, force_plot; linewidth=2, color=:blue, linestyle=:solid)
save("data/plots/load_displacement_rlm.png", fig_load)

println("载荷-位移曲线已保存至 data/plots/load_displacement_rlm.png。")

# 保存数据以供后续分析
using JLD2
@save "data/jld2/rlm_results.jld2" disp force
println("RLM 仿真数据已保存至 data/jld2/rlm_results.jld2。")
