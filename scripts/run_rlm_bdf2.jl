using PffSAV
using Ferrite
using FerriteGmsh
using CairoMakie
using JLD2

mat = MaterialParams(;
    dim=2,
    E=25840.0,
    ν=0.18,
    gc=2.7,
    l=0.015,
    k=1e-8
)

# 注意：对于 BDF2 格式，必须使用合理的物理时间步长 Δt！
# 论文中推荐 M=1.0，Δt 需要使得 1 / Δt 能体现出适当的粘性，或者足够小以捕捉裂纹演化。
# 比如 500 个增量步，如果是准静态加载，可以设 T_final = 1.0, Δt = 1.0 / 500 = 0.002
rlm_params = RLMParams(;
    Δt=0.002, 
    M=500.0,  
    ϵ=1e-5, 
    θ=1e-10   
)

setup = setup_tension(
    msh_file="data/mesh/l_shape.msh",
    final_displacement=0.8,
    fixed_face="top",
    dir=2,
)

println("网格准备完成。位移自由度: ", ndofs(setup.dh_u), " 相场自由度: ", ndofs(setup.dh_d))
println("开始 BDF2-RLM 仿真...")

disp, force = solve_rlm_bdf2(setup, mat, rlm_params; n_steps=500)

disp_plot = sign(setup.final_displacement) * disp
force_plot = sign(setup.final_displacement) * force

peak_idx = argmax(abs.(force_plot))
peak_force = force_plot[peak_idx]
peak_disp = disp_plot[peak_idx]
println("峰值载荷: F_max = $(round(peak_force, digits=4)) N @ ū = $(round(peak_disp, digits=4)) mm")

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Displacement (mm)", ylabel = "Force (N)", title = "Load-Displacement Curve (RLM BDF2)")
lines!(ax, disp_plot, abs.(force_plot), color = :blue, linewidth = 2)

mkpath("data/plots")
save("data/plots/load_displacement_rlm_bdf2.png", fig)
println("载荷-位移曲线已保存至 data/plots/load_displacement_rlm_bdf2.png。")

mkpath("data/jld2")
@save "data/jld2/rlm_bdf2_results.jld2" disp force
println("RLM BDF2 仿真数据已保存至 data/jld2/rlm_bdf2_results.jld2。")
