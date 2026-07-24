# RLM–BDF2 实现计划（Mode-I、Miehe 谱分解、原始能量直接闭合）

## 0. 本计划的使用方式

本计划面向接手实现的 Agent。开始编码前必须先读：

1. `docs/RLM-BDF2_Mode_I_Miehe.md`：公式的唯一理论依据；
2. `src/rlm/config.jl`：现有配置、状态、诊断和结果类型；
3. `src/rlm/assembly.jl`：常矩阵、外载、非线性力和能量；
4. `src/rlm/scalar_solver.jl`：现有 BDF1–QM 二次求根；
5. `src/rlm/solver_bdf1.jl`：四分支、事务式提交、能量和输出；
6. `test/runtests.jl`：当前回归基线。

每完成一个阶段就运行该阶段列出的测试。不要一次性改完所有文件后再排错。

本计划只规划代码，不改变理论模型。若代码与本文公式不一致，以
`docs/RLM-BDF2_Mode_I_Miehe.md` 为准；若理论文档自身出现矛盾，停止修改并向用户确认。

---

## 1. 目标与第一版范围

### 1.1 必须实现

在不破坏现有 RLM–BDF1–QM 的前提下，新增固定步长 RLM–BDF2 求解路径：

- 二维平面应变；
- Miehe 拉压谱分解；
- AT2 裂纹能；
- 准静态位移平衡；
- 黏性相场 BDF2 时间离散；
- 非线性力使用二阶显式外推状态；
- 四个仿射线性分支；
- 首选的“原始非线性能量一维直接闭合”；
- BDF1 原始能量直接闭合启动步；
- 两层已接受状态的事务式管理和失败回滚；
- BDF2 两步 \(G\)-能量、代数外功、黏性耗散、数值耗散和闭合残差；
- 单元测试、短算例脚本和必要文档更新。

公开入口应为：

```julia
problem = build_rlm_problem(config)
result = solve_rlm_bdf2(problem)
```

同时公开一个便于单步测试、且绝不修改输入历史状态的函数：

```julia
trial = compute_rlm_bdf2_trial(problem, state_nm1, state_n, t_next)
```

### 1.2 第一版明确不做

- 不实现变步长 BDF2；
- 不在单步失败后局部减半 `dt` 再继续；
- 不施加 \(d_t\ge0\)、\(d^{n+1}\ge d^n\) 或 \(0\le d\le1\)；
- 不裁剪损伤场；
- 不加入惯性；
- 不加入不可逆主动集或历史场；
- 不把 Staggered 改成 BDF2；
- 不删除或改变现有 `solve_rlm_bdf1` 的数值行为；
- 不把 BDF2–QM 作为第一版完成条件；
- 不在标量根失败时接受“使 \(|F|\) 最小但 \(F\ne0\)”的点。

理论文档只完整推导了固定步长 BDF2。虽然其中提到根失败时减小
\(\Delta t\)，第一版的安全解释是：本步完整回滚并停止，用户用更小的统一
`dt` 从 \(t=0\) 重新运行。禁止直接插入一个局部半步，否则后续仍使用固定步长
BDF2 系数会变成未经推导的变步长算法。

### 1.3 闭合路线的优先级

第一版只把以下路线作为必须交付：

1. 第一步：BDF1 四分支 + 原始能量直接闭合；
2. 第二步起：BDF2 四分支 + 原始能量直接闭合。

BDF2–QM 可在上述路线完全通过后作为独立第二阶段增加。不要为了复用当前
BDF1–QM 的二次根代码而把 BDF2 默认实现成 QM。

---

## 2. 当前代码事实与不能踩的坑

### 2.1 可直接复用

- `miehe_response_2d`：Miehe 正负能量和应力；
- `assemble_rlm_nonlinear_forces!`：在给定状态组装
  \(\mathcal N_u,\mathcal N_d\)；
- `rlm_nonlinear_energy`：计算原始 \(\mathcal E_1(U)\)；
- `K_u`、`M_d`、`K_AT2`；
- 位移约束、体力、牵引和物理时间历史；
- BDF1 的位移基准分支处理、反力自由度和输出框架；
- `RLMStepFailure` 事务式失败机制；
- 现有稳定二次求根仅供 BDF1–QM 和将来的 BDF2–QM 使用。

### 2.2 不能原样复用

当前 `RLMProblem.K_d` 和 `factor_d` 在 `build_rlm_problem` 中按

\[
K_d^{\mathrm{BDF1}}=K_{\mathrm{AT2}}+\frac{\eta}{\Delta t}M_d
\]

构造。BDF2 主步必须使用

\[
K_d^{\mathrm{BDF2}}=K_{\mathrm{AT2}}+\frac{3\eta}{2\Delta t}M_d.
\]

因此不能只改 BDF2 的右端后继续调用当前 `factor_d`。另一方面，启动步仍需要
BDF1 矩阵。实现时必须同时提供两套相场矩阵/分解：

- `K_d`、`factor_d`：继续表示 BDF1，保持现有 BDF1 代码兼容；
- `K_d_bdf2`、`factor_d_bdf2`：新增，专供 BDF2 主步。

### 2.3 历史状态

不要把 `RLMState` 改成“内部自带两层历史”的可变对象。保留它表示单个已接受
时刻的状态：

```julia
RLMState(u, d, q)
```

在 `solve_rlm_bdf2` 中显式维护：

```text
state_nm1 = (U^{n-1}, q^{n-1})
state_n   = (U^n,     q^n)
trial     = 尚未接受的 n+1 候选
```

只有整个候选通过根残差、场方程残差和 \(G\)-能量检查后，才执行：

```text
state_nm1 <- state_n
state_n   <- accepted trial
```

任意失败都必须保持两层已接受状态不变。

### 2.4 `RLMTrial` 的风险

当前 `RLMTrial` 字段很多，且使用长位置参数构造。不要继续向这个位置构造器添加
BDF2 字段。新增 `Base.@kwdef struct RLMBDF2Trial`，全部字段按名字构造；现有
`RLMTrial` 留给 BDF1，避免无意改变 BDF1。

---

## 3. 文件级改动总览

| 文件 | 必须改动 |
|---|---|
| `src/rlm/config.jl` | 新增直接求根配置、BDF2 trial、扩展诊断字段；校验配置 |
| `src/rlm/assembly.jl` | 构造 BDF2 相场矩阵/分解；新增仿射方向能量缓存与求值 |
| `src/rlm/scalar_solver.jl` | 保留二次求根；新增一维扫描、Brent、偶重根保守处理 |
| `src/rlm/solver_bdf2.jl` | 新增启动步、BDF2 主步、状态提交、\(G\)-能量和输出 |
| `src/PffSAV.jl` | include 新文件并导出 BDF2 API/配置 |
| `scripts/run_rlm_bdf2.jl` | 新增独立短算例/正式算例入口 |
| `test/runtests.jl` | 新增求根、分支、时间二阶、能量和回滚测试 |
| `STATE.md` | 实现完成并通过测试后才更新“已实现/验证”状态 |
| `PROJECT_STRUCTURE.md` | 实现完成后加入 `solver_bdf2.jl` 和脚本入口 |

不要重命名或删除现有 BDF1 文件和 API。

`src/PffSAV.jl` 中应在 `solver_bdf1.jl` 之后 include `solver_bdf2.jl`。这样
BDF2 可以复用已经定义的通用小工具；反向 include 会造成未定义符号。若发现要复用
的函数其实带有 BDF1 专属语义，应先把它移到 `assembly.jl`、`scalar_solver.jl`
或新的通用辅助文件并保持 BDF1 回归，而不是让 BDF2 隐式依赖错误公式。

---

## 4. 阶段 A：先做无行为变化的基础设施

### A1. 扩展 `RLMProblem`

在 `src/rlm/assembly.jl`：

1. 保留现有 `K_d`、`factor_d` 的含义；
2. 新增：

```julia
K_d_bdf2::SparseMatrixCSC{Float64, Int}
factor_d_bdf2
```

3. 在 `build_rlm_problem` 中计算：

```julia
eta_over_dt = viscosity / dt
K_d = K_AT2 + eta_over_dt * M_d
K_d_bdf2 = K_AT2 + 1.5 * eta_over_dt * M_d
factor_d = cholesky(Symmetric(K_d))
factor_d_bdf2 = cholesky(Symmetric(K_d_bdf2))
```

不要用 `3 * eta_over_dt`；正确系数是 \(3\eta/(2\Delta t)\)。

### A2. 直接求根配置

在 `src/rlm/config.jl` 新增独立配置，例如：

```julia
Base.@kwdef struct RLMDirectRootConfig
    q_min::Float64 = 0.1
    q_max::Float64 = 2.0
    initial_radius::Float64 = 0.05
    radius_growth::Float64 = 2.0
    max_expansions::Int = 8
    scan_points::Int = 65
    max_iterations::Int = 100
    q_abs_tolerance::Float64 = 1.0e-12
    q_rel_tolerance::Float64 = 1.0e-10
    near_zero_factor::Float64 = 10.0
end
```

并作为 `RLMConfig` 的新字段：

```julia
direct_root::RLMDirectRootConfig = RLMDirectRootConfig()
```

约束：

- `q_min > 0`；
- `q_max > q_min`；
- `scan_points >= 3`，最好要求奇数，以便稳定包含中点；
- `initial_radius > 0`；
- `radius_growth > 1`；
- 次数均为非负/正整数；
- 所有浮点数有限；
- 所有容差非负。

标量方程归一化残差阈值继续复用
`RLMToleranceConfig.scalar_residual`，不要再定义第二个含义相同的阈值。

### A3. BDF2 trial 与诊断字段

新增 `RLMBDF2Trial`，至少包含：

```text
四分支：u_a, u_b, d_a, d_b
外推状态：u_bar, d_bar
冻结力：n_u_bar, n_d_bar
新状态：u, d, q
历史能量：phi_nm1, phi_n, phi_next
标量量：c0, c1, q_prediction
求根统计：root_count, function_evaluations, iterations
最终括区：bracket_left, bracket_right, bracket_width
标量残差：scalar_residual_raw, scalar_residual_normalized
场残差：displacement_residual, phase_residual
能量：internal_energy, relaxed_internal_energy, g_energy
外功与耗散：external_work, viscous_dissipation, numerical_dissipation
闭合：g_balance_residual_raw, g_balance_residual_normalized
反力与损伤诊断：reaction_rlm, reaction_phys, phase_increment,
                  phase_relative_increment, healing, min_d, max_d
步骤类型：startup_bdf1 或 bdf2
```

给共享的 `RLMDiagnostic` 增加对应的 BDF2 字段，默认均为 `NaN`、`0` 或空字符串，
确保现有 BDF1 的 `_diagnostic` 不必填它们也能构造。建议新增：

```text
scheme
closure
q_prediction
q_prediction_error
root_count
root_function_evaluations
root_iterations
bracket_left
bracket_right
bracket_width
scalar_residual_raw
displacement_equilibrium_residual
g_energy
g_energy_minus_relaxed_internal
g_balance_residual
g_balance_relative_residual
```

现有 `scalar_residual` 继续表示归一化标量残差。BDF2 中现有 QM 专属字段
`phi_star`、`g_star`、`curvature`、`majorant_margin`、`A/B/C/discriminant`
保持 `NaN`，不要伪造含义。

### A4. 阶段 A 验收

运行完整现有测试：

```powershell
julia --project=. -e "using Pkg; Pkg.test()"
```

此时 BDF1 的所有测试结果必须不变。若失败，先修复兼容性，不要继续写 BDF2。

---

## 5. 阶段 B：实现通用的一维直接求根器

### B1. 不新增外部依赖

当前项目没有 `Roots.jl`。第一版在 `src/rlm/scalar_solver.jl` 内实现最小、可测试的
Brent 求根和扫描逻辑，不修改 `Project.toml`/`Manifest.toml`。

二次求根 `solve_rlm_quadratic` 必须保留，不能改成 Brent。

### B2. 求根结果类型

新增关键字构造的结果类型，例如 `RLMDirectRootResult`：

```text
success
code
message
q
f_q
normalized_residual
candidates
candidate_residuals
root_count
function_evaluations
iterations
bracket_left
bracket_right
bracket_width
```

失败必须有可机器判断的 `code`，至少区分：

- 非有限函数值；
- 搜索区间非法；
- 未找到括区或近零候选；
- Brent 达到最大迭代；
- 找到候选但最终归一化残差不合格。

### B3. 归一化残差

不要只检查 `abs(F(q))`。调用方应同时提供能量尺度函数或在构造残差时返回尺度。
BDF2 主步可使用：

\[
\begin{aligned}
E_{\mathrm{scale}}(q)=\epsilon
&+3|\phi_n(q)|
+3|\alpha-c_1|q^2
+|c_0q|\\
&+4|\mathcal E_1(U^n)|
+|\mathcal E_1(U^{n-1})|
+4\alpha(q^n)^2+\alpha(q^{n-1})^2 .
\end{aligned}
\]

启动步使用同样原则，但系数按启动方程取值。最终接受条件：

\[
\frac{|F(q)|}{E_{\mathrm{scale}}(q)}
\le \texttt{tolerances.scalar_residual}.
\]

分母中的 \(\epsilon\) 使用
`tolerances.scalar_denominator_epsilon`。

### B4. 搜索顺序

求根器按以下固定顺序工作：

1. 计算并验收 `q = 1`；若已满足归一化残差，立即返回；
2. 将预测值 `q_pred` 投影到 `[q_min, q_max]`；
3. 以 `q_pred` 为中心，从 `initial_radius` 开始按 `radius_growth` 扩展局部区间；
4. 检查相邻采样点是否异号；
5. 若局部扩展未找到，使用包含 `q_min`、`q_max`、`1` 和 `q_pred` 的全区间粗扫描；
6. 对每个异号子区间单独调用 Brent；
7. 去重候选根；
8. 每个根重新计算原始方程及归一化残差；
9. 在合格正根中选择离 `q_pred` 最近者；距离相同时选择离 `1` 最近者。

不要只返回扫描发现的第一个根，这会破坏连续物理解支跟踪。

### B5. Brent 的硬性要求

- 只在端点异号或端点本身合格时调用；
- 始终保持根被括在区间内；
- 插值步骤不安全时退回二分；
- 同时使用绝对和相对 `q` 容差；
- 达到最大迭代但残差不合格时返回失败；
- 所有函数值必须检查 `isfinite`；
- 最终结果必须重新走统一验收逻辑。

### B6. 偶重根/近零极小值

只靠异号扫描会漏掉偶重根。第一版至少做到“保守发现、绝不冒充”：

1. 在粗扫描中找离散的 \(|F|\) 局部极小点；
2. 只有当其归一化残差小于
   `near_zero_factor * scalar_residual` 时才作为精化种子；
3. 使用区间保护的 Newton；Newton 导数不可用或步长出界时退回割线/二分式保护步；
4. 精化后仍必须满足正式 `scalar_residual`；
5. 若不满足，返回无合格根，而不是接受局部极小点。

调用方应提供 \(F'(q)\)。BDF2 中：

\[
F_n'(q)=3\phi_n'(q)+6(\alpha-c_1)q-c_0.
\]

启动步中：

\[
F_0'(q)=\phi_0'(q)+2(\alpha-c_1)q-c_0.
\]

### B7. 求根器单元测试

至少测试：

1. 单个简单根；
2. 根正好为 `q=1` 的快速路径；
3. 同一区间多个根，按 `q_pred` 选择；
4. 偶重根，例如 \((q-1.2)^2=0\)；
5. 无根；
6. 函数返回 `NaN`；
7. 极窄括区；
8. 达到最大迭代；
9. 原始残差小但归一化残差/区间条件不满足时不接受。

阶段 B 只测试纯标量函数，不依赖 Ferrite 网格。

---

## 6. 阶段 C：仿射方向能量缓存

直接闭合会多次计算

\[
\phi(q)=\mathcal E_1(U_a+qU_b).
\]

不能在每次 Brent 迭代中重新求解全局线性系统。建议在
`src/rlm/assembly.jl` 新增本步只构造一次的缓存，例如
`RLMAffineEnergyCache`。对每个积分点缓存：

- 求积权重；
- \(d_a\) 和 \(d_b\) 的积分点值；
- \(\epsilon(u_a)\) 和 \(\epsilon(u_b)\)；
- 计算 Miehe 响应所需的材料和容差由 `problem` 提供。

然后实现：

```julia
rlm_affine_nonlinear_energy(problem, cache, q)
rlm_affine_nonlinear_energy_derivative(problem, cache, q)
```

求值公式：

\[
d(q)=d_a+qd_b,\qquad
\epsilon(q)=\epsilon_a+q\epsilon_b,
\]

\[
\phi(q)=\int_\Omega
[g_\kappa(d(q))-1]\psi_0^+(\epsilon(q))\,dx,
\]

\[
\phi'(q)=\int_\Omega
\left[
g_\kappa'(d(q))d_b\psi_0^+(\epsilon(q))
+(g_\kappa(d(q))-1)\sigma_0^+(\epsilon(q)):\epsilon_b
\right]dx.
\]

要求：

- 使用与 `assemble_rlm_nonlinear_forces!`、`rlm_nonlinear_energy` 完全相同的
  `cellvalues`、积分点和权重；
- 不缓存依赖 `q` 的谱分解结果；
- 用普通 `rlm_nonlinear_energy(problem, u_a + q*u_b, d_a + q*d_b)` 对多个
  `q` 值交叉验证；
- 用中心差分验证 \(\phi'(q)\)，但中心差分只用于测试，不用于正式导数；
- Miehe 主应变过零附近放宽导数差分测试或避开该点，不要错误要求全局 \(C^2\)。

---

## 7. 阶段 D：BDF1 原始能量启动步

BDF2 需要 \(U^0,U^1\)。不要调用现有 `compute_rlm_bdf1_trial`，因为它使用 QM
闭合。新增 BDF2 文件内的私有函数，例如：

```julia
_compute_rlm_bdf1_direct_startup_trial(problem, state_0, t_1)
```

### D1. 冻结非线性力

在 \(U^0\) 组装：

\[
N_u^0=\mathcal N_u(U^0),\qquad
N_d^0=\mathcal N_d(U^0).
\]

### D2. 四分支

位移分支与现有 BDF1 相同：

\[
K_u u_a^1=f_{\mathrm{ext}}(t_1),\qquad
K_u u_b^1=-N_u^0,
\]

其中 \(u_a^1\) 施加 \(t_1\) 的非齐次 Dirichlet 条件，\(u_b^1\) 是齐次条件。

相场分支必须使用 BDF1 的 `factor_d`：

\[
K_d^{\mathrm{BDF1}}d_a^1
=\frac{\eta}{\Delta t}M_dd^0,\qquad
K_d^{\mathrm{BDF1}}d_b^1=-N_d^0.
\]

### D3. 启动标量方程

\[
c_0^0=N_u^0\cdot(u_a^1-u^0)+N_d^0\cdot(d_a^1-d^0),
\]

\[
c_1^0=N_u^0\cdot u_b^1+N_d^0\cdot d_b^1,
\]

\[
F_0(q)=\phi_0(q)+(\alpha-c_1^0)q^2-c_0^0q
-\mathcal E_1(U^0)-\alpha(q^0)^2.
\]

预测值使用 `q_pred = q^0 = 1`。用阶段 B 的直接求根器。

### D4. 启动步验证

接受前检查：

\[
c_1^0\approx
-u_b^{1T}K_uu_b^1-d_b^{1T}K_d^{\mathrm{BDF1}}d_b^1\le0,
\]

自由位移自由度上的平衡：

\[
K_uu^1+q^1N_u^0-f_{\mathrm{ext}}(t_1)=0,
\]

相场平衡：

\[
\frac{\eta}{\Delta t}M_d(d^1-d^0)
+K_{\mathrm{AT2}}d^1+q^1N_d^0=0.
\]

启动失败时最终 `RLMResult.state` 必须仍是 \(U^0\)，诊断记录为未接受。

---

## 8. 阶段 E：BDF2 四分支与直接闭合

新增 `src/rlm/solver_bdf2.jl`。

### E1. 单步函数签名与纯事务要求

```julia
function compute_rlm_bdf2_trial(
    problem::RLMProblem,
    state_nm1::RLMState,
    state_n::RLMState,
    t_next::Real,
)
```

函数内不得修改 `state_nm1.u/d`、`state_n.u/d` 或两者的 `q`。所有数组均新建，
失败通过 `RLMStepFailure` 抛出。

检查：

```julia
t_next ≈ t_n + dt
```

时间索引由外层明确传入/计算，不能用“当前调用次数”猜测时间。

### E2. 外推与冻结力

\[
\bar u^{n+1}=2u^n-u^{n-1},\qquad
\bar d^{n+1}=2d^n-d^{n-1}.
\]

只组装一次：

\[
\bar N_u^{n+1}=\mathcal N_u(\bar U^{n+1}),\qquad
\bar N_d^{n+1}=\mathcal N_d(\bar U^{n+1}).
\]

禁止在 Brent 的每个 `q` 上重新组装场方程用的 \(\bar N\)。标量迭代中的
\(\phi(q)\) 是候选新状态的原始能量，与冻结力的外推状态不是同一个对象。

### E3. 四分支

先把约束更新到 \(t_{n+1}\)，把外力更新到 \(t_{n+1}\)。

位移：

\[
K_uu_a^{n+1}=f_{\mathrm{ext}}(t_{n+1}),
\qquad
K_uu_b^{n+1}=-\bar N_u^{n+1}.
\]

相场必须使用 `factor_d_bdf2`：

\[
K_d^{\mathrm{BDF2}}d_a^{n+1}
=\frac{\eta}{2\Delta t}M_d(4d^n-d^{n-1}),
\]

\[
K_d^{\mathrm{BDF2}}d_b^{n+1}=-\bar N_d^{n+1}.
\]

最终场：

\[
u^{n+1}=u_a^{n+1}+q^{n+1}u_b^{n+1},\qquad
d^{n+1}=d_a^{n+1}+q^{n+1}d_b^{n+1}.
\]

重构后再次对 `u_next` 应用当前非齐次位移约束，以消除舍入误差；`u_b` 的受约束
自由度必须为零。

### E4. 公共标量

\[
\begin{aligned}
c_0={}&
\bar N_u\cdot(3u_a-4u^n+u^{n-1})\\
&+\bar N_d\cdot(3d_a-4d^n+d^{n-1}),
\end{aligned}
\]

\[
c_1=\bar N_u\cdot u_b+\bar N_d\cdot d_b.
\]

必须检查分支恒等式：

\[
c_1\approx
-u_b^TK_uu_b-d_b^TK_d^{\mathrm{BDF2}}d_b\le0.
\]

沿用现有 `branch_identity`、`c1_abs`、`c1_rel` 容差。若显著不符，抛出
`RLMStepFailure(:c1_branch_identity, ...)` 或 `:positive_c1`。

### E5. BDF2 原始能量残差

先直接求积并保存：

\[
\phi_n=\mathcal E_1(U^n),\qquad
\phi_{n-1}=\mathcal E_1(U^{n-1}).
\]

建立：

\[
\begin{aligned}
F_n(q)={}&3\phi(q)+3(\alpha-c_1)q^2-c_0q\\
&-4\phi_n+\phi_{n-1}
-4\alpha(q^n)^2+\alpha(q^{n-1})^2.
\end{aligned}
\]

预测值：

\[
q_{\mathrm{pred}}=2q^n-q^{n-1}.
\]

调用阶段 B 的求根器。找到多个根时按预测值选择，不按根的绝对大小选择。接受后
保存精确的：

\[
\phi_{n+1}=\phi(q^{n+1}).
\]

### E6. 场方程残差

定义：

\[
\Delta_Bu=3u^{n+1}-4u^n+u^{n-1},
\]

\[
\Delta_Bd=3d^{n+1}-4d^n+d^{n-1}.
\]

位移残差：

\[
r_u=K_uu^{n+1}+q^{n+1}\bar N_u-f_{\mathrm{ext}}(t_{n+1}).
\]

只在 `free_dofs(problem.ch_u)` 上检查平衡；受约束自由度用于反力。

相场残差：

\[
r_d=\frac{\eta}{2\Delta t}M_d\Delta_Bd
+K_{\mathrm{AT2}}d^{n+1}+q^{n+1}\bar N_d.
\]

两者都要按各项范数之和归一化。不要只记录、不验收。若超出配置容差，候选不提交。
可新增专门的 `linear_equilibrium` 容差，或明确复用现有
`branch_identity`；优先新增语义清楚的容差字段。

---

## 9. 阶段 F：BDF2 外功与 \(G\)-能量

这部分是实现正确性的核心，不能用 BDF1 的累计能量公式代替。

### F1. 线性能量范数辅助函数

实现私有纯函数：

\[
\|U\|_{\mathcal L}^2=u^TK_uu+d^TK_{\mathrm{AT2}}d.
\]

注意相场部分使用 `K_AT2`，不是含黏性质量项的 `K_d_bdf2`。

### F2. 两步 \(G\)-能量

对一对已接受状态 `(current, previous)` 计算：

\[
\begin{aligned}
\mathscr G^n={}&
\frac14\left[
\|U^n\|_{\mathcal L}^2+
\|2U^n-U^{n-1}\|_{\mathcal L}^2
\right]\\
&+\frac32\mathcal E_1(U^n)
-\frac12\mathcal E_1(U^{n-1})\\
&+\alpha\left[
\frac32(q^n)^2-\frac12(q^{n-1})^2-1
\right].
\end{aligned}
\]

主步中：

```text
g_old = G(state_n, state_nm1)
g_new = G(state_next, state_n)
```

`internal_energy` 和 `relaxed_internal_energy` 仍按单时刻物理含义输出，不要把它们
重命名为 `g_energy`。同时输出：

\[
\mathscr G^n-\widetilde{\mathcal H}(U^n,q^n),
\]

用于时间加密检查。

### F3. 一致 Dirichlet 提升

新增辅助函数生成三个时刻的纯边界位移向量：

```text
uD_nm1 at t_{n-1}
uD_n   at t_n
uD_np1 at t_{n+1}
```

每个向量先置零，再 `update!` 对应时刻的约束并 `apply!`。构造：

\[
z_B=3u_D^{n+1}-4u_D^n+u_D^{n-1}.
\]

辅助函数结束前必须把 `problem.ch_u` 恢复到 \(t_{n+1}\)，防止后续输出或约束使用
了旧时间。

### F4. 代数 BDF2 外功

使用冻结的 \(\bar N_u\)，定义完整机械残差（包括受约束自由度）：

\[
r_{\mathrm{mech}}=K_uu^{n+1}+q^{n+1}\bar N_u-f_{\mathrm{ext}}(t_{n+1}).
\]

计算：

\[
W_{\mathrm{ext},B}^{n+1}
=\frac12\left[
f_{\mathrm{ext}}(t_{n+1})^T\Delta_Bu
+r_{\mathrm{mech}}^Tz_B
\right].
\]

不要使用 `reaction_phys` 代替这里的 RLM 一致反力；物理反力仍可单独诊断。

### F5. 两个耗散项

二阶差分：

\[
\Delta^2u=u^{n+1}-2u^n+u^{n-1},\qquad
\Delta^2d=d^{n+1}-2d^n+d^{n-1}.
\]

数值耗散：

\[
D_{\mathrm{num}}=
\frac14\left[
(\Delta^2u)^TK_u\Delta^2u+
(\Delta^2d)^TK_{\mathrm{AT2}}\Delta^2d
\right].
\]

黏性耗散：

\[
D_{\mathrm{vis}}=
\frac{\eta}{4\Delta t}
(\Delta_Bd)^TM_d\Delta_Bd.
\]

不要照搬 BDF1 中的 `0.5` 系数和 `eta/dt * ||d^{n+1}-d^n||²`。

### F6. 精确 \(G\)-能量闭合

原始能量直接闭合没有 QM 缺口：

\[
r_G=
\mathscr G^{n+1}-\mathscr G^n
-W_{\mathrm{ext},B}^{n+1}
+D_{\mathrm{num}}+D_{\mathrm{vis}}.
\]

归一化尺度至少包含上述五项绝对值之和与 epsilon。要求：

```text
abs(r_G) <= energy_balance_abs + energy_balance_rel * scale
```

超限必须使本步失败并回滚，不能只记 warning。

### F7. 反力

- `reaction_rlm`：用 `r_mech` 在加载分量自由度求和；
- `reaction_phys`：在最终 \(U^{n+1}\) 重新组装真实
  \(\mathcal N_u(U^{n+1})\)，然后用
  `K_u*u_next + n_u_phys - f_ext` 求和。

两者含义不同，均保留。

---

## 10. 阶段 G：外层 `solve_rlm_bdf2`

### G1. 初始化

与 BDF1 一样构造 \(U^0\)：

- `u = 0` 后施加 \(t=0\) Dirichlet；
- `d = initial_damage`；
- `q = 1`；
- 写入 initial diagnostic/VTK（若配置要求）。

要求 BDF2 至少有两个时间区间：

```text
nsteps = final_time / dt
nsteps >= 2
```

这个约束只在 `solve_rlm_bdf2` 检查，不要让 BDF1 的单步算例失效。

### G2. 启动

调用阶段 D 的直接闭合启动函数得到 \(U^1\)。只有成功后才设置：

```text
state_nm1 = state_0
state_n = state_1
```

启动诊断明确写：

```text
scheme = "BDF1-startup"
closure = "direct-original-energy"
```

不要为启动步伪造 BDF2 \(G\)-平衡；可以记录其 BDF1 直接闭合能量关系，或将
BDF2 专属 `g_balance_residual` 留为 `NaN`。从 \(n=1\to2\) 起才检查 BDF2
\(G\)-平衡。

### G3. 主循环

```text
for step in 2:nsteps
    t_next = step * dt（末步用 final_time 消除浮点漂移）
    trial = compute_rlm_bdf2_trial(problem, state_nm1, state_n, t_next)
    trial 全部验收后：
        state_nm1 = state_n
        state_n = RLMState(trial.u, trial.d, trial.q)
        写诊断/VTK/CSV
end
```

每个物理时间步只推进一次状态。Brent 的多次函数求值不是多次物理时间推进。

### G4. 失败

捕获 `RLMStepFailure` 后：

1. 写一条 `accepted=false` 的失败诊断；
2. 状态仍返回最后两层中较新的已接受 `state_n`；
3. `success=false`、`completed=false`；
4. 刷新 CSV；
5. 消息明确包含失败码、步号和“candidate rolled back”；
6. 不写失败候选的 VTK；
7. 不更新累计量。

若失败发生在启动步，返回 \(U^0\)。

### G5. 输出

新增 `scripts/run_rlm_bdf2.jl`，不要把 BDF2 开关塞进
`scripts/run_rlm_bdf1.jl`。BDF2 默认输出目录不得与 BDF1 相同，例如：

```text
data/sims/rlm_bdf2
data/jld2/rlm_bdf2_results.jld2
data/plots/rlm_bdf2_*.png
```

CSV 至少包含：

- `step/time/dt/scheme/closure/accepted/status`；
- 位移、两种反力、损伤范围；
- `q`、`q-1`、预测值和预测误差；
- 根数、函数调用次数、括区宽度、标量归一化残差；
- 单时刻物理能量、松弛能量、\(G\)-能量；
- 单步 BDF2 外功、两种耗散、\(G\)-闭合残差；
- 位移与相场方程残差。

当前 `write_rlm_diagnostics` 会遍历 `RLMDiagnostic` 的全部字段，扩展类型后可自然
包含新字段；但 `write_rlm_time_history` 使用固定字段元组，不会自动加入 BDF2
字段。实现者必须二选一：

1. 扩展这个固定字段元组，并确认 BDF1 CSV 仍能写出；或
2. 新增 `write_rlm_bdf2_time_history`，由 BDF2 外层单独调用。

不要以“完整 diagnostics.csv 里已经有”为由遗漏精简 time-history 文件中的核心
BDF2 诊断。

启动步的 `external_work` 是 BDF1 量，后续主步的 `external_work` 是
\(W_{\mathrm{ext},B}\)，二者不能不加区分地累计后再与单一能量公式比较。建议：

- 保留逐步 `external_work`，由 `scheme` 说明含义；
- 新增 `cumulative_bdf2_external_work`，从第一个 BDF2 主步开始累加；
- BDF2 的累计黏性/数值耗散也从第一个 BDF2 主步开始；
- 启动步若需要累计量，使用明确的 `startup_*` 字段或只记录单步值。

进度条文案必须是 `RLM-BDF2-direct`，不能仍显示 `RLM-BDF1-QM`。

---

## 11. 阶段 H：测试矩阵

所有数值测试使用小网格、关闭 CSV/VTK/进度输出。

### H1. BDF1 回归

完整现有测试必须继续通过，重点确认：

- `solve_rlm_bdf1` 的步数、状态和诊断不变；
- BDF1 仍使用 `factor_d`，没有误用 BDF2 分解；
- BDF1–QM 的二次求根和 majorant 检查不变。

### H2. 两套相场矩阵

直接检查：

\[
K_d-K_{\mathrm{AT2}}=\frac{\eta}{\Delta t}M_d,
\]

\[
K_d^{\mathrm{BDF2}}-K_{\mathrm{AT2}}
=\frac{3\eta}{2\Delta t}M_d.
\]

### H3. 仿射能量缓存

对多个正 `q`：

```julia
cached_phi ≈ rlm_nonlinear_energy(problem, u_a + q*u_b, d_a + q*d_b)
```

并验证缓存导数与方向有限差分一致。

### H4. BDF1 启动步

验证输入 `state_0` 未被修改，并检查：

- 位移自由自由度平衡；
- 相场 BDF1 平衡；
- `c1` 分支恒等式；
- 标量直接闭合残差；
- `q` 位于区间内。

### H5. BDF2 单步四分支

手工准备已接受的 `state_0/state_1`，调用
`compute_rlm_bdf2_trial`，检查：

- 两个输入状态逐元素未变；
- \(\bar U=2U^1-U^0\)；
- 相场基准右端确实是
  \(\eta/(2dt)M_d(4d^1-d^0)\)；
- 位移和相场残差；
- `c1` 恒等式和非正性；
- 标量原始能量残差；
- 根选择接近 \(2q^1-q^0\)。

### H6. 精确 \(G\)-能量闭合

在非零位移历史以及可选体力/牵引下，至少走到 `step=2`，检查：

```text
abs(g_balance_residual) <= configured tolerance
```

同时独立重算 \(z_B\)、外功、两种耗散和 \(G\)-能量，避免测试只是重复调用被测辅助函数。

### H7. 均匀黏性衰减的二阶时间收敛

复用当前 BDF1 测试的无载、均匀初始损伤问题。解析解：

\[
d(t)=d_0\exp\left(-\frac{G_c}{\eta\ell}t\right).
\]

用相同最终时间计算 `dt`、`dt/2`、`dt/4` 的误差，估计：

\[
p=\log_2(e_{\Delta t}/e_{\Delta t/2}).
\]

要求两次估计接近 2；建议使用不过度脆弱的区间，例如 `1.7 < p < 2.3`，并确认
网格误差在该均匀问题中不会主导。不能只断言“细步长误差更小”。

同时保留现有 BDF1 一阶测试，形成明确对照。

### H8. 载荷历史

使用非零且彼此独立的位移、体力、牵引历史，验证 BDF2 主步的外力取
\(t_{n+1}\)，Dirichlet 提升使用三个正确时刻。

在分段线性历史的折点附近不做严格二阶阶数断言；时间阶测试使用光滑/常系数问题。

### H9. 失败回滚

至少覆盖：

1. 启动步根区间排除物理解：返回 \(U^0\)；
2. BDF2 单步无合格根：`state_nm1/state_n` 均未修改；
3. 人为收紧 \(G\)-能量容差触发失败：累计外功/耗散不增加；
4. 失败诊断 `accepted=false`，无失败候选 VTK。

### H10. 完整测试

```powershell
julia --project=. -e "using Pkg; Pkg.test()"
```

再运行一个只有数步、关闭大规模输出的 BDF2 smoke test。不要在自动测试中使用
正式细网格或十万步算例。

---

## 12. 建议实施顺序与每阶段停止点

1. **基础设施**：双相场分解 + 配置 + 类型；跑全部旧测试。
2. **纯标量求根器**：不接 FEM；跑求根单测。
3. **仿射能量缓存**：与原能量逐点对照；跑缓存/导数测试。
4. **BDF1 直接启动**：只做一步；检查三类残差。
5. **BDF2 单步**：先不写外层循环；检查四分支和根。
6. **\(G\)-能量**：独立重算并通过闭合测试。
7. **外层求解与回滚**：短程 2–4 步。
8. **二阶收敛测试**。
9. **脚本和文档**。
10. **完整回归与短算例**。

若某阶段未通过，不要继续到下一阶段。

---

## 13. 可选第二阶段：BDF2–QM

只有直接闭合完全通过后才考虑。实现时新增显式配置，例如：

```text
closure = :direct_original_energy   # 默认
closure = :qm                       # 可选
```

BDF2–QM 必须使用理论文档第 5 节的系数：

\[
A=3(\alpha+S/2-c_1),\quad
B=3(g^*-S)-c_0,
\]

\[
C=3(\phi^*-g^*+S/2)-4\phi_n+\phi_{n-1}
-4\alpha(q^n)^2+\alpha(q^{n-1})^2.
\]

不能直接复用 BDF1 的 \(A,B,C\)。根选择使用
\(2q^n-q^{n-1}\)，最终必须检查 QM 上界。\(G\)-平衡中额外加入：

\[
\frac32\varepsilon_{\mathrm{QM}},
\qquad
\varepsilon_{\mathrm{QM}}=\Pi(q^{n+1})-\mathcal E_1(U^{n+1})\ge0.
\]

该可选阶段不应改变直接闭合的默认行为和测试。

---

## 14. 完成

提醒用户人工测试是否满足以下条件：

- BDF1 启动使用原始能量直接闭合，不是现有 QM 闭合；
- BDF2 主步使用正确的 \(3\eta/(2dt)\) 矩阵和历史右端；
- 非线性力只在 \(2U^n-U^{n-1}\) 冻结一次；
- 标量方程使用真实 \(\mathcal E_1(U_a+qU_b)\)；
- 根满足区间、有限性、连续分支和归一化残差检查；
- 两层历史仅在整步成功后提交；
- \(G\)-能量闭合残差达到配置容差；
- 失败回滚测试通过；
- `scripts/run_rlm_bdf2.jl` 的算例可运行；
- 输出目录不会覆盖 BDF1 结果；
- `PROJECT_STRUCTURE.md` 和 `STATE.md` 与最终代码一致。

交接时必须报告：

1. 修改了哪些文件；
2. 采用的是直接闭合还是另加了 QM；
3. 运行了哪些测试及结果；
4. 最大标量残差、最大场残差和最大 \(G\)-平衡残差；
5. 是否仍存在根括区、谱分解低正则性或大时间步失效问题。
