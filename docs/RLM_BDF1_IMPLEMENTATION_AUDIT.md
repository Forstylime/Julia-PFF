# Miehe–RLM-PC–BDF1 实现审核与设计基线

本文是 `RLM_Mode_I_Miehe.md` 到本项目 Ferrite.jl 离散实现的对应关系审核。第一阶段只实现文档第 9–16、19–21 节的 BDF1/RLM-PC，不实现 CN、RLM-Q、裂纹不可逆约束或事后区间投影。本文不替换原理论文档，也不改变其中的符号、能量或标量二次方程。

## 1. 弱式到代数系统的逐项审核

取独立的有限元空间

\[
V_h=\operatorname{span}\{\boldsymbol\phi_i\}_{i=1}^{N_u},
\qquad
W_h=\operatorname{span}\{\chi_i\}_{i=1}^{N_d}.
\]

位移和相场使用两个独立 `DofHandler`；位移每个 Q1 节点两个自由度，相场每个 Q1 节点一个自由度。两个处理器共享同一网格和单元顺序，但全局编号互不混合。

定义

\[
(K_u)_{ij}=a_u(\boldsymbol\phi_j,\boldsymbol\phi_i),
\quad
(M_d)_{ij}=(\chi_j,\chi_i),
\quad
(K_{AT2})_{ij}=a_d(\chi_j,\chi_i),
\]

其中

\[
K_d=\frac{1}{M\Delta t}M_d+K_{AT2}.
\]

外载向量和冻结非线性力向量定义为

\[
(f_{ext})_i=\ell_{ext}(\boldsymbol\phi_i),
\quad
(n_u^n)_i=\mathcal N_u^n(\boldsymbol\phi_i),
\quad
(n_d^n)_i=\mathcal N_d^n(\chi_i).
\]

因此四个仿射分支的未消元代数式严格对应为

\[
K_u u_a=f_{ext},\qquad K_u u_b=-n_u^n,
\]

\[
K_d d_a=\frac{1}{M\Delta t}M_d d^n,
\qquad K_d d_b=-n_d^n.
\]

预测态、预测能量和标量内积为

\[
u_*=u_a+u_b,\quad d_*=d_a+d_b,
\quad P^{n+1}=\mathcal E_1(u_*,d_*),
\]

\[
c_0^n=(n_u^n)^T(u_a-u^n)+(n_d^n)^T(d_a-d^n),
\]

\[
c_1^n=(n_u^n)^Tu_b+(n_d^n)^Td_b.
\]

代回文档第 12 节标量关系后得到的系数没有改写：

\[
A_n=\alpha-c_1^n,
\qquad B_n=-c_0^n,
\qquad
C_n=P^{n+1}-P^n-\alpha(q^n)^2.
\]

这与文档中的二次方程

\[
A_n(q^{n+1})^2+B_nq^{n+1}+C_n=0
\]

完全一致。响应分支还给出

\[
c_1^n=-u_b^TK_uu_b-d_b^TK_dd_b\leq0,
\]

这里第二项已经同时包含 \((M\Delta t)^{-1}M_d\) 和 \(K_{AT2}\)。程序会同时用向量泛函定义和该负二次型恒等式检查装配符号。

原始能量与 RLM 代理能量分别按

\[
\mathcal E(u,d)=\tfrac12u^TK_uu+\tfrac12d^TK_{AT2}d-f_{ext}^Tu+\mathcal E_1(u,d),
\]

\[
\widetilde{\mathcal E}_{PC}^n=
\tfrac12(u^n)^TK_uu^n+
\tfrac12(d^n)^TK_{AT2}d^n-
f_{ext}^Tu^n+P^n+\alpha[(q^n)^2-1]
\]

计算。`Pⁿ` 保存的是上一步实际预测态的 \(\mathcal E_1\)，不会被 \(\mathcal E_1(U^n)\) 替换。

## 2. 二维平面应变与谱分解约定

使用三维各向同性材料的 Lamé 常数

\[
\lambda=\frac{E\nu}{(1+\nu)(1-2\nu)},
\qquad
\mu=\frac{E}{2(1+\nu)},
\]

并取平面应变 \(\epsilon_{zz}=\epsilon_{xz}=\epsilon_{yz}=0\)。零的第三主应变不改变文档给出的二维面内 \(\psi_0^\pm\) 或面内应力表达式，因此在正交点对面内 2×2 对称应变做解析谱分解即可；不使用平面应力等效 \(\lambda\)。

数值实现不保存特征向量历史。对 2×2 对称矩阵用迹、偏量和 `hypot` 构造两个主应变及谱张量，避免近重根时显式特征向量方向抖动。主应变在

\[
|\epsilon_a|\le \epsilon_{zero,abs}
+\epsilon_{zero,rel}\max_b|\epsilon_b|
\]

时视为正分支的数值零。为保持有限元分裂
\(\psi_0=\psi_0^++\psi_0^-\) 与常矩阵 \(K_u\) 一致，该容差带内的微小量保留在互补分支中，而不是从总应变中删除。默认值显式放在容差配置中；这是对文档 Macaulay 括号在浮点数上的分支约定，不改变其解析定义。近重根另用独立的绝对/相对阈值。

RLM 使用文档原式

\[
g_\kappa(d)=(1-\kappa)(1-d)^2+\kappa,
\qquad
g_\kappa'(d)=-2(1-\kappa)(1-d).
\]

现有 Staggered 模块中的 `(1-d)^2 + k` 不是这个归一化函数，不能复用；RLM 模块将独立实现文档定义，也不会静默改变旧模块。

## 3. Dirichlet 消元方案

设自由自由度为 `F`、Dirichlet 自由度为 `C`，当前固定载荷阶段的边界值为 `g`。基准位移分支采用

\[
(K_u)_{FF}(u_a)_F=(f_{ext})_F-(K_u)_{FC}g,
\qquad (u_a)_C=g.
\]

响应分支必须采用齐次边界：

\[
(K_u)_{FF}(u_b)_F=-(n_u^n)_F,
\qquad (u_b)_C=0.
\]

实现使用 Ferrite 的 `RHSData` 从同一份未约束 `K_u` 产生两类右端，然后只对约束后的常矩阵分解一次。绝不先给 `u_b` 施加非零边界再相减。由于相场为齐次自然 Neumann 条件，没有相场 Dirichlet 消元；\(G_c/\ell>0\) 和质量项使 `K_d` 正定。

## 4. 模块结构

- `src/rlm/config.jl`：材料、网格、加载、时间、容差、输出显式配置，以及状态、诊断、失败结果类型。
- `src/rlm/miehe2d.jl`：稳健二维谱分解、\(\psi_0^\pm\)、\(\sigma_0^\pm\)、\(g_\kappa\) 与导数。
- `src/rlm/assembly.jl`：`K_u`、`M_d`、`K_AT2`、外载、`n_u`、`n_d`、能量与 FE 范数组装。
- `src/rlm/scalar_solver.jl`：不改写系数的稳定二次求根、候选根筛选和归一化残差。
- `src/rlm/solver_bdf1.jl`：四分支、预测态、Load–Relax、事务式接受/回滚、CSV/VTK 输出。
- `scripts/run_rlm_bdf1.jl`：根目录 `l_shape.msh` 的可重复入口。
- `test/`：本构、求根、装配恒等式、一步能量律、回滚和微型 Load–Relax 测试。

## 5. 仍需选择或原文未唯一规定的实现细节

下列内容不是理论文档唯一确定的；第一阶段采用右列的显式选择，并全部暴露在配置中。

| 细节 | 第一阶段选择 |
|---|---|
| L 形网格的单元类型/阶次 | 延续当前项目：二维四节点 Q1 四边形；导入后若不是该类型，立即报错，不猜测转换。 |
| 默认力学边界条件 | 延续现有 L 形脚本语义：`top` 两分量固定，`right` 的第 2 分量递增位移；边界名和方向可配置。 |
| 固定/加载边界共享角点 | Ferrite 对重复自由度采用后加入条件；配置 `overlap_policy=:loaded` 显式规定默认由加载值优先，也可取 `:fixed`。不依赖隐含的添加顺序。 |
| 体力与 Neumann 力 | 入口显式取零，因此 \(f_{ext}=0\)；代码中的能量和代数定义仍保留 `f_ext`，不把 Dirichlet 反力错误计入 \(\ell_{ext}\)。 |
| 初始损伤/预制裂纹 | 默认 `d⁰=0`。几何缺口由网格表达；第一阶段不增加 `d=1` 裂纹 Dirichlet 条件。 |
| 新载荷级的初始位移 | 复制上一个载荷级最终位移并覆盖新 Dirichlet 值，使其属于新仿射空间；然后按文档重置 `q=1`、`P=E₁(U)`。 |
| 迁移率与旧参数 `η` 的关系 | 新配置直接给出 \(M>0\)，不静默假设 \(M=1/\eta\)。 |
| \(\alpha\) 的尺度 | 直接以有能量量纲的固定 `alpha` 配置；每个固定载荷阶段内不改变。入口中的数值只是可编辑算例值，不宣称普适。 |
| 正交阶次 | Q1 四边形采用 2×2 Gauss；可配置但不得低于 2。 |
| 主应变过零/近重根容差 | 绝对加相对尺度阈值，见第 2 节；数值写入诊断配置。 |
| “显著负判别式”的尺度 | `tol_D = disc_abs + disc_rel*(B²+|4AC|)`；`D∈[-tol_D,0)` 才置零，更负则失败。 |
| \(|A|,|B|\) 退化阈值 | 分别使用绝对加相对尺度阈值；尽管本 BDF1 理论应有 `A≥alpha>0`，仍实现文档规定的线性退化路径。 |
| 正根判定 | 要求有限且 `q > positive_root_tol`；不接受零根。多个合格正根时才按 `|q-qⁿ|` 最小选择。 |
| 标量残差 | 使用文档第 21 节归一化残差，分母的 `ε` 取配置值；超过 `scalar_residual_tol` 即失败。 |
| `c₁≤0` 数值检查 | 允许配置的舍入阈值；显著为正视为装配/边界错误并回滚。 |
| 内层收敛 | 同时要求 FE \(L^2\) 相场相对增量小于 `phase_tol` 且 `|q-1|<q_tol`，并允许配置最少/最多松弛步数。 |
| 愈合量 | 在正交点计算 \(\|\min(d^{n+1}-d^n,0)\|_{L^2}\)，不采用节点裁剪近似。 |
| 失败处置 | 一步的候选量只存在局部变量中；根、残差或能量一致性检查失败时不写回 `u,d,q,P`，记录失败诊断并终止返回 `success=false`。 |
| 跨载荷能量 | 每个载荷级记录 reset 前后跳量；不宣称跨载荷级单调。只检查固定载荷内的代理能量律。 |
| 输出符号 | 输出 `q_minus_one=q-1`（带符号）并可由其取绝对值；`phase_increment` 与 `healing` 均为非负 FE \(L^2\) 范数。 |

## 6. 稳定求根与回滚流程

一般二次情形严格使用文档第 14 节：

\[
\widehat q=-\tfrac12\left(B+\operatorname{copysign}(\sqrt D,B)\right),
\qquad q_1=\widehat q/A,
\qquad q_2=C/\widehat q.
\]

`B=0` 时 `copysign` 取正号。若 \(\widehat q=0\)，候选根用不发生消减的等价分支补全，但每个候选仍必须代回原方程检查，不改变原方程。线性/恒等/矛盾退化情形按文档逐项报告。

显著负判别式、没有合格正根、所有根均未通过原式残差、`c₁` 显著为正，或离散代理能量恒等式残差超限，都会触发同一步回滚；不会把判别式、相场或根静默截断到可接受区间。

## 7. 测试计划

1. **Miehe 本构单元测试**：单轴拉伸、单轴压缩、纯剪切、静水状态、零应变、近重根和旋转不变性；核对 \(\epsilon=\epsilon_++\epsilon_-\)、\(\psi_0=\psi_0^++\psi_0^-\) 和面内应力。
2. **退化函数测试**：严格核对 `g(0)=1`、`g(1)=kappa` 及解析导数。
3. **标量求根测试**：唯一正根、双正根、近零判别式、显著负判别式、线性退化、无正根、NaN/Inf 和残差拒绝；包含严重消减误差算例。
4. **矩阵/向量装配测试**：小型 Q1 网格上检查对称性、正定性、`c₁` 负二次型恒等式和四分支边界值。
5. **BDF1 一步测试**：核对四分支重构后的两条离散场方程、标量方程和第 15 节一步代理能量恒等式。
6. **回滚测试**：用故意不合格的根/容差触发失败，逐位比较失败前后的 `u,d,q,P`。
7. **Load–Relax 集成测试**：生成的小网格上执行两个位移载荷级，检查每级 reset、固定载荷内代理能量不增、诊断字段完整且不进行不可逆/区间投影。
8. **L 形烟雾测试**：当根目录 `l_shape.msh` 存在时运行短配置，验证四个命名边界、VTK/CSV 输出和可重复入口；网格缺失时给出明确导入提示。

## 8. 可重复运行

`Project.toml` 精确约束 Julia 1.12.5、Ferrite 1.4.1、FerriteGmsh 1.3.0、
Tensors 1.17.1 与 WriteVTK 1.22.0，其余传递依赖由 `Manifest.toml` 锁定。
把 `l_shape.msh` 放到项目根目录后运行：

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. scripts/run_rlm_bdf1.jl
```

测试入口为：

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

默认结果写到 `data/sims/rlm_bdf1/`：`diagnostics.csv` 保存每个 reset、
接受步和失败回滚行，`load_XXXX.vtu` 保存每个配置输出载荷级的位移与相场。
