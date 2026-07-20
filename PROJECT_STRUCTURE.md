# 项目结构与快速导航

> 这份文件是给首次进入仓库的人和智能体看的“地图”。
> 目标是在不通读全部源码的情况下，先知道项目做什么、从哪里运行、修改某类逻辑应该看哪些文件。

## 1. 项目一句话

这是一个基于 Julia + Ferrite 的二维有限元相场断裂项目，当前包含三类求解路径：

- `Staggered`：位移场和损伤相场交替求解，支持准静态和黏性 BDF1；
- `RLM-BDF1`：Miehe 应变分裂 + RLM 辅助变量的真实时间黏性相场断裂求解器；
- `RLM-CH`：RLM-PC Cahn–Hilliard 方程示例。

主要研究对象是位移 `u`、损伤相场 `d` 以及相关的弹性能、裂纹表面能和黏性耗散。

## 2. 先读什么

推荐按下面顺序建立全局认识：

1. 本文件：了解目录、入口和调用链；
2. `STATE.md`：了解已经实现的功能、当前验证情况和已知限制；
3. `PLAN.md`：了解近期开发目标和设计约束；
4. `src/PffSAV.jl`：确认模块的 include 顺序；
5. 与当前任务对应的求解器文件；
6. `scripts/` 中的运行脚本和 `test/runtests.jl`。

不要一开始从 `Manifest.toml` 或 `data/` 开始阅读：前者是依赖锁定文件，后者主要是网格、结果和图像数据。

## 3. 根目录

```text
PFF/
├── src/                 核心 Julia 模块源码
├── scripts/             可直接运行的算例/绘图脚本
├── test/                自动化测试
├── data/                输入网格、仿真输出、JLD2 结果和图像
├── docs/                主题说明和实现审计文档
├── .agents/             智能体相关的项目约定或工作资料
├── Project.toml         包名、Julia 版本和直接依赖
├── Manifest.toml        完整依赖版本锁定；通常不要手工编辑
├── STATE.md             当前实现状态、验证结果和限制
├── PLAN.md              开发计划与验收标准
├── 结构.md              本文件：项目地图和阅读导航
└── LICENSE              许可证
```

## 4. `src/` 分层

```text
src/
├── PffSAV.jl
├── physics/              材料、本构、能量和物理状态
├── fem/                  网格、自由度、边界约束和有限元组装
├── rlm/                  RLM-BDF1 的配置、组装、Miehe 分裂和标量求解
├── solvers/              面向算例的高层时间推进器
└── utils/                反力、驱动力等通用辅助函数
```

### 4.1 模块入口：`src/PffSAV.jl`

该文件定义 `PffSAV` 模块，并按以下顺序加载源码：

```text
physics → fem → rlm → solvers → utils
```

因此新增类型或函数时，要注意它所依赖的类型必须在更早的文件中定义。当前模块没有单独的 API 层，源码中定义的模块成员会直接成为 `PffSAV` 的可访问符号。

### 4.2 `src/physics/`：物理模型

| 文件 | 作用 | 主要内容 |
|---|---|---|
| `types.jl` | 基础物理参数 | `MaterialParams` 等材料参数 |
| `constitutive.jl` | 本构计算 | 平面应力/应变、三维模式、谱分解、正负弹性能分裂、应力 |
| `energies.jl` | 相场断裂能量 | 弹性能、裂纹表面能及其积分 |
| `ch_types.jl` | Cahn–Hilliard 状态 | `CHParams`、`CHState` |
| `ch_energies.jl` | Cahn–Hilliard 能量 | 双稳势及相关函数 |

修改材料定律、能量公式或积分点计算时，优先从这里开始；不要直接在脚本中复制公式。

### 4.3 `src/fem/`：有限元基础设施

| 文件 | 作用 |
|---|---|
| `setup.jl` | 结构相场算例的网格、位移/损伤自由度和约束；入口是 `setup_tension` |
| `assembly.jl` | 结构相场的位移和损伤方程组装；核心是 `assemble_u!`、`assemble_d!` |
| `ch_setup.jl` | Cahn–Hilliard 网格、自由度和约束；入口是 `setup_ch` |
| `ch_assembly.jl` | Cahn–Hilliard 左端、右端和 RLM 积分组装 |

这一层只负责“如何离散和组装”，高层时间循环位于 `src/solvers/` 或 `src/rlm/solver_bdf1.jl`。

### 4.4 `src/rlm/`：RLM-BDF1 求解器

这是目前最独立、配置最完整的一条求解路径。

| 文件 | 作用 |
|---|---|
| `config.jl` | 配置、状态、诊断和结果类型；包括 `RLMConfig`、`RLMState`、`RLMResult` |
| `miehe2d.jl` | 二维平面应变的 Miehe 谱分解和本构响应 |
| `assembly.jl` | 建立 `RLMProblem`、组装常量矩阵/外力/非线性力和能量 |
| `scalar_solver.jl` | RLM 标量二次方程的稳定求根和残差检查 |
| `solver_bdf1.jl` | 真实时间 BDF1 步进、试算/提交状态、诊断和输出 |

典型调用链：

```text
RLMConfig
  → build_rlm_problem(config)
  → solve_rlm_bdf1(problem)
  → RLMResult + diagnostics + CSV/VTK
```

修改 RLM 时间推进、失败回滚、外载历史或能量记账时，重点阅读 `config.jl`、`assembly.jl` 和 `solver_bdf1.jl`。

### 4.5 `src/solvers/`：高层求解器

| 文件 | 作用 |
|---|---|
| `staggered.jl` | 结构相场的外层时间推进和位移/损伤交替迭代；`eta = 0` 为准静态，`eta > 0` 为黏性 BDF1 |
| `rlm_ch.jl` | Cahn–Hilliard 的 RLM-PC 时间推进 |

`staggered.jl` 的关键状态边界是：上一物理时间步的已接受 `d_old`，不能被时间步内的交替迭代覆盖；只有收敛后才提交新状态。

## 5. 运行入口

从仓库根目录运行：

```powershell
# 全部测试
julia --project=. -e "using Pkg; Pkg.test()"

# 交错相场：修改脚本顶部的 ETA、N_STEPS、DT 等参数
julia --project=. scripts/run_staggered.jl

# RLM-BDF1：修改脚本中的材料、载荷历史和输出配置
julia --project=. scripts/run_rlm_bdf1.jl

# Cahn–Hilliard RLM-PC
julia --project=. scripts/run_ch_rlm.jl
```

脚本通常负责参数、算例选择、绘图和输出；算法实现应放在 `src/`，不要把核心求解逻辑继续堆进脚本。

## 6. 数据目录约定

```text
data/
├── mesh/     Gmsh 输入网格及 Ferrite 序列化缓存（`.jls`）
├── sims/     VTK 等时步场输出
├── jld2/     仿真历史、能量和诊断结果
└── plots/    由脚本生成的图片
```

通常可复用已有网格；网格缓存文件由 `setup.jl` 或 `ch_setup.jl` 根据 `.msh` 文件生成。生成结果时优先写入对应的 `data/` 子目录，不要把输出散落到根目录。

## 7. 按任务定位文件

| 要做的事情 | 首先看 | 然后看 |
|---|---|---|
| 改材料参数或本构公式 | `src/physics/types.jl`、`constitutive.jl` | `src/physics/energies.jl`、相关测试 |
| 改网格、自由度或边界条件 | `src/fem/setup.jl` | `src/fem/assembly.jl`、运行脚本 |
| 改准静态/黏性 Staggered | `src/solvers/staggered.jl` | `src/fem/assembly.jl`、`scripts/run_staggered.jl` |
| 改 RLM-BDF1 时间推进 | `src/rlm/solver_bdf1.jl` | `src/rlm/config.jl`、`src/rlm/assembly.jl` |
| 改 RLM 的 Miehe 分裂 | `src/rlm/miehe2d.jl` | `src/physics/constitutive.jl`、测试 |
| 改 RLM 标量求根 | `src/rlm/scalar_solver.jl` | `test/runtests.jl` |
| 改 Cahn–Hilliard | `src/fem/ch_setup.jl`、`ch_assembly.jl` | `src/solvers/rlm_ch.jl`、`scripts/run_ch_rlm.jl` |
| 改输出、反力或能量诊断 | 对应求解器文件 | `src/utils/utils_fun.jl`、`data/` |
| 判断功能是否已经完成 | `STATE.md` | `PLAN.md`、`docs/` |

## 8. 修改代码时必须保持的约束

- 保持 `src/PffSAV.jl` 的依赖加载顺序，除非同时调整所有依赖关系；
- 时间相关逻辑使用真实物理时间，不要把“时间步编号”误当作时间；
- 已接受状态和当前试算状态必须分开，尤其是 `d_old`、历史场和 RLM 的 `q`；
- 修改 `eta = 0` 的路径时，要确认准静态回归行为没有被改变；
- 新增物理公式、边界条件或状态字段时，同时补充测试或更新 `STATE.md`；
- 不要手工修改 `Manifest.toml`，依赖变更应通过 Julia 的包管理器完成；
- 大型输出文件属于结果数据，不是算法源码；修改算法时先确认是否需要重新生成它们。

## 9. 给智能体的最小工作协议

接手任务时，先回答四个问题再改代码：

1. 任务属于 `physics`、`fem`、`rlm`、`solvers` 还是 `scripts`？
2. 它影响哪条调用链：`Staggered`、`RLM-BDF1` 还是 `RLM-CH`？
3. 哪些状态是“上一时间步已接受”的，哪些是“当前步试算”的？
4. 最小验证是什么：单元测试、模块加载、完整测试，还是一个短算例？

完成修改后，至少运行与改动直接相关的验证，并在交接说明中写明：改了什么、验证了什么、仍有哪些限制。

## 10. 文档维护规则

当目录、入口、核心状态或求解器行为发生变化时，优先更新本文件对应的小节；当“已实现/未实现/验证状态”发生变化时，同时更新 `STATE.md`。本文件描述稳定结构，不记录每次提交的临时细节。
