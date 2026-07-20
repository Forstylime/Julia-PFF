# Miehe 谱分解下的 Mode-I 准静态黏性相场断裂 RLM 方法

## 0. 研究范围

本文研究采用 **Miehe 应变谱分解**和 AT2 裂纹表面能的相场断裂模型，并将松弛拉格朗日乘子方法（Relaxed Lagrange Multiplier, RLM）推广到由以下两部分组成的真实时间耦合系统：

1. 位移场的准静态平衡；

2. 相场变量的黏性 $L^2$ 梯度流。

本文忽略位移惯性，但保留损伤演化的真实时间尺度。外载和 Dirichlet 位移允许随物理时间 $t$ 变化；每个时间步只推进一次耦合状态，不在每个载荷水平把损伤内迭代到稳态。因而本模型描述的是**准静态力学平衡下的率相关黏性损伤/断裂**，而不是含应力波传播的动态断裂。

> [!important] 真实时间与松弛时间
> 时间变量 $t$ 解释为物理时间，$\Delta t$ 由时间分辨率和精度要求确定，黏性参数必须具有可标定的物理量纲。若在每个载荷水平反复推进直到 $d_t\approx0$，算法就退化为 Load--Relax 稳态求解，并失去加载速率信息。

本文暂时忽略裂纹不可逆约束：

$$
d_t\geq 0,
\qquad
d^{n+1}\geq d^n,
$$

以及区间约束

$$
0\leq d\leq1.
$$

因此，当前模型严格来说是一个**可恢复的黏性相场损伤模型**。它可以用于建立 RLM 的能量稳定性、标量可解性和时间收敛理论，但不能直接代表完整的不可逆准静态断裂过程。

> [!warning] 物理范围
> Miehe 谱分解能够抑制纯压缩状态下的不合理损伤，但它并不严格隔离 Mode I。纯剪切会产生正、负主应变，因此正主应变部分仍可能驱动损伤。
>
> 更准确的表述是：
>
> **Miehe-split phase-field fracture under Mode-I loading**，
>
> 而不是“只允许 Mode-I 断裂的本构模型”。

> [!abstract] 理论审校结论
> 在时变体力、Neumann 载荷和 Dirichlet 位移下，连续与离散能量关系应写成“内部能量变化 = 外功输入 - 黏性耗散 - 数值耗散”。只有外功为零时，它才退化为单调耗散律。在线性分支精确求解、载荷时间求积一致且所选标量根可接受的前提下，BDF1/CN 的相应功—能平衡在形式上成立。
>
> 需要特别区分两件事：RLM 的一般梯度流理论是现有结果；把它推广到“准静态位移平衡 + 损伤梯度流”的耦合代数—微分系统，则仍需单独证明。本文后面的能量推导可以作为该证明的框架，但不能直接视为已有 RLM 定理的无条件推论。
>
> $q$ 和 $P$ 是沿真实时间连续传递的数值状态，不应在普通时间步重置。重置只允许作为显式记录能量跳量的数值重启。本文仍未加入裂纹不可逆约束，因此卸载和循环加载只能用于算法验证，不能直接作为不可逆断裂预测。

---

## 1. 基本变量与函数空间

设材料区域为

$$
\Omega\subset\mathbb R^m,
\qquad
m=2\ \text{或}\ 3.
$$

位移场为

$$
\mathbf u:\Omega\rightarrow\mathbb R^m,
$$

相场损伤变量为

$$
d:\Omega\rightarrow\mathbb R,
$$

并采用约定：

$$
d=0
\quad\text{表示完整材料},
$$

$$
d=1
\quad\text{表示完全断裂}.
$$

边界划分为

$$
\partial\Omega
=
\Gamma_D\cup\Gamma_N,
\qquad
\Gamma_D\cap\Gamma_N=\varnothing.
$$

定义位移试探空间和测试空间：

$$
V_{\mathbf u_D(t)}
=
\left\{
\mathbf u\in H^1(\Omega;\mathbb R^m):
\mathbf u=\mathbf u_D(t)
\text{ on }\Gamma_D
\right\},
$$

$$
V_0
=
\left\{
\mathbf v\in H^1(\Omega;\mathbb R^m):
\mathbf v=\mathbf0
\text{ on }\Gamma_D
\right\}.
$$

相场空间取

$$
W=H^1(\Omega).
$$

相场采用自然边界条件：

$$
\nabla d\cdot\mathbf n=0
\qquad
\text{on }\partial\Omega.
$$

---

## 2. Miehe 应变谱分解

### 2.1 小变形应变

小变形应变张量定义为

$$
\boldsymbol\epsilon(\mathbf u)
=
\frac12
\left(
\nabla\mathbf u+\nabla\mathbf u^T
\right).
$$

对称应变张量的谱分解为

$$
\boldsymbol\epsilon
=
\sum_{a=1}^{m}
\epsilon_a
\mathbf n_a\otimes\mathbf n_a,
$$

其中：

- $\epsilon_a$ 为主应变；

- $\mathbf n_a$ 为对应主方向。

定义 Macaulay 括号：

$$
\langle x\rangle_+
=
\frac{x+|x|}{2}
=
\max(x,0),
$$

$$
\langle x\rangle_-
=
\frac{x-|x|}{2}
=
\min(x,0).
$$

正、负应变张量为

$$
\boxed{
\boldsymbol\epsilon_+
=
\sum_{a=1}^{m}
\langle\epsilon_a\rangle_+
\mathbf n_a\otimes\mathbf n_a,
}
$$

$$
\boxed{
\boldsymbol\epsilon_-
=
\sum_{a=1}^{m}
\langle\epsilon_a\rangle_-
\mathbf n_a\otimes\mathbf n_a.
}
$$

因此

$$
\boldsymbol\epsilon
=
\boldsymbol\epsilon_+
+
\boldsymbol\epsilon_-.
$$

---

### 2.2 正负应变能

对于各向同性线弹性材料，未损伤应变能为

$$
\psi_0(\boldsymbol\epsilon)
=
\frac{\lambda}{2}
\left(\operatorname{tr}\boldsymbol\epsilon\right)^2
+
\mu\,
\boldsymbol\epsilon:\boldsymbol\epsilon.
$$

Miehe 谱分解定义

$$
\boxed{
\psi_0^+(\boldsymbol\epsilon)
=
\frac{\lambda}{2}
\left\langle
\operatorname{tr}\boldsymbol\epsilon
\right\rangle_+^2
+
\mu\,
\boldsymbol\epsilon_+:
\boldsymbol\epsilon_+,
}
$$

$$
\boxed{
\psi_0^-(\boldsymbol\epsilon)
=
\frac{\lambda}{2}
\left\langle
\operatorname{tr}\boldsymbol\epsilon
\right\rangle_-^2
+
\mu\,
\boldsymbol\epsilon_-:
\boldsymbol\epsilon_-.
}
$$

二者满足

$$
\psi_0
=
\psi_0^+
+
\psi_0^-.
$$

相应的积极与消极应力为

$$
\boxed{
\boldsymbol\sigma_0^+
=
\lambda
\left\langle
\operatorname{tr}\boldsymbol\epsilon
\right\rangle_+
\mathbf I
+
2\mu\boldsymbol\epsilon_+,
}
$$

$$
\boxed{
\boldsymbol\sigma_0^-
=
\lambda
\left\langle
\operatorname{tr}\boldsymbol\epsilon
\right\rangle_-
\mathbf I
+
2\mu\boldsymbol\epsilon_-.
}
$$

> [!note] 光滑性
> 标量函数
>
> $$
> x\mapsto\langle x\rangle_+^2
> $$
>
> 是 $C^1$ 的，因此 Miehe 能量具有足以支持一阶变分和能量链式法则的光滑性。
>
> 但它一般不是全局 $C^2$ 的；真正的非光滑位置是主应变穿过零的状态。仅仅出现非零重根或主方向基底切换，并不会使各向同性谱能量本身失去一阶可微性，因为谱函数与特征向量的具体选取无关。不过，数值实现二阶谱切线时仍需用重根极限或除差公式，不能直接对单个特征向量求导。
>
> RLM 中显式冻结非线性变分力，因此第一版算法不需要完整的一致谱 Hessian。

---

## 3. 退化函数与 AT2 裂纹能量

建议采用归一化残余刚度退化函数

$$
\boxed{
g_\kappa(d)
=
(1-\kappa)(1-d)^2+\kappa,
}
$$

其中

$$
0<\kappa\ll1.
$$

它满足

$$
g_\kappa(0)=1,
\qquad
g_\kappa(1)=\kappa.
$$

其一阶导数为

$$
\boxed{
g_\kappa'(d)
=
-2(1-\kappa)(1-d).
}
$$

AT2 裂纹表面能为

$$
\boxed{
\mathcal G_c(d)
=
\int_\Omega
\left[
\frac{G_c}{2\ell}d^2
+
\frac{G_c\ell}{2}
|\nabla d|^2
\right]dx.
}
$$

其中：

- $G_c$ 为临界能量释放率；

- $\ell$ 为相场正则化长度。

---

## 4. 时变载荷、内部能量与总势能

设时变体力和 Neumann 载荷对应的线性泛函为

$$
\ell_{\mathrm{ext}}(t;\mathbf u)
=
\int_\Omega
\mathbf b(t)\cdot\mathbf u\,dx
+
\int_{\Gamma_N}
\bar{\mathbf t}(t)\cdot\mathbf u\,ds.
$$

定义内部弹性—断裂能

$$
\boxed{
\begin{aligned}
\mathcal H(\mathbf u,d)
={}&
\int_\Omega
\left[
g_\kappa(d)
\psi_0^+
\bigl(\boldsymbol\epsilon(\mathbf u)\bigr)
+
\psi_0^-
\bigl(\boldsymbol\epsilon(\mathbf u)\bigr)
\right]dx
\\
&+
\int_\Omega
\left[
\frac{G_c}{2\ell}d^2
+
\frac{G_c\ell}{2}
|\nabla d|^2
\right]dx.
\end{aligned}
}
$$

相应的时变总势能为

$$
\boxed{
\mathcal E(t;\mathbf u,d)
=
\mathcal H(\mathbf u,d)
-
\ell_{\mathrm{ext}}(t;\mathbf u).
}
$$

在力控制下，$\mathcal E$ 便于写平衡方程；在真实时间能量验证中，$\mathcal H$ 更直接，因为其变化等于外功减耗散。Dirichlet 边界输入功由反力给出，不能包含在 $\ell_{\mathrm{ext}}$ 中。

有损应力为

$$
\boxed{
\boldsymbol\sigma(\mathbf u,d)
=
g_\kappa(d)\boldsymbol\sigma_0^+(\mathbf u)
+
\boldsymbol\sigma_0^-(\mathbf u).
}
$$

---

## 5. 原始准静态—真实时间黏性系统

### 5.1 位移平衡

对任意
$$
\mathbf v\in V_0,
$$
准静态平衡弱式为
$$
\boxed{
\int_\Omega
\boldsymbol\sigma(\mathbf u,d):
\boldsymbol\epsilon(\mathbf v)\,dx
=
\ell_{\mathrm{ext}}(t;\mathbf v).
}
$$

---

### 5.2 相场化学势

相场变分导数为
$$
\mu_d
=
\frac{\delta\mathcal H}{\delta d}
=
\frac{\delta\mathcal E}{\delta d}.
$$
其弱式为
$$
\boxed{
\begin{aligned}
(\mu_d,\xi)
={}&
\int_\Omega
g_\kappa'(d)
\psi_0^+
\bigl(\boldsymbol\epsilon(\mathbf u)\bigr)
\,\xi\,dx
\\
&+
\int_\Omega
\left[
\frac{G_c}{\ell}d\xi
+
G_c\ell\nabla d\cdot\nabla\xi
\right]dx,
\end{aligned}
}
$$

对任意
$$
\xi\in W.
$$
相场采用黏性 $L^2$ 梯度流：
$$
\boxed{
\eta d_t=-\mu_d,
\qquad
M=\eta^{-1},
}
$$
其中
$$
\eta>0
$$
为损伤黏性参数，$M>0$ 为迁移率。若 $d$ 无量纲且 $\mu_d$ 按单位体积能量计，则 $\eta$ 具有“应力乘时间”的量纲；有限元弱式采用总能量量纲时应按厚度和积分测度保持一致。$\eta$ 或 $M$ 决定真实的损伤时间尺度，不能仅作为任意的数值加速参数。

等价弱式为
$$
\boxed{
\eta(d_t,\xi)
+
(\mu_d,\xi)
=0.
}
$$

---

### 5.3 原始功—能平衡

设 $\mathbf r_D$ 为与当前应力、体力和 Neumann 载荷一致的 Dirichlet 边界反力。外功率定义为

$$
\boxed{
\mathcal P_{\mathrm{ext}}(t)
=
\ell_{\mathrm{ext}}(t;\mathbf u_t)
+
\int_{\Gamma_D}
\mathbf r_D\cdot\dot{\mathbf u}_D\,ds.
}
$$

在位移始终满足准静态平衡且解足够光滑时，内部能量满足

$$
\boxed{
\frac{d}{dt}\mathcal H(\mathbf u,d)
=
\mathcal P_{\mathrm{ext}}(t)
-
\eta\|d_t\|_{L^2}^2
=
\mathcal P_{\mathrm{ext}}(t)
-
M\|\mu_d\|_{L^2}^2.
}
$$

等价地，总势能满足

$$
\boxed{
\frac{d}{dt}\mathcal E(t;\mathbf u,d)
=
\int_{\Gamma_D}\mathbf r_D\cdot\dot{\mathbf u}_D\,ds
-
\partial_t\ell_{\mathrm{ext}}(t;\mathbf u)
-
\eta\|d_t\|_{L^2}^2.
}
$$

因此真实时间加载下不应检查“总能量单调下降”，而应检查累计外功、内部能量和黏性耗散之间的闭合误差。固定载荷和固定 Dirichlet 数据只是 $\mathcal P_{\mathrm{ext}}=0$ 的特殊情形。

---

## 6. RLM 所需的线性—非线性能量分裂

### 6.1 常系数二次部分

定义未损伤弹性双线性形式

$$
\boxed{
a_u(\mathbf u,\mathbf v)
=
\int_\Omega
\boldsymbol\epsilon(\mathbf v):
\mathbb C_0:
\boldsymbol\epsilon(\mathbf u)\,dx.
}
$$

定义相场双线性形式

$$
\boxed{
a_d(d,\xi)
=
\int_\Omega
\left[
\frac{G_c}{\ell}d\xi
+
G_c\ell\nabla d\cdot\nabla\xi
\right]dx.
}
$$

令

$$
\boxed{
\mathcal E_0(t;\mathbf u,d)
=
\frac12a_u(\mathbf u,\mathbf u)
+
\frac12a_d(d,d)
-
\ell_{\mathrm{ext}}(t;\mathbf u).
}
$$

由于

$$
\frac12a_u(\mathbf u,\mathbf u)
=
\int_\Omega
\left[
\psi_0^+
+
\psi_0^-
\right]dx,
$$

原始能量可精确分解为

$$
\boxed{
\mathcal E(t;\mathbf u,d)
=
\mathcal E_0(t;\mathbf u,d)+\mathcal E_1(\mathbf u,d),
}
$$

其中耦合非线性能量为

$$
\boxed{
\mathcal E_1(\mathbf u,d)
=
\int_\Omega
\left[
g_\kappa(d)-1
\right]
\psi_0^+
\bigl(\boldsymbol\epsilon(\mathbf u)\bigr)
\,dx.
}
$$

> [!abstract] 分裂的核心意义
> 该分裂将：
>
> - 未损伤弹性刚度；
>
> - AT2 裂纹表面能；
>
> - 时变外载势能（只改变右端）
>
>
> 全部放入二次/线性部分 $\mathcal E_0$。这里“常系数”指双线性算子不随时间变化；载荷向量允许逐步更新。
>
> Miehe 谱分解和位移—相场耦合全部进入 $\mathcal E_1$。
>
> 因此，全局线性矩阵不依赖于当前损伤变量，也不依赖于谱分解结果。

---

## 7. 非线性变分力

定义

$$
\mathcal N_u(\mathbf u,d;\mathbf v)
=
\delta_{\mathbf u}
\mathcal E_1(\mathbf u,d)[\mathbf v],
$$

$$
\mathcal N_d(\mathbf u,d;\xi)
=
\delta_d
\mathcal E_1(\mathbf u,d)[\xi].
$$

其显式形式为

$$
\boxed{
\mathcal N_u(\mathbf u,d;\mathbf v)
=
\int_\Omega
\left[
g_\kappa(d)-1
\right]
\boldsymbol\sigma_0^+(\mathbf u):
\boldsymbol\epsilon(\mathbf v)\,dx,
}
$$

以及

$$
\boxed{
\mathcal N_d(\mathbf u,d;\xi)
=
\int_\Omega
g_\kappa'(d)
\psi_0^+
\bigl(\boldsymbol\epsilon(\mathbf u)\bigr)
\,\xi\,dx.
}
$$

注意：$\mathcal N_u$ 和 $\mathcal N_d$ 必须来自同一个非线性能量 $\mathcal E_1$

不能分别构造彼此不一致的位移驱动力与相场驱动力，否则 RLM 标量能量方程无法形成严格闭环。

---

## 8. 连续 RLM 重构

引入全局标量乘子

$$
q(t),
\qquad
q(0)=1,
$$

以及松弛参数

$$
\alpha>0.
$$

其中 $q$ 无量纲，而 $\alpha$ 必须具有能量量纲，因为 $\alpha(q^2-1)$ 与 $\mathcal E$ 相加。跨算例比较时宜写成

$$
\alpha=\widehat\alpha E_{\mathrm{ref}},
$$

其中 $\widehat\alpha$ 无量纲，$E_{\mathrm{ref}}$ 是固定参考能量。

### 8.1 位移方程

求

$$
\mathbf u(t)\in V_{\mathbf u_D(t)},
$$

使得

$$
\boxed{
a_u(\mathbf u,\mathbf v)
+
q\,
\mathcal N_u(\mathbf u,d;\mathbf v)
=
\ell_{\mathrm{ext}}(t;\mathbf v),
\qquad
\forall\mathbf v\in V_0.
}
$$

---

### 8.2 相场方程

$$
\boxed{
\eta(d_t,\xi)
+
(\mu,\xi)
=0,
}
$$

$$
\boxed{
(\mu,\xi)
=
a_d(d,\xi)
+
q\,
\mathcal N_d(\mathbf u,d;\xi),
\qquad
\forall\xi\in W.
}
$$

---

### 8.3 RLM 标量关系

$$
\boxed{
\frac{d}{dt}
\mathcal E_1(\mathbf u,d)
+
\alpha
\frac{d}{dt}q^2
=
q
\left[
\mathcal N_u(\mathbf u,d;\mathbf u_t)
+
\mathcal N_d(\mathbf u,d;d_t)
\right].
}
$$

---

### 8.4 连续 RLM 功—能平衡

定义

$$
\boxed{
\widetilde{\mathcal H}(\mathbf u,d,q)
=
\mathcal H(\mathbf u,d)
+
\alpha(q^2-1).
}
$$

时变载荷和时变 Dirichlet 位移下，

$$
\boxed{
\frac{d}{dt}
\widetilde{\mathcal H}(\mathbf u,d,q)
=
\mathcal P_{\mathrm{ext}}(t)
-
\eta\|d_t\|_{L^2}^2
=
\mathcal P_{\mathrm{ext}}(t)
-
M\|\mu\|_{L^2}^2.
}
$$

当 $\mathcal P_{\mathrm{ext}}=0$ 时，才恢复原来的单调耗散结论。真实时间计算应累计 $\mathcal P_{\mathrm{ext}}$，并监测功—能残差。

连续层面，

$$
q(t)\equiv1
$$

是与原始模型一致的解分支。

当

$$
q=1
$$

时，重构系统恢复原始 Miehe 相场模型。

上述推导是形式能量推导。对非齐次时变 Dirichlet 数据，应取提升 $\mathbf z(t)$，写成 $\mathbf u=\mathbf w+\mathbf z(t)$、$\mathbf w_t\in V_0$，并把 $\mathbf z_t$ 对应的项解释为 Dirichlet 反力功。直接把 $\mathbf u_t$ 当作 $V_0$ 测试函数会漏掉边界输入功。

---

## 9. 首选格式：Miehe–RLM-PE(Predicted Energy)–BDF1

BDF1 是当前最适合首先建立完整理论的格式，因为它具有：

- 一步能量律；

- 简单仿射分解；

- 简单二次标量方程；

- 清晰的唯一正根充分条件；

- 不需要两步 $G$-稳定修正能量。

---

### 9.1 外推非线性力

BDF1 中取

$$
\overline U^{n+1}=U^n,
\qquad
U^n=(\mathbf u^n,d^n).
$$

定义冻结的非线性线性泛函

$$
\mathcal N_u^n(\mathbf v)
=
\mathcal N_u(\mathbf u^n,d^n;\mathbf v),
$$

$$
\mathcal N_d^n(\xi)
=
\mathcal N_d(\mathbf u^n,d^n;\xi).
$$

---

### 9.2 BDF1 真实时间离散方程

令 $t_{n+1}=t_n+\Delta t$。$\Delta t$ 是物理时间步，不是固定载荷下的内循环编号。每个时间步只完成一次以下更新。

位移方程为

$$
\boxed{
a_u(\mathbf u^{n+1},\mathbf v)
+
q^{n+1}
\mathcal N_u^n(\mathbf v)
=
\ell_{\mathrm{ext}}(t_{n+1};\mathbf v),
\qquad
\mathbf u^{n+1}=\mathbf u_D(t_{n+1})
\ \text{on }\Gamma_D.
}
$$

相场方程为

$$
\boxed{
\frac{\eta}{\Delta t}
(d^{n+1}-d^n,\xi)
+
a_d(d^{n+1},\xi)
+
q^{n+1}
\mathcal N_d^n(\xi)
=0
}
$$

---

## 10. 仿射分解

令

$$
\boxed{
\mathbf u^{n+1}
=
\mathbf u_a^{n+1}
+
q^{n+1}
\mathbf u_b^{n+1},
}
$$

$$
\boxed{
d^{n+1}
=
d_a^{n+1}
+
q^{n+1}
d_b^{n+1}.
}
$$

### 10.1 位移分支

基准分支满足

$$
\boxed{
a_u(\mathbf u_a^{n+1},\mathbf v)
=
\ell_{\mathrm{ext}}(t_{n+1};\mathbf v).
}
$$

其中

$$
\mathbf u_a^{n+1}
=
\mathbf u_D(t_{n+1})
\qquad\text{on }\Gamma_D.
$$

乘子响应分支满足

$$
\boxed{
a_u(\mathbf u_b^{n+1},\mathbf v)
=
-\mathcal N_u^n(\mathbf v).
}
$$

其中

$$
\mathbf u_b^{n+1}
=
\mathbf0
\qquad\text{on }\Gamma_D.
$$

---

### 10.2 相场分支

基准分支满足

$$
\boxed{
\frac{\eta}{\Delta t}
(d_a^{n+1}-d^n,\xi)
+
a_d(d_a^{n+1},\xi)
=0
}
$$

乘子响应分支满足

$$
\boxed{
\frac{\eta}{\Delta t}
(d_b^{n+1},\xi)
+
a_d(d_b^{n+1},\xi)
=
-\mathcal N_d^n(\xi).
}
$$

---

### 10.3 常系数矩阵

位移矩阵对应

$$
K_u
\leftrightarrow
a_u.
$$

相场矩阵对应

$$
K_d
\leftrightarrow
\frac{\eta}{\Delta t}
(\cdot,\cdot)
+
a_d(\cdot,\cdot).
$$

在固定时间步长下：

- $K_u$ 不随时间改变；

- $K_d$ 不随时间改变；

- 两个矩阵均可在时间循环外预分解；

- 每步只重新组装非线性右端。

时变载荷下，$\mathbf u_a^{n+1}$ 的矩阵不变，但右端和非齐次边界值随 $t_{n+1}$ 更新，因此一般每个时间步都要求解。若 $\Gamma_D$ 的自由度集合固定，Dirichlet 消元后的矩阵仍可复用；若边界分区本身随时间改变，则需重新构造相应矩阵和约束。

Miehe 谱分解只影响 $\mathcal N_u^n,\mathcal N_d^n$，不影响全局矩阵。

---

## 11. RLM-PE 预测非线性能量

定义 $q=1$ 的预测态：

$$
\boxed{
\mathbf u_*^{n+1}
=
\mathbf u_a^{n+1}
+
\mathbf u_b^{n+1},
}
$$

$$
\boxed{
d_*^{n+1}
=
d_a^{n+1}
+
d_b^{n+1}.
}
$$

即

$$
U_*^{n+1}
=
U_a^{n+1}+U_b^{n+1}.
$$

定义预测非线性能量

$$
\boxed{
P^{n+1}
=
\mathcal E_1(U_*^{n+1}).
}
$$

显式写为

$$
\boxed{
P^{n+1}
=
\int_\Omega
\left[
g_\kappa(d_*^{n+1})-1
\right]
\psi_0^+
\bigl(
\boldsymbol\epsilon(\mathbf u_*^{n+1})
\bigr)
\,dx.
}
$$

程序必须保存上一时刻实际使用的预测能量

$$
P^n.
$$

不能在下一步将其任意替换为

$$
\mathcal E_1(U^n),
$$

否则离散 RLM 代理能量无法严格望远镜相消。

初始时取

$$
q^0=1,
$$

$$
P^0=\mathcal E_1(U^0).
$$

---

## 12. BDF1 标量能量方程

离散 RLM 标量关系定义为

$$
\boxed{
\begin{aligned}
&
P^{n+1}-P^n
+
\alpha
\left[
(q^{n+1})^2-(q^n)^2
\right]
\\
&=
q^{n+1}
\left[
\mathcal N_u^n
(\mathbf u^{n+1}-\mathbf u^n)
+
\mathcal N_d^n
(d^{n+1}-d^n)
\right].
\end{aligned}
}
$$

定义

$$
\boxed{
c_0^n
=
\mathcal N_u^n
(\mathbf u_a^{n+1}-\mathbf u^n)
+
\mathcal N_d^n
(d_a^{n+1}-d^n),
}
$$

以及

$$
\boxed{
c_1^n
=
\mathcal N_u^n
(\mathbf u_b^{n+1})
+
\mathcal N_d^n
(d_b^{n+1}).
}
$$

由仿射分解可得

$$
\begin{aligned}
&
\mathcal N_u^n
(\mathbf u^{n+1}-\mathbf u^n)
+
\mathcal N_d^n
(d^{n+1}-d^n)
\\
&=
c_0^n
+
q^{n+1}c_1^n.
\end{aligned}
$$

因此 $q^{n+1}$ 满足二次方程

$$
\boxed{
A_n(q^{n+1})^2
+
B_nq^{n+1}
+
C_n
=
0,
}
$$

其中

$$
\boxed{
A_n
=
\alpha-c_1^n,
}
$$

$$
\boxed{
B_n
=
-c_0^n,
}
$$

$$
\boxed{
C_n
=
P^{n+1}-P^n
-
\alpha(q^n)^2.
}
$$

时变载荷不直接出现在标量关系中，因为外载属于线性部分；它通过当前基准分支 $\mathbf u_a^{n+1}$ 进入 $c_0^n$ 和 $P^{n+1}$。只要载荷是预先给定、与 $q$ 无关的死载荷，二次结构保持不变。若采用随变形方向变化的 follower load 或其他状态相关外力，则必须把相应非线性势和变分力纳入统一分裂并重新推导标量方程。

---

## 13. 标量方程的可解性

判别式为

$$
\boxed{
D_n
=
B_n^2-4A_nC_n.
}
$$

先利用两个响应分支的定义确定 $c_1^n$ 的符号。分别取测试函数 $\mathbf v=\mathbf u_b^{n+1}$ 和 $\xi=d_b^{n+1}$，得到

$$
\mathcal N_u^n(\mathbf u_b^{n+1})
=
-a_u(\mathbf u_b^{n+1},\mathbf u_b^{n+1}),
$$

以及

$$
\mathcal N_d^n(d_b^{n+1})
=
-\frac{\eta}{\Delta t}\|d_b^{n+1}\|_{L^2}^2
-a_d(d_b^{n+1},d_b^{n+1}).
$$

因此

$$
\boxed{
c_1^n
=
-a_u(\mathbf u_b^{n+1},\mathbf u_b^{n+1})
-\frac{\eta}{\Delta t}\|d_b^{n+1}\|_{L^2}^2
-a_d(d_b^{n+1},d_b^{n+1})
\leq0.
}
$$

只要 $\alpha>0$，就自动有

$$
A_n=\alpha-c_1^n\geq\alpha>0.
$$

所以 BDF1 情形不需要额外假设 $\alpha>c_1^n$。在 $q^n\neq0$ 时，唯一正根的一个简单充分条件只剩下 $C_n<0$，即

$$
\boxed{
\alpha
>
\frac{P^{n+1}-P^n}{(q^n)^2}
}
$$

结合 $\alpha>0$，也可写为

$$
\boxed{
\alpha>
\max\left\{
0,\,
\frac{P^{n+1}-P^n}{(q^n)^2}
\right\}.
}
$$

此时 $A_nC_n<0$，故 $D_n>0$，两个实根异号，从而恰有一个正根。注意这只是容易检查的**充分条件**，不是实根存在的必要条件；若它不满足，仍应直接检查判别式和候选根。

真实时间加载速率越高，$P^{n+1}-P^n$ 的单步变化通常越大，因此固定 $\alpha$ 的充分条件可能比固定载荷测试更严格。这是根可解性与时间分辨率共同作用的结果，不能通过每步重置 $q$ 来规避；应优先减小 $\Delta t$，并重新检查 $q$ 一致性和功—能残差。

> [!important] 关于 $\alpha$
> 在一次完整的真实时间计算中，理论分析应使用固定的
>
> $$
> \alpha>0.
> $$
>
> 如果每一步任意改变 $\alpha$，离散 RLM 能量中的 $\alpha[(q^n)^2-1]$ 将不能直接望远镜相消；除非把 $\alpha$ 的变化项显式计入能量平衡，否则原有证明失效。
>
> 实际工作中应按参考能量无量纲化，并通过先验估计或参数扫描寻找一个同时满足：
>
> - 二次方程可解；
>
> - $q$ 接近 1；
>
> - 浮点条件良好
>
>
> 的固定参数。增大 $\alpha$ 通常会把根推向 $\pm q^n$，但“越大越好”并不自动成立，仍需检查尺度和舍入误差。

---

## 14. 稳定二次求根

不建议直接同时使用标准公式

$$
q_\pm
=
\frac{-B\pm\sqrt{D}}{2A},
$$

因为当

$$
B^2\gg4AC
$$

时可能发生严重消减误差。

一般情形建议先计算

$$
\boxed{
\widehat q
=
-\frac12
\left(
B+\operatorname{copysign}(\sqrt D,B)
\right),
}
$$

其中约定 $\operatorname{copysign}(\sqrt D,0)=+\sqrt D$。随后令

$$
\boxed{
q_1=\frac{\widehat q}{A},
\qquad
q_2=\frac{C}{\widehat q}.
}
$$

这避免了原公式在 $B=0$ 时因 $\operatorname{sign}(0)=0$ 而产生 $\widehat q=0$ 和除零错误。

在第 13 节的 $A>0,\ C<0$ 条件下，正根还可直接稳定地计算为

$$
\boxed{
q_+
=
\begin{cases}
\dfrac{-B+\sqrt D}{2A}, & B<0,\\[6pt]
\dfrac{2C}{-B-\sqrt D}, & B\geq0.
\end{cases}
}
$$

求根流程：

1. 若
    $$
    D\in[-\varepsilon_D,0),
    $$
    则视为浮点误差并令
    $$
    D=0;
    $$

2. 若
    $$
    |A|<\varepsilon_A,
    $$
    则退化为线性方程 $Bq+C=0$；若同时 $|B|<\varepsilon_B$，还必须单独判断 $C$ 是否也在容差内；

3. 若 $D<-\varepsilon_D$，报告“无实根”，不能把它静默截断为零；

4. 删除 NaN、无穷根以及不满足原二次方程残差容差的候选根；

5. 若存在唯一正根，则选取该根；

6. 只有在理论上允许多个候选实根时，才选择
    $$
    |q^{n+1}-q^n|
    $$
    最小者；不要用“最接近”规则覆盖第 13 节已经保证的唯一正根。

---

## 15. BDF1 离散功—能平衡

真实时间计算优先使用不含载荷势的离散内部代理能量

$$
\boxed{
\widetilde{\mathcal H}_{\mathrm{PE}}^n
=
\frac12a_u(\mathbf u^n,\mathbf u^n)
+
\frac12a_d(d^n,d^n)
+P^n
+\alpha\left[(q^n)^2-1\right].
}
$$

记

$$
\Delta\mathbf u^n=\mathbf u^{n+1}-\mathbf u^n,
\qquad
\Delta d^n=d^{n+1}-d^n.
$$

在离散层面定义 $t_{n+1}$ 时刻的 Dirichlet 反力泛函

$$
\boxed{
\mathcal R_D^{n+1}(\mathbf z)
=
a_u(\mathbf u^{n+1},\mathbf z)
+q^{n+1}\mathcal N_u^n(\mathbf z)
-\ell_{\mathrm{ext}}(t_{n+1};\mathbf z),
}
$$

其中 $\mathbf z$ 是边界增量 $\mathbf u_D(t_{n+1})-\mathbf u_D(t_n)$ 的任意有限元提升。由于平衡残差在 $V_0$ 上为零，$\mathcal R_D^{n+1}(\mathbf z)$ 与提升的内部延拓无关。BDF1 端点外功定义为

$$
\boxed{
W_{\mathrm{ext}}^{n+1}
=
\ell_{\mathrm{ext}}(t_{n+1};\Delta\mathbf u^n)
+
\mathcal R_D^{n+1}(\mathbf z).
}
$$

若线性分支可解、标量二次方程选取了可接受实根且外功采用上述一致离散，则

$$
\boxed{
\begin{aligned}
&\widetilde{\mathcal H}_{\mathrm{PE}}^{n+1}
-\widetilde{\mathcal H}_{\mathrm{PE}}^n
\\
={}&W_{\mathrm{ext}}^{n+1}
-\frac12a_u(\Delta\mathbf u^n,\Delta\mathbf u^n)
-\frac12a_d(\Delta d^n,\Delta d^n)
-\frac{\eta}{\Delta t}\|\Delta d^n\|_{L^2}^2.
\end{aligned}
}
$$

最后一项逼近真实黏性耗散

$$
\int_{t_n}^{t_{n+1}}\eta\|d_t\|_{L^2}^2\,dt,
$$

而两个二次型是 BDF1 的数值耗散。只有 $W_{\mathrm{ext}}^{n+1}=0$ 时，才有 $\widetilde{\mathcal H}_{\mathrm{PE}}^{n+1}\leq\widetilde{\mathcal H}_{\mathrm{PE}}^n$。这里的“无条件稳定”是指功—能估计本身不要求时间步长限制，但不表示任意大 $\Delta t$ 都能准确解析加载历程、起裂时刻或损伤时间尺度；同时仍要求标量方程存在可接受实根。

---

## 16. RLM 代理能量与物理内部能量的一致性

原始物理内部能量为

$$
\mathcal H(U^n)
=
\frac12a_u(\mathbf u^n,\mathbf u^n)
+
\frac12a_d(d^n,d^n)
+
\mathcal E_1(U^n).
$$

二者差值为

$$
\boxed{
\widetilde{\mathcal H}_{\mathrm{PE}}^n
-
\mathcal H(U^n)
=
P^n-\mathcal E_1(U^n)
+
\alpha
\left[
(q^n)^2-1
\right].
}
$$

由于

$$
U^n
=
U_a^n+q^nU_b^n,
$$

而

$$
U_*^n
=
U_a^n+U_b^n,
$$

有

$$
\boxed{
U_*^n-U^n
=
(1-q^n)U_b^n.
}
$$

若 $\mathcal E_1$ 在数值解邻域局部 Lipschitz，则

$$
\left|
P^n-\mathcal E_1(U^n)
\right|
\leq
L_n
|1-q^n|
|U_b^n|.
$$

因此

$$
\boxed{
\begin{aligned}
\left|
\widetilde{\mathcal H}_{\mathrm{PE}}^n
-
\mathcal H(U^n)
\right|
\leq{}&
L_n
|1-q^n|
|U_b^n|
\\
&+
\alpha
\left|
(q^n)^2-1
\right|.
\end{aligned}
}
$$

只要

$$
q^n\rightarrow1,
$$

且

$$
U_b^n
$$

保持有界，RLM 内部代理能量便收敛到原始物理内部能量。总势能版本只需在二者中同时减去 $\ell_{\mathrm{ext}}(t_n;\mathbf u^n)$，但真实时间功—能验证仍应优先使用 $\mathcal H$。

> [!tip] 数值验证重点
> 数值实验不应只画两条能量曲线，而应同时输出：
>
> $$
> |q^n-1|,
> $$
>
> $$
> |P^n-\mathcal E_1(U^n)|,
> $$
>
> $$
> \left|
> \widetilde{\mathcal H}_{\mathrm{PE}}^n
> -
> \mathcal H(U^n)
> \right|.
> $$
>
> 这样才能验证 RLM 代理能量趋近物理内部能量的实际机制。

---

## 17. 二阶扩展：Miehe–RLM-PE–CN

CN 比 BDF2 更适合作为第一个二阶格式，因为 CN 可以保持一步能量律，而 BDF2 通常只能得到包含两步历史的 $G$-修正能量。

这里的“二阶”首先是形式精度结论。严格二阶收敛还要求解具有足够高的时间正则性，并要求预测能量和所选 $q$ 根分支保持二阶一致；若主应变在时间步内穿过零，还需额外分析 Miehe 谱力的有限光滑性。

### 17.1 中点和外推

定义

$$
U^{n+\frac12}
=
\frac12
\left(
U^{n+1}+U^n
\right),
$$

$$
q^{n+\frac12}
=
\frac12
\left(
q^{n+1}+q^n
\right).
$$

二阶显式外推为

$$
\boxed{
\overline U^{n+\frac12}
=
\frac32U^n
-
\frac12U^{n-1}.
}
$$

在该外推状态上冻结

$$
\mathcal N_u^{n+\frac12},
\qquad
\mathcal N_d^{n+\frac12}.
$$

对时变载荷定义与 CN 一致的中点数据

$$
\ell_{\mathrm{ext}}^{n+\frac12}(\mathbf v)
\approx
\ell_{\mathrm{ext}}(t_{n+\frac12};\mathbf v),
\qquad
\mathbf u_D^{n+\frac12}
=
\frac12\left[\mathbf u_D(t_{n+1})+\mathbf u_D(t_n)\right].
$$

---

### 17.2 CN 离散方程

位移方程：

$$
\boxed{
a_u(
\mathbf u^{n+\frac12},
\mathbf v
)
+
q^{n+\frac12}
\mathcal N_u^{n+\frac12}(\mathbf v)
=
\ell_{\mathrm{ext}}^{n+\frac12}(\mathbf v).
}
$$

相场方程：

$$
\boxed{
\eta\left(
\frac{d^{n+1}-d^n}{\Delta t},
\xi
\right)
+
a_d(
d^{n+\frac12},
\xi
)
+
q^{n+\frac12}
\mathcal N_d^{n+\frac12}(\xi)
=0
}
$$

标量关系：

$$
\boxed{
\begin{aligned}
&
P^{n+1}-P^n
+
\alpha
\left[
(q^{n+1})^2-(q^n)^2
\right]
\\
&=
q^{n+\frac12}
\left[
\mathcal N_u^{n+\frac12}
(\mathbf u^{n+1}-\mathbf u^n)
+
\mathcal N_d^{n+\frac12}
(d^{n+1}-d^n)
\right].
\end{aligned}
}
$$

---

### 17.3 CN 仿射分支

仍取

$$
\mathbf u^{n+1}
=
\mathbf u_a^{n+1}
+
q^{n+1}\mathbf u_b^{n+1},
$$

$$
d^{n+1}
=
d_a^{n+1}
+
q^{n+1}d_b^{n+1}.
$$

位移基准分支：

$$
\boxed{
\begin{aligned}
a_u(\mathbf u_a^{n+1},\mathbf v)
={}&
2\ell_{\mathrm{ext}}^{n+\frac12}(\mathbf v)
-
a_u(\mathbf u^n,\mathbf v)
\\
&-
q^n
\mathcal N_u^{n+\frac12}(\mathbf v).
\end{aligned}
}
$$

并施加

$$
\mathbf u_a^{n+1}=\mathbf u_D(t_{n+1})
\qquad\text{on }\Gamma_D.
$$

位移乘子分支：

$$
\boxed{
a_u(\mathbf u_b^{n+1},\mathbf v)
=
-\mathcal N_u^{n+\frac12}(\mathbf v).
}
$$

其中 $\mathbf u_b^{n+1}=\mathbf0$ on $\Gamma_D$。

相场基准分支：

$$
\boxed{
\begin{aligned}
&
\frac{\eta}{\Delta t}
(d_a^{n+1},\xi)
+
\frac12
a_d(d_a^{n+1},\xi)
\\
={}&
\frac{\eta}{\Delta t}
(d^n,\xi)
-
\frac12
a_d(d^n,\xi)
-
\frac12q^n
\mathcal N_d^{n+\frac12}(\xi).
\end{aligned}
}
$$

相场乘子分支：

$$
\boxed{
\frac{\eta}{\Delta t}
(d_b^{n+1},\xi)
+
\frac12
a_d(d_b^{n+1},\xi)
=
-\frac12
\mathcal N_d^{n+\frac12}(\xi).
}
$$

---

### 17.4 CN 二次方程

定义

$$
c_0^n
=
\mathcal N_u^{n+\frac12}
(\mathbf u_a^{n+1}-\mathbf u^n)
+
\mathcal N_d^{n+\frac12}
(d_a^{n+1}-d^n),
$$

$$
c_1^n
=
\mathcal N_u^{n+\frac12}
(\mathbf u_b^{n+1})
+
\mathcal N_d^{n+\frac12}
(d_b^{n+1}).
$$

则 $q^{n+1}$ 满足

$$
A_n(q^{n+1})^2
+
B_nq^{n+1}
+
C_n
=
0,
$$

其中

$$
\boxed{
A_n
=
\alpha-\frac12c_1^n,
}
$$

$$
\boxed{
B_n
=
-\frac12
\left(
c_0^n+q^nc_1^n
\right),
}
$$

$$
\boxed{
C_n
=
P^{n+1}
-
P^n
-
\alpha(q^n)^2
-
\frac12q^nc_0^n.
}
$$

CN 响应分支同样给出

$$
\mathcal N_u^{n+\frac12}(\mathbf u_b^{n+1})
=
-a_u(\mathbf u_b^{n+1},\mathbf u_b^{n+1}),
$$

$$
\mathcal N_d^{n+\frac12}(d_b^{n+1})
=
-\frac{2\eta}{\Delta t}\|d_b^{n+1}\|_{L^2}^2
-a_d(d_b^{n+1},d_b^{n+1}),
$$

所以 $c_1^n\leq0$，从而

$$
A_n=\alpha-\frac12c_1^n\geq\alpha>0.
$$

若 $q^n\neq0$，保证 CN 二次方程存在唯一正根的一个简单充分条件为

$$
\boxed{
\alpha>
\max\left\{
0,\,
\frac{
P^{n+1}-P^n-\frac12q^nc_0^n
}{
(q^n)^2
}
\right\}.
}
$$

这只是充分条件；条件不满足时仍可通过 $D_n\geq0$ 和候选根筛选判断可解性。

---

### 17.5 CN 离散功—能平衡

定义中点反力泛函

$$
\mathcal R_D^{n+\frac12}(\mathbf z)
=
a_u(\mathbf u^{n+\frac12},\mathbf z)
+q^{n+\frac12}\mathcal N_u^{n+\frac12}(\mathbf z)
-\ell_{\mathrm{ext}}^{n+\frac12}(\mathbf z),
$$

以及中点外功

$$
W_{\mathrm{ext,CN}}^{n+1}
=
\ell_{\mathrm{ext}}^{n+\frac12}(\Delta\mathbf u^n)
+\mathcal R_D^{n+\frac12}(\mathbf z).
$$

则 CN 满足

$$
\boxed{
\widetilde{\mathcal H}_{\mathrm{PE}}^{n+1}
-
\widetilde{\mathcal H}_{\mathrm{PE}}^n
=
W_{\mathrm{ext,CN}}^{n+1}
-\frac{\eta}{\Delta t}
\|d^{n+1}-d^n\|_{L^2}^2.
}
$$

与 BDF1 相比，CN 不包含额外的位移和相场二次数值耗散项；在外功为零时才退化为单调耗散律。

首步可使用 BDF1 启动。

---

## 18. 可选扩展：方向二次 RLM-Q

RLM-PE 使用

$$
P^{n+1}
=
\mathcal E_1(U_*^{n+1})
$$

作为预测能量。

为了提高代理能量与真实非线性能量的一致性，可以沿仿射方向构造二次代理。

定义

$$
U^{n+1}(q)
=
U_a^{n+1}+qU_b^{n+1},
$$

$$
U_*^{n+1}
=
U^{n+1}(1).
$$

定义方向一阶导数

$$
\boxed{
G_n
=
D\mathcal E_1(U_*^{n+1})
[U_b^{n+1}].
}
$$

显式形式为

$$
\boxed{
\begin{aligned}
G_n
={}&
\int_\Omega
\left[
g_\kappa(d_*)-1
\right]
\boldsymbol\sigma_0^+(\mathbf u_*):
\boldsymbol\epsilon(\mathbf u_b)\,dx
\\
&+
\int_\Omega
g_\kappa'(d_*)
d_b
\psi_0^+(\mathbf u_*)\,dx.
\end{aligned}
}
$$

定义方向二次代理

$$
\boxed{
\begin{aligned}
P_Q^{n+1}(q)
={}&
\mathcal E_1(U_*^{n+1})
+
(q-1)G_n
\\
&+
\frac{S_n}{2}(q-1)^2.
\end{aligned}
}
$$

其中

$$
S_n\geq0
$$

为方向稳定化参数。

由于 Miehe 谱能量的精确二阶导数涉及谱投影导数和重根处理，第一版 RLM-Q 不建议使用精确谱 Hessian。可以采用：

- 固定标量稳定化；

- 数值方向差分；

- 方向曲率上界；

- 自适应非负 $S_n$。

RLM-Q 仍产生关于 $q^{n+1}$ 的二次方程。

> [!warning] RLM-Q 尚未闭合
> 仅规定 $S_n\geq0$ 并不足以自动得到代理能量一致性、二次方程可解性或唯一正根。要形成完整格式，还必须明确上一时刻保存的是哪一个 $P_Q^n$，并在同一个标量关系中使用它；随后重新推导二次系数和功—能律。因此，本节目前只能视为研究设想，不能与前面的 RLM-PE 结论等价使用。

---

## 19. 真实时间加载：Time-Marching Workflow

给定物理时间区间和时间网格

$$
0=t_0<t_1<\cdots<t_N=T,
$$

以及连续或分段光滑的加载历史

$$
\mathbf b(t),
\qquad
\bar{\mathbf t}(t),
\qquad
\mathbf u_D(t).
$$

实际 Mode-I 位移控制可写成

$$
\mathbf u_D(t)=\lambda(t)\widehat{\mathbf u}_D,
$$

其中 $\lambda(t)$ 是具有真实时间单位的加载程序，例如恒定位移速率、保持、卸载或循环段。

### 19.1 初值与辅助变量

构造满足 $\mathbf u_D(0)$ 的初始平衡位移 $\mathbf u^0$，给定 $d^0$，并仅在计算初始时刻设置

$$
\boxed{
q^0=1,
\qquad
P^0=\mathcal E_1(\mathbf u^0,d^0).
}
$$

随后 $\mathbf u^n,d^n,q^n,P^n$ 都作为真实时间状态连续传递。普通时间步禁止执行

$$
q^n\leftarrow1,
\qquad
P^n\leftarrow\mathcal E_1(U^n),
$$

否则会人为改变 RLM 动力学和功—能记账。若因重启、容错或重新网格化必须重新初始化，应单独记录代理能量跳量。

### 19.2 单步推进，而非固定载荷内循环

从 $t_n$ 到 $t_{n+1}$：

1. 计算 $t_{n+1}$ 的 BDF1 载荷，或 $t_{n+\frac12}$ 的 CN 载荷；

2. 由 $U^n$（CN 还使用 $U^{n-1}$）冻结非线性力；

3. 求解四个仿射分支和一个标量二次方程；

4. 恢复 $U^{n+1}$，累计该步外功与黏性耗散；

5. 直接进入 $t_{n+2}$。

每个物理时间步只执行一次 RLM 时间更新。线性求解器为了达到代数残差容差所做的迭代不属于“损伤内循环”；但不能在保持 $t_{n+1}$ 和载荷不变的情况下反复执行多个 RLM 时间步直到 $d_t\approx0$，否则会把率相关模型重新变成稳态优化器。

### 19.3 时间步、黏性与加载速率

$\eta$、$\Delta t$ 和加载时间尺度必须分开处理：

- $\eta$ 是材料/正则化模型参数，应通过松弛试验、加载速率试验或目标裂纹速度标定；

- $\Delta t$ 是数值分辨率参数，应通过时间步收敛确定；

- $\lambda(t)$ 决定外部加载速率，是算例输入而不是迭代计数。

即使能量格式对 $\Delta t$ 无条件稳定，过大的时间步仍会错过起裂时刻、峰值反力和快速损伤阶段。宜限制单步位移增量和单步损伤增量，并在

$$
\|d^{n+1}-d^n\|_{L^2}
$$

或功—能残差突增时减小 $\Delta t$。改变 $\Delta t$ 会改变相场有效矩阵 $K_d$，通常需要重新分解；$\alpha$ 原则上保持固定。

### 19.4 准静态假设的适用范围

本路线没有位移动能项，要求加载和损伤的主要时间尺度显著慢于结构的弹性波传播/振动时间尺度。若关注冲击、应力波或高速裂纹传播，应改用包含

$$
\rho\mathbf u_{tt}
$$

和动能的动态模型，并重新推导 RLM 格式；不能仅把本节的 $t$ 改名为真实时间来替代惯性。

---

## 20. FEM 实现 Workflow

### Step 0：预处理

组装并分解：

$$
K_u,
\qquad
K_d.
$$

存储：

- $\mathbf u^n$；$d^n$；$q^n$；$P^n$；
- CN 所需的 $U^{n-1}$。
- 物理参数 $\eta$、加载历史和当前物理时间；
- 累计外功、累计黏性耗散和累计数值耗散。

只在 $t_0$ 初始化 $q^0=1$ 和 $P^0=\mathcal E_1(U^0)$。若采用固定 $\Delta t$，$K_u,K_d$ 可在时间循环外预分解；自适应改变 $\Delta t$ 时需更新并重新分解 $K_d$。

---

### Step 1：更新时间数据并构造外推状态

令 $t_{n+1}=t_n+\Delta t$，组装 BDF1 所需的

$$
\ell_{\mathrm{ext}}(t_{n+1};\cdot),
\qquad
\mathbf u_D(t_{n+1}),
$$

或 CN 所需的中点数据。随后构造外推状态：

BDF1：

$$
\overline U^{n+1}=U^n.
$$

CN：

$$
\overline U^{n+\frac12}
=
\frac32U^n-\frac12U^{n-1}.
$$

---

### Step 2：正交点谱分解

在每个正交点计算

$$
\boldsymbol\epsilon(\overline{\mathbf u}),
$$

并求解对称特征值问题

$$
\boldsymbol\epsilon
\mathbf n_a
=
\epsilon_a\mathbf n_a.
$$

构造

$$
\boldsymbol\epsilon_+,
\qquad
\boldsymbol\epsilon_-,
$$

$$
\psi_0^+,
\qquad
\psi_0^-,
$$

$$
\boldsymbol\sigma_0^+,
\qquad
\boldsymbol\sigma_0^-.
$$

> [!tip] 谱分解实现
> 不建议存储单个特征向量作为历史变量。
>
> 特征向量存在：
>
> - 任意符号；
>
> - 排序变化；
>
> - 重根不唯一性。
>
>
> 每次只需由当前应变构造谱张量
>
> $$
> \boldsymbol\epsilon_+,
> \quad
> \boldsymbol\epsilon_-.
> $$

---

### Step 3：组装非线性右端

组装

$$
\mathcal N_u^n,
$$

$$
\mathcal N_d^n.
$$

位移右端的单元贡献为

$$
-\int_{\Omega_e}
\left[
g_\kappa(\bar d)-1
\right]
\boldsymbol\sigma_0^+(\bar{\mathbf u})
:
\boldsymbol\epsilon(\mathbf v)\,dx.
$$

相场右端的单元贡献为

$$
-\int_{\Omega_e}
g_\kappa'(\bar d)
\psi_0^+(\bar{\mathbf u})
\,\xi\,dx.
$$

具体正负号根据仿射分支方程的右端定义保持一致。

---

### Step 4：求解四个仿射分支

求解

$$
\mathbf u_a,
\qquad
\mathbf u_b,
\qquad
d_a,
\qquad
d_b.
$$

由于矩阵已经分解，每一步主要成本为四次前代/回代。$\mathbf u_a^{n+1}$ 使用当前时刻载荷和 Dirichlet 数据，$\mathbf u_b^{n+1}$ 使用齐次 Dirichlet 边界。

---

### Step 5：构造预测态

$$
\mathbf u_*
=
\mathbf u_a+\mathbf u_b,
$$

$$
d_*
=
d_a+d_b.
$$

---

### Step 6：计算预测非线性能量

在预测态正交点重新计算

$$
\psi_0^+(\mathbf u_*),
$$

并积分

$$
P^{n+1}
=
\int_\Omega
[g_\kappa(d_*)-1]
\psi_0^+(\mathbf u_*)\,dx.
$$

注意：

$$
P^{n+1}
$$

不能通过自由度向量直接点乘得到，必须在正交点进行空间积分。

---

### Step 7：计算标量系数

计算

$$
c_0^n,
\qquad
c_1^n,
$$

随后组装

$$
A_n,
\qquad
B_n,
\qquad
C_n.
$$

同时检查理论符号

$$
c_1^n\leq0,
\qquad
A_n\geq\alpha>0.
$$

若该检查明显失败，通常意味着响应分支右端符号、边界条件或内积装配不一致。

---

### Step 8：解二次方程

计算判别式，进行浮点容差处理，并选择物理根：

$$
q^{n+1}.
$$

---

### Step 9：恢复物理场

$$
\boxed{
\mathbf u^{n+1}
=
\mathbf u_a
+
q^{n+1}\mathbf u_b,
}
$$

$$
\boxed{
d^{n+1}
=
d_a
+
q^{n+1}d_b.
}
$$

---

### Step 10：保存状态

保存

$$
q^{n+1},
\qquad
P^{n+1},
\qquad
U^{n+1}.
$$

更新历史变量进入下一步。

同一步还必须由当前反力和位移增量计算 $W_{\mathrm{ext}}^{n+1}$，并累计

$$
W_{\mathrm{ext}}^{0\to n+1},
\qquad
\mathcal D_{\mathrm{vis}}^{0\to n+1},
\qquad
\mathcal D_{\mathrm{num}}^{0\to n+1}.
$$

随后直接推进到下一个物理时间步；不检查 $d$ 是否已经在当前载荷下松弛到稳态，也不重置 $q^{n+1},P^{n+1}$。

---

## 21. 必须记录的诊断量

每一步至少输出：

$$
t_n,
\qquad
\Delta t_n,
\qquad
\lambda(t_n),
\qquad
\dot\lambda(t_n),
$$

$$
q^n,
$$

$$
|q^n-1|,
$$

$$
D_n=B_n^2-4A_nC_n,
$$

$$
c_1^n,
\qquad
A_n-\alpha,
$$

$$
r_q^n
=
\frac{
\left|A_n(q^{n+1})^2+B_nq^{n+1}+C_n\right|
}{
|A_n|(q^{n+1})^2+|B_n||q^{n+1}|+|C_n|+\varepsilon
},
$$

$$
\mathcal H(U^n),
\qquad
\mathcal E(t_n;U^n),
$$

$$
\widetilde{\mathcal H}_{\mathrm{PE}}^n,
$$

$$
P^n-\mathcal E_1(U^n),
$$

$$
\left|
\widetilde{\mathcal H}_{\mathrm{PE}}^n
-
\mathcal H(U^n)
\right|,
$$

$$
|d^{n+1}-d^n|_{L^2},
$$

$$
\min_\Omega d^n,
\qquad
\max_\Omega d^n.
$$

真实时间计算还必须输出当前反力、单步外功和耗散：

$$
\mathbf R_D^n,
\qquad
W_{\mathrm{ext}}^n,
\qquad
\mathcal D_{\mathrm{vis}}^n
=
\frac{\eta}{\Delta t_n}\|d^n-d^{n-1}\|_{L^2}^2,
$$

以及累计功—能残差

$$
\boxed{
r_{\mathrm{energy}}^n
=
\widetilde{\mathcal H}_{\mathrm{PE}}^n
-\widetilde{\mathcal H}_{\mathrm{PE}}^0
-\sum_{j=1}^nW_{\mathrm{ext}}^j
+\sum_{j=1}^n
\left(
\mathcal D_{\mathrm{vis}}^j
+\mathcal D_{\mathrm{num}}^j
\right).
}
$$

对 CN 取 $\mathcal D_{\mathrm{num}}^j=0$；对 BDF1 使用第 15 节的两个二次型。应同时报告归一化残差，避免总能量尺度过大掩盖装配错误。

由于当前没有不可逆约束，还应记录愈合量：

$$
\boxed{
H_{\mathrm{heal}}^{n+1}
=
\left|
\min(
d^{n+1}-d^n,
0
)
\right|_{L^2}.
}
$$

不要简单执行

$$
d^{n+1}
\leftarrow
\max(d^{n+1},d^n),
$$

因为这种事后投影通常会破坏已建立的 RLM 离散能量恒等式。

---

## 22. 建议建立的数学结果（尚待严格证明）

本节列出的是研究目标，而不是已经由前述代数推导自动得到的定理。特别是，现有一般 RLM 结果不能替代对准静态位移约束、非齐次边界提升和 Miehe 谱非线性的专门分析。

### 引理 1：线性分支唯一可解

在以下条件下：

- $\Gamma_D$ 具有正测度；

- $\mathbb C_0$ 正定；

- $G_c>0$；

- $\ell>0$；

- $\eta>0$；

- $\Delta t>0$；

证明四个仿射线性分支唯一可解。

---

### 引理 2：仿射降维

证明

$$
U^{n+1}
=
U_a^{n+1}
+
q^{n+1}U_b^{n+1}.
$$

---

### 引理 3：标量二次结构

证明 RLM-PE 标量方程关于 $q^{n+1}$ 至多为二次方程。

---

### 定理 1：BDF1 离散功—能稳定性

证明

$$
\widetilde{\mathcal H}_{\mathrm{PE}}^{n+1}
-
\widetilde{\mathcal H}_{\mathrm{PE}}^n
=
W_{\mathrm{ext}}^{n+1}
-
\mathcal D_{\mathrm{vis}}^{n+1}
-
\mathcal D_{\mathrm{num}}^{n+1}.
$$

并证明离散反力功与非齐次 Dirichlet 提升选取无关。零外功时的单调耗散只是该定理的推论。

---

### 定理 2：唯一正根充分条件

在

$$
A_n>0,
\qquad
C_n<0
$$

条件下，证明标量方程存在唯一正根。

对本文 BDF1/CN 响应分支，已经可由分支方程直接证明 $c_1^n\leq0$，因此 $A_n>0$ 不需要再假设 $c_1^n$ 的统一上界。真正需要研究的是：$q^n$ 是否一致远离零，以及下列组合能否由固定 $\alpha$ 控制：

$$
P^{n+1}-P^n
$$

和 CN 中的

$$
P^{n+1}-P^n-\frac12q^nc_0^n.
$$

从而给出固定 $\alpha$ 的先验充分条件。

---

### 定理 3：RLM 代理能量一致性

在固定 $\alpha$、选取与 $q=1$ 连续相连的根分支、$q^n$ 一致远离零以及预测能量一致的条件下，研究是否可证明

$$
\max_{0\leq n\leq N}
|q^n-1|
\rightarrow0
\qquad
\text{as }\Delta t\rightarrow0,
$$

以及

$$
\max_n
\left|
\widetilde{\mathcal H}_{\mathrm{PE}}^n
-
\mathcal H(U^n)
\right|
\rightarrow0.
$$

不能仅由离散功—能稳定性推出 $q^n\to1$；还需要一致性估计和根分支控制。

---

### 定理 4：BDF1 时间一阶收敛

在解与时变载荷足够光滑、准静态平衡解关于 $d$ 稳定、非齐次边界提升时间正则、离散初值一致且非线性谱力满足所需 Lipschitz 估计时，证明例如

$$
\|\mathbf u(t_n)-\mathbf u^n\|_{H^1}
+
\|d(t_n)-d^n\|_{L^2}
+
|q(t_n)-q^n|
\leq
C\Delta t.
$$

---

### 定理 5：CN 二阶收敛

CN 的二阶结论还要求载荷和 Dirichlet 数据采用二阶一致的中点求积，并需要更强的时间正则性。由于 Miehe 谱能量一般只有 $C^1$，若解轨道发生主应变过零，经典二阶局部截断误差证明可能失效。可在“活跃谱分支在时间步内不切换”或采用光滑化谱分裂的附加条件下证明

$$
\|\mathbf u(t_n)-\mathbf u^n\|_{H^1}
+
\|d(t_n)-d^n\|_{L^2}
+
|q(t_n)-q^n|
\leq
C\Delta t^2.
$$

---

## 23. 数值验证路线

### 23.1 制造解测试

使用平滑外载和不会产生尖锐局部化的参数，验证：

- BDF1 一阶；

- CN 二阶；

- $q-1$ 的时间收敛；

- RLM 代理能量误差的时间收敛。

- 时变 Dirichlet/Neumann 数据下累计功—能残差的收敛。

---

### 23.2 真实时间斜坡—保持加载

施加“线性位移斜坡 + 位移保持”的开口历史。在斜坡段验证外功输入，在保持段验证损伤黏性松弛，并检查

$$
\Delta\widetilde{\mathcal H}_{\mathrm{PE}}
+\mathcal D_{\mathrm{vis}}
+\mathcal D_{\mathrm{num}}
-W_{\mathrm{ext}}
\approx0.
$$

同时观察损伤局部化、反力随时间的松弛、$q$ 是否接近 1，以及 RLM 代理内部能量是否跟踪物理内部能量。固定载荷稳态松弛只作为该算例长保持时间极限的对照，不再作为主算法。

---

### 23.3 $\alpha$ 参数分析

比较不同固定 $\alpha$：

- 判别式；

- 根的正性；

- $\max|q-1|$；

- 能量误差；

- 条件数；

- 最大稳定时间步。

这里“最大稳定时间步”应与“满足物理精度的最大时间步”分开报告；后者通常更严格。

---

### 23.4 Miehe 与 Amor 对比

在相同网格、材料参数和加载条件下比较：

- 裂纹路径；

- 起裂载荷；

- 纯压缩响应；

- 剪切敏感性；

- 每步 RHS 组装时间；

- 总计算时间；

- $q$ 偏离程度。

Miehe 的全局线性求解成本与 Amor 基本相同，主要额外成本来自正交点特征分解。

另需在相同 $\eta$ 和相同加载历史下比较峰值反力时刻、起裂时刻和保持段松弛曲线。

---

### 23.5 单边缺口拉伸

采用经典 Mode-I 单边缺口拉伸算例，输出：

- 反力—位移曲线；

- 位移—时间、反力—时间和损伤—时间曲线；

- 损伤场；

- 积极弹性能；

- 消极弹性能；

- 裂纹表面能；

- 外部功；

- 累计黏性耗散与 BDF1 数值耗散；

- 原始总势能；

- RLM 代理内部能量；

- $q$ 历史。

至少使用三组 $\Delta t$ 做时间步收敛，并使用多组 $\eta$ 和加载速率区分材料率效应与数值时间步效应。$\eta$ 改变而 $\Delta t$ 收敛后响应仍改变，才是模型的率相关性；仅随 $\Delta t$ 改变的差异属于离散误差。

由于当前没有不可逆性，卸载和复杂变幅载荷结果只能作为算法测试，不能作为严格断裂预测。

---

## 24. 当前模型的局限

> [!warning] 不可逆性缺失
> 当前最大的物理缺陷是允许
>
> $$
> d_t<0.
> $$
>
> 后续需要研究：
>
> - 障碍问题；
>
> - primal–dual active set；
>
> - variational inequality；
>
> - slack variable；
>
> - 单向梯度流
>
>
> 与 RLM 的耦合。

> [!warning] 区间约束缺失
> 对原始连续梯度流，把强形式写成
>
> $$
> d_t
> =
> M G_c\ell\,\Delta d
> +
> 2M(1-\kappa)(1-d)\psi_0^+
> -
> M\frac{G_c}{\ell}d,
> $$
>
> 在足够正则、齐次 Neumann 边界和 $0\leq d(\cdot,0)\leq1$ 下，可以形式上用最大值原理说明区间不变性。但本文的显式冻结 RLM 离散格式并不继承这个最大值原理，因此数值上仍不保证
>
> $$
> 0\leq d\leq1.
> $$
>
> 不应直接裁剪 $d$，除非重新证明裁剪后的能量性质或采用与能量律兼容的界保持离散化。

> [!warning] Miehe 谱切线
> 第一版算法只需要：
>
> $$
> \psi_0^+,
> \qquad
> \boldsymbol\sigma_0^+.
> $$
>
> 若未来使用精确 RLM-Q 曲率或 Newton 法，则必须处理：
>
> - 主应变重根；
>
> - 谱投影导数；
>
> - 主应变过零；
>
> - 一致切线的数值正则化。
>

> [!warning] 准静态范围与惯性缺失
> 时变位移和载荷会向系统输入能量，必须计算反力功并验证功—能平衡。当前位移方程没有 $\rho\mathbf u_{tt}$ 和动能，因此不能描述冲击、弹性波或高速动态断裂；这不是通过减小 $\Delta t$ 就能修复的离散误差，而是模型层面的假设。

> [!warning] 黏性标定
> 把 $t$ 解释为真实时间后，$\eta=1/M$ 必须具有一致物理量纲并接受实验或目标时间尺度标定。未经标定的 $\eta$ 只能用于参数化率相关研究，不能给出可信的绝对起裂时间或裂纹速度。

---

## 25. 推荐研究顺序

$$
\boxed{
\begin{aligned}
&\text{Miehe--RLM-PE--BDF1}
\\
&\Downarrow
\\
&\text{时变载荷下的准静态黏性功--能律与严格适定性}
\\
&\Downarrow
\\
&\text{唯一正根与离散功--能稳定性}
\\
&\Downarrow
\\
&\text{RLM 代理内部能量向物理内部能量收敛}
\\
&\Downarrow
\\
&\text{Miehe--RLM-PE--CN}
\\
&\Downarrow
\\
&\text{黏性参数标定与加载速率验证}
\\
&\Downarrow
\\
&\text{不可逆障碍约束扩展}
\\
&\Downarrow
\\
&\text{含位移惯性的动态断裂扩展}.
\end{aligned}
}
$$

当前最适合作为完整理论和程序原型的是：

$$
\boxed{
\text{Miehe--RLM-PE--BDF1 for time-dependent quasi-static viscous fracture}.
}
$$

在标量方程选到可接受实根、各线性分支精确可解的条件下，它具有：

- 物理上合理的拉压谱分解；

- 常系数全局矩阵；

- 每步四个线性分支；

- 一个标量二次方程；

- 简单的唯一正根充分条件；

- 可直接验证的一步离散功—能恒等式；

- 清晰的原始能量一致性指标。

但“时间一阶收敛、$q\to1$、$\eta$ 的物理标定、不可逆断裂的有效性以及忽略惯性的适用范围”仍是彼此独立、尚需证明或扩展的问题。

---

## 26. 主要参考文献

1. C. Miehe, F. Welschinger, M. Hofacker, “Thermodynamically consistent phase-field models of fracture: Variational principles and multi-field FE implementations,” *International Journal for Numerical Methods in Engineering*, 83 (2010), 1273–1311. [DOI](https://doi.org/10.1002/nme.2861)

2. Q. Cheng, C. Liu, J. Shen, “A new Lagrange multiplier approach for gradient flows,” *Computer Methods in Applied Mechanics and Engineering*, 367 (2020), 113070. [DOI](https://doi.org/10.1016/j.cma.2020.113070)

3. X. Jing, J. Zhao, “Relaxed Lagrange Multiplier (RLM) Schemes for Phase Field Models Preserving the Relaxed Original Energy Dissipation Law,” arXiv:2607.00355v1 (2026). [arXiv](https://arxiv.org/abs/2607.00355)

> [!note] 文献边界
> 第 3 篇是 2026 年的预印本。本文对“准静态弹性平衡 + 损伤梯度流”的推广、四分支求解以及 Miehe 谱分裂的专门结论，应视为基于该框架的新推导，而不是原文已经覆盖的定理。
