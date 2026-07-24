# Miehe 谱分解下的 Mode-I 黏性相场断裂 RLM

本文总结论文草稿 PDF 中的模型与算法。核心是：在准静态位移平衡下，让相场损伤按真实时间进行黏性 $L^2$ 梯度流；用 Miehe 谱分解只退化拉伸能；再用松弛拉格朗日乘子（relaxed Lagrange multiplier, RLM）把每个物理时间步的耦合更新化为四个线性分支加一个标量二次方程。

> [!important]
> 这里的“Mode-I”是加载方式（例如缺口试件的张开型拉伸），不是严格只允许 Mode-I 破坏的本构约束。Miehe 分解抑制纯压缩损伤，但剪切状态仍可能含正主应变并驱动损伤。

## 1. 模型范围与变量

材料区域为 $\Omega\subset\mathbb R^m$，$m=2$ 或 $3$。位移和相场分别为

$$
\mathbf u:\Omega\to\mathbb R^m,
\qquad d:\Omega\to\mathbb R,
$$

约定 $d=0$ 为完整材料，$d=1$ 为完全断裂。边界分解为

$$
\partial\Omega=\Gamma_D\cup\Gamma_N,
\qquad \Gamma_D\cap\Gamma_N=\varnothing.
$$

$$
V_{\mathbf u_D(t)}=\{\mathbf u\in H^1(\Omega;\mathbb R^m):\mathbf u=\mathbf u_D(t)\text{ on }\Gamma_D\},
$$

$$
V_0=\{\mathbf v\in H^1(\Omega;\mathbb R^m):\mathbf v=0\text{ on }\Gamma_D\},
\qquad W=H^1(\Omega).
$$

相场采用自然边界条件

$$
\nabla d\cdot\mathbf n=0\quad\text{on }\partial\Omega.
$$

模型忽略位移惯性，但保留损伤的物理时间尺度。因此它是“准静态力学 + 率相关黏性损伤”，不是含应力波的动态断裂。$\Delta t$ 是物理时间分辨率，不是固定载荷下的伪时间内循环步长；每个物理时间步只推进一次耦合状态。

当前模型不施加

$$
d_t\ge 0,\qquad d^{n+1}\ge d^n,\qquad 0\le d\le1.
$$

所以它是可恢复的黏性损伤模型，卸载时可能愈合，不能直接视为完整的不可逆断裂模型。

## 2. Miehe 谱分解与 AT2 能量

小变形应变为

$$
\boldsymbol\epsilon(\mathbf u)=\frac12(\nabla\mathbf u+\nabla\mathbf u^T).
$$

令 $epsilon_a,\mathbf n_a$ 为主应变和主方向，则

$$
\boldsymbol\epsilon=\sum_{a=1}^m\epsilon_a\,\mathbf n_a\otimes\mathbf n_a.
$$

Macaulay 括号定义为

$$
\langle x\rangle_+=\max(x,0),
\qquad \langle x\rangle_-=\min(x,0).
$$

$$
\boldsymbol\epsilon_+=\sum_a\langle\epsilon_a\rangle_+\mathbf n_a\otimes\mathbf n_a,
\qquad
\boldsymbol\epsilon_-=\sum_a\langle\epsilon_a\rangle_-\mathbf n_a\otimes\mathbf n_a.
$$

因此 $\boldsymbol\epsilon=\boldsymbol\epsilon_++\boldsymbol\epsilon_-$. 对各向同性线弹性材料（Lamé 常数 $\lambda,\mu$），

$$
\psi_0^+(\boldsymbol\epsilon)
=\frac\lambda2\langle\operatorname{tr}\boldsymbol\epsilon\rangle_+^2
+\mu\,\boldsymbol\epsilon_+:\boldsymbol\epsilon_+,
$$

$$
\psi_0^-(\boldsymbol\epsilon)
=\frac\lambda2\langle\operatorname{tr}\boldsymbol\epsilon\rangle_-^2
+\mu\,\boldsymbol\epsilon_-:\boldsymbol\epsilon_-.
$$

有 $\psi_0=\psi_0^++\psi_0^-$，相应应力为

$$
\boldsymbol\sigma_0^+=\lambda\langle\operatorname{tr}\boldsymbol\epsilon\rangle_+\mathbf I+2\mu\boldsymbol\epsilon_+,
$$

$$
\boldsymbol\sigma_0^-=\lambda\langle\operatorname{tr}\boldsymbol\epsilon\rangle_-\mathbf I+2\mu\boldsymbol\epsilon_-.
$$

取残余刚度退化函数

$$
g_\kappa(d)=(1-\kappa)(1-d)^2+\kappa,
\qquad 0<\kappa\ll1,
$$

$$
g_\kappa(0)=1,\qquad g_\kappa(1)=\kappa,\qquad
g_\kappa'(d)=-2(1-\kappa)(1-d).
$$

AT2 裂纹表面能为

$$
\mathcal G_c(d)=\int_\Omega\left(\frac{G_c}{2\ell}d^2+\frac{G_c\ell}{2}|\nabla d|^2\right)\,dx,
$$

其中 $G_c$ 是临界能量释放率，$\ell$ 是正则化长度。

## 3. 原始准静态—黏性系统

体力和 Neumann 载荷定义外力泛函

$$
\ell_{\mathrm{ext}}(t;\mathbf u)=\int_\Omega\mathbf b(t)\cdot\mathbf u\,dx
+\int_{\Gamma_N}\bar{\mathbf t}(t)\cdot\mathbf u\,ds.
$$

内部弹性—断裂能和总势能为

$$
\mathcal H(\mathbf u,d)=\int_\Omega[g_\kappa(d)\psi_0^+(\boldsymbol\epsilon(\mathbf u))+\psi_0^-(\boldsymbol\epsilon(\mathbf u))]\,dx+\mathcal G_c(d),
$$

$$
\mathcal E(t;\mathbf u,d)=\mathcal H(\mathbf u,d)-\ell_{\mathrm{ext}}(t;\mathbf u).
$$

退化应力为

$$
\boldsymbol\sigma(\mathbf u,d)=g_\kappa(d)\boldsymbol\sigma_0^+(\mathbf u)+\boldsymbol\sigma_0^-(\mathbf u).
$$

### 位移平衡

$$
\boxed{
\int_\Omega\boldsymbol\sigma(\mathbf u,d):\boldsymbol\epsilon(\mathbf v)\,dx
=\ell_{\mathrm{ext}}(t;\mathbf v),\quad\forall\mathbf v\in V_0.
}
$$

### 相场化学势与梯度流

化学势 $\mu_d=\delta\mathcal H/\delta d$ 由

$$
(\mu_d,\xi)=\int_\Omega g_\kappa'(d)\psi_0^+(\boldsymbol\epsilon(\mathbf u))\xi\,dx
+\int_\Omega\left(\frac{G_c}{\ell}d\xi+G_c\ell\nabla d\cdot\nabla\xi\right)dx
$$

定义。黏性 $L^2$ 梯度流为

$$
\boxed{\eta(d_t,\xi)+(\mu_d,\xi)=0,\quad\forall\xi\in W,\qquad \eta>0.}
$$

等价地，$\eta d_t=-\mu_d$，迁移率 $M=\eta^{-1}$。$\eta$ 决定真实损伤时间尺度，不能只当作数值加速参数。

若 $\mathbf r_D$ 是一致的 Dirichlet 反力，则外功率为

$$
\mathcal P_{\mathrm{ext}}(t)=\ell_{\mathrm{ext}}(t;\mathbf u_t)+\int_{\Gamma_D}\mathbf r_D\cdot\dot{\mathbf u}_D\,ds.
$$

连续功—能关系为

$$
\boxed{\frac{d}{dt}\mathcal H(\mathbf u,d)=\mathcal P_{\mathrm{ext}}(t)-\eta\|d_t\|_{L^2(\Omega)}^2.}
$$

所以时变加载时应检查“内部能量变化 = 外功 - 黏性耗散”，而不是要求能量单调下降。只有外功率为零时才有单调耗散。

## 4. 线性—非线性能量分裂

定义未损伤弹性和 AT2 双线性形式

$$
a_u(\mathbf u,\mathbf v)=\int_\Omega\boldsymbol\epsilon(\mathbf v):\mathbb C_0:\boldsymbol\epsilon(\mathbf u)\,dx,
$$

$$
a_d(d,\xi)=\int_\Omega\left(\frac{G_c}{\ell}d\xi+G_c\ell\nabla d\cdot\nabla\xi\right)dx.
$$

由于 $\frac12a_u(\mathbf u,\mathbf u)=\int_\Omega(\psi_0^++\psi_0^-)dx$，总势能有精确分裂

$$
\mathcal E=\mathcal E_0+\mathcal E_1,
$$

$$
\mathcal E_0(t;\mathbf u,d)=\frac12a_u(\mathbf u,\mathbf u)+\frac12a_d(d,d)-\ell_{\mathrm{ext}}(t;\mathbf u),
$$

$$
\boxed{\mathcal E_1(\mathbf u,d)=\int_\Omega[g_\kappa(d)-1]\psi_0^+(\boldsymbol\epsilon(\mathbf u))\,dx.}
$$

非线性变分力必须来自同一个 $\mathcal E_1$：

$$
\mathcal N_u(\mathbf u,d;\mathbf v)=\int_\Omega[g_\kappa(d)-1]\boldsymbol\sigma_0^+(\mathbf u):\boldsymbol\epsilon(\mathbf v)\,dx,
$$

$$
\mathcal N_d(\mathbf u,d;\xi)=\int_\Omega g_\kappa'(d)\psi_0^+(\boldsymbol\epsilon(\mathbf u))\xi\,dx.
$$

这样，全球矩阵只含未损伤弹性、AT2 和黏性质量项，不依赖当前 $d$ 或积分点谱分解；Miehe 非线性只进入右端和能量评估。

## 5. 连续 RLM 重构

引入无量纲乘子 $q(t)$ 和固定松弛参数 $\alpha>0$，其中 $\alpha$ 具有能量量纲，通常写成 $\alpha=\widehat\alpha E_{\mathrm{ref}}$。重构系统为

$$
a_u(\mathbf u,\mathbf v)+q\,\mathcal N_u(\mathbf u,d;\mathbf v)=\ell_{\mathrm{ext}}(t;\mathbf v),
$$

$$
\eta(d_t,\xi)+a_d(d,\xi)+q\,\mathcal N_d(\mathbf u,d;\xi)=0.
$$

标量 RLM 关系为

$$
\boxed{
\frac{d}{dt}\mathcal E_1(\mathbf u,d)+\alpha\frac{d}{dt}q^2
=q[\mathcal N_u(\mathbf u,d;\mathbf u_t)+\mathcal N_d(\mathbf u,d;d_t)].
}
$$

定义松弛原始内部能量

$$
\boxed{\widetilde{\mathcal H}(\mathbf u,d,q)=\mathcal H(\mathbf u,d)+\alpha(q^2-1).}
$$

形式上有

$$
\boxed{\frac{d}{dt}\widetilde{\mathcal H}=\mathcal P_{\mathrm{ext}}-\eta\|d_t\|_{L^2}^2.}
$$

连续分支 $q\equiv1$ 恢复原始 Miehe 模型。时变非齐次 Dirichlet 数据须通过提升 $\mathbf u=\mathbf w+\mathbf z(t)$ 处理，否则会漏掉反力功。上述耦合系统的专门收敛定理仍需单独证明，不能直接当作一般 RLM 梯度流定理的现成推论。

## 6. 首选离散格式：BDF1–QM

PDF 的首选格式是带方向二次上界（quadratic majorant, QM）的 BDF1。它保留常系数线性系统、一步能量估计和标量二次闭合。令 $t_{n+1}=t_n+\Delta t$，并在 $U^n=(\mathbf u^n,d^n)$ 处冻结非线性力：

$$
\mathcal N_u^n(\mathbf v)=\mathcal N_u(\mathbf u^n,d^n;\mathbf v),
\qquad
\mathcal N_d^n(\xi)=\mathcal N_d(\mathbf u^n,d^n;\xi).
$$

单步方程为

$$
a_u(\mathbf u^{n+1},\mathbf v)+q^{n+1}\mathcal N_u^n(\mathbf v)=\ell_{\mathrm{ext}}(t_{n+1};\mathbf v),
$$

$$
\frac\eta{\Delta t}(d^{n+1}-d^n,\xi)+a_d(d^{n+1},\xi)+q^{n+1}\mathcal N_d^n(\xi)=0.
$$

### 6.1 仿射分支

写成

$$
\mathbf u^{n+1}=\mathbf u_a^{n+1}+q^{n+1}\mathbf u_b^{n+1},
\qquad
d^{n+1}=d_a^{n+1}+q^{n+1}d_b^{n+1}.
$$

四个线性问题是

$$
a_u(\mathbf u_a^{n+1},\mathbf v)=\ell_{\mathrm{ext}}(t_{n+1};\mathbf v),
\qquad \mathbf u_a^{n+1}=\mathbf u_D(t_{n+1})\text{ on }\Gamma_D,
$$

$$
a_u(\mathbf u_b^{n+1},\mathbf v)=-\mathcal N_u^n(\mathbf v),
\qquad \mathbf u_b^{n+1}=0\text{ on }\Gamma_D,
$$

$$
\frac\eta{\Delta t}(d_a^{n+1}-d^n,\xi)+a_d(d_a^{n+1},\xi)=0,
$$

$$
\frac\eta{\Delta t}(d_b^{n+1},\xi)+a_d(d_b^{n+1},\xi)=-\mathcal N_d^n(\xi).
$$

位移矩阵 $K_u\leftrightarrow a_u$，相场矩阵

$$
K_d\leftrightarrow \frac\eta{\Delta t}(\cdot,\cdot)+a_d(\cdot,\cdot).
$$

固定 $\Delta t$ 时二者可预分解并复用；改变 $\Delta t$ 通常需要重新分解 $K_d$。时变载荷只改变基准分支右端。

### 6.2 方向二次上界

沿仿射方向定义

$$
\phi_n(q)=\mathcal E_1(U_a^{n+1}+qU_b^{n+1}),
$$

$$
\phi_n^*=\phi_n(1),
\qquad g_n^*=\phi_n'(1)=D\mathcal E_1(U_*^{n+1})[U_b^{n+1}],
$$

其中 $U_*^{n+1}=U_a^{n+1}+U_b^{n+1}$ 是 $q=1$ 预测态。给定包含候选根的区间 $I_n$，取 $S_n\ge0$ 并定义

$$
\Pi_n(q)=\phi_n^*+g_n^*(q-1)+\frac{S_n}{2}(q-1)^2.
$$

QM 的必要验收条件是

$$
\boxed{\phi_n(q)\le\Pi_n(q),\qquad q\in I_n.}
$$

仅有 $S_n\ge0$ 不够；若 $\phi_n'$ 在 $I_n$ 上 Lipschitz 常数为 $L_n$，则 $S_n\ge L_n$ 可保证上界。实现上可用方向有限差分估计初值，再在候选根处直接检查上界；不满足时增大 $S_n$，只需重解标量方程，不需重做四个线性分支。

### 6.3 标量二次方程

用 $\Pi_n$ 替代 $\phi_n$ 的标量闭合为

$$
\Pi_n(q^{n+1})-\mathcal E_1(U^n)+\alpha[(q^{n+1})^2-(q^n)^2]
=q^{n+1}[\mathcal N_u^n(\mathbf u^{n+1}-\mathbf u^n)+\mathcal N_d^n(d^{n+1}-d^n)].
$$

定义

$$
c_0^n=\mathcal N_u^n(\mathbf u_a^{n+1}-\mathbf u^n)+\mathcal N_d^n(d_a^{n+1}-d^n),
$$

$$
c_1^n=\mathcal N_u^n(\mathbf u_b^{n+1})+\mathcal N_d^n(d_b^{n+1}).
$$

则 $q^{n+1}$ 满足

$$
\boxed{A_n(q^{n+1})^2+B_nq^{n+1}+C_n=0,}
$$

其中

$$
A_n=\alpha+\frac{S_n}{2}-c_1^n,
\qquad B_n=g_n^*-S_n-c_0^n,
$$

$$
\boxed{C_n=\phi_n^*-g_n^*+\frac{S_n}{2}-\mathcal E_1(U^n)-\alpha(q^n)^2.}
$$

对于响应分支，有

$$
c_1^n=-a_u(\mathbf u_b^{n+1},\mathbf u_b^{n+1})
-\frac\eta{\Delta t}\|d_b^{n+1}\|_{L^2}^2-a_d(d_b^{n+1},d_b^{n+1})\le0,
$$

故

$$
A_n\ge\alpha+\frac{S_n}{2}>0.
$$

判别式为 $D_n=B_n^2-4A_nC_n$。当 $q^n\ne0$ 时，易检查的唯一正根充分条件是

$$
C_n<0
\quad\Longleftarrow\quad
\boxed{
\alpha>
\max\left\{0,\frac{\phi_n^*-g_n^*+S_n/2-\mathcal E_1(U^n)}{(q^n)^2}\right\}.
}
$$

这只是充分条件，不是必要条件；实际仍须检查 $D_n\ge0$、根是否属于 $I_n$、QM 上界以及原二次方程残差。若回溯增大 $S_n$ 后仍无可接受根，应减小 $\Delta t$，不能把不满足上界的根当成稳定根。

二次方程宜用抗消去公式求解：

$$
\widehat q=-\frac12\left(B_n+\operatorname{copysign}(\sqrt{D_n},B_n)\right),
\qquad q_1=\frac{\widehat q}{A_n},\quad q_2=\frac{C_n}{\widehat q}.
$$

数值上只在 $D_n\in[-\varepsilon_D,0)$ 时截断为零；$D_n< -\varepsilon_D$ 应报告无实根。删除非有限根和多项式残差过大的根；有唯一正根保证时选正根，只有存在多个可接受根时才用最接近 $q^n$ 的规则。

## 7. BDF1 功—能不等式

定义离散松弛原始内部能量

$$
\boxed{\widetilde{\mathcal H}_{\mathrm{QM}}^n=\mathcal H(U^n)+\alpha[(q^n)^2-1].}
$$

令 $\Delta\mathbf u^n=\mathbf u^{n+1}-\mathbf u^n$，$\Delta d^n=d^{n+1}-d^n$。取 $\mathbf z$ 为 Dirichlet 边界增量的任意有限元提升，定义

$$
\mathcal R_D^{n+1}(\mathbf z)=a_u(\mathbf u^{n+1},\mathbf z)+q^{n+1}\mathcal N_u^n(\mathbf z)-\ell_{\mathrm{ext}}(t_{n+1};\mathbf z),
$$

$$
W_{\mathrm{ext}}^{n+1}=\ell_{\mathrm{ext}}(t_{n+1};\Delta\mathbf u^n)+\mathcal R_D^{n+1}(\mathbf z).
$$

在线性分支求解充分准确、所选根属于 $I_n$ 且满足 QM 上界、外功离散一致时，

$$
\boxed{
\begin{aligned}
\widetilde{\mathcal H}_{\mathrm{QM}}^{n+1}-\widetilde{\mathcal H}_{\mathrm{QM}}^n
\le{}&W_{\mathrm{ext}}^{n+1}
-\frac12a_u(\Delta\mathbf u^n,\Delta\mathbf u^n)
\\
&-\frac12a_d(\Delta d^n,\Delta d^n)
-\frac\eta{\Delta t}\|\Delta d^n\|_{L^2}^2.
\end{aligned}}
$$

最后一项是该步黏性耗散，前两项是 BDF1 数值耗散。该估计对 $\Delta t$ 无条件稳定，但“无条件稳定”不等于任意大时间步都能解析真实加载；它也不保证标量方程一定存在可接受根。只有 $W_{\mathrm{ext}}^{n+1}=0$ 时，才得到能量单调下降。

与原始物理能量的差异在 QM 方案中是显式的：

$$
\widetilde{\mathcal H}_{\mathrm{QM}}^n-\mathcal H(U^n)=\alpha[(q^n)^2-1].
$$

因此应持续监测 $|q^n-1|$，而不必再引入 PE-BDF1 中的 $P^n-\mathcal E_1(U^n)$ 作为稳定性误差。

## 9. 参数、适用范围与验证重点

- $\eta$ 是材料/正则化参数，应由松弛或加载速率试验标定；$\Delta t$ 是数值分辨率，应做时间步收敛；$\lambda(t)$ 是物理加载历史。三者不能混为一谈。
- 即使能量估计对时间步无条件稳定，过大的 $\Delta t$ 仍可能错过起裂、峰值反力和快速损伤阶段；可按位移增量、损伤增量或功—能残差自适应减小时间步。
- 准静态假设要求加载和损伤时间尺度慢于弹性波传播/结构振动时间尺度。冲击、高速断裂需要加入 $\rho\mathbf u_{tt}$ 和动能，并重新推导 RLM。
- Miehe 谱能量通常只有 $C^1$；主应变过零或重根附近的二阶谱切线需要极限/除差公式。首版 QM 实现只需计算 $\psi_0^+$、$\boldsymbol\sigma_0^+$ 和方向能量，不需要完整谱 Hessian。
- 验证应至少覆盖：制造解时间收敛、升载—保持过程的功—能闭合、$\alpha$ 敏感性、Miehe 与 Amor 分解对比，以及单边缺口 Mode-I 拉伸试验。应区分 $\eta$ 引起的真实率效应和 $\Delta t$ 引起的离散误差。

## 10. 主要局限与后续方向

当前没有 $d_t\ge0$、$0\le d\le1$ 的显式约束，也没有惯性项。若要得到不可逆断裂，需要把 RLM 与障碍问题、变分不等式、primal–dual active set、slack variable 或单向梯度流相容地耦合；不能依赖事后裁剪。

## 11. 参考文献

1. C. Miehe, F. Welschinger, M. Hofacker, “Thermodynamically consistent phase-field models of fracture: Variational principles and multi-field FE implementations,” *International Journal for Numerical Methods in Engineering*, 83 (2010), 1273–1311. [DOI](https://doi.org/10.1002/nme.2861)
2. Q. Cheng, C. Liu, J. Shen, “A new Lagrange multiplier approach for gradient flows,” *Computer Methods in Applied Mechanics and Engineering*, 367 (2020), 113070. [DOI](https://doi.org/10.1016/j.cma.2020.113070)
3. X. Jing, J. Zhao, “Relaxed Lagrange Multiplier (RLM) Schemes for Phase Field Models Preserving the Relaxed Original Energy Dissipation Law,” arXiv:2607.00355 (2026). [arXiv](https://arxiv.org/abs/2607.00355)

> [!note]
> PDF 对“准静态位移平衡 + 黏性相场梯度流 + Miehe 谱分解”的耦合 RLM 结论是本文档所依据的专门构造；它不应被误解为上述 RLM 文献已经直接证明的完整耦合收敛定理。
