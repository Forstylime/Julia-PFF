# Miehe 谱分解下的 Mode-I 黏性相场断裂 RLM–BDF2

本文在 [[RLM-BDF1_Mode_I_Miehe]] 的基础上构造二阶时间离散。连续模型、Miehe 拉压谱分解、AT2 裂纹表面能以及线性—非线性能量分裂均保持不变；主要改动是把相场时间导数改为 BDF2、把非线性力改为二阶外推，并把一步能量估计改写为两步 $G$-稳定估计。

本文给出两种共享同一组四个线性分支的标量闭合。首选方案直接代入原始非线性能量 $\mathcal E_1(U_a+qU_b)$，得到一维非线性方程并用 Brent 等标量算法求根，从而不引入 QM 代理缺口；备选方案使用方向二次上界，将标量问题化为一元二次方程。

> [!important]
> 这里的“原始能量直接闭合”是指标量方程逐次计算真实的 $\mathcal E_1$，而不是方向上界 $\Pi_n$。BDF2 稳定量仍是含两层历史的 $G$-能量，并非单时刻的 $\mathcal H(U^n)$。这里的“二阶”首先指固定时间步下的形式二阶一致性；主应变过零或发生重根时，Miehe 谱分解的低正则性可能降低实际时间收敛阶。

## 1. 保持不变的连续模型

材料区域为 $\Omega\subset\mathbb R^m$，其中 $m=2$ 或 $3$。位移为 $\mathbf u:\Omega\to\mathbb R^m$，相场为 $d:\Omega\to\mathbb R$，并约定 $d=0$ 表示完整材料、$d=1$ 表示完全断裂。位移试探与测试空间为$$
V_{\mathbf u_D(t)}=\{\mathbf u\in H^1(\Omega;\mathbb R^m):\mathbf u=\mathbf u_D(t)\text{ on }\Gamma_D\},
\qquad
V_0=\{\mathbf v\in H^1(\Omega;\mathbb R^m):\mathbf v=0\text{ on }\Gamma_D\},
$$
相场空间为 $W=H^1(\Omega)$，自然边界条件为 $\nabla d\cdot\mathbf n=0$。

小变形应变及其谱分解为$$
\boldsymbol\epsilon(\mathbf u)=\frac12(\nabla\mathbf u+\nabla\mathbf u^T),
\qquad
\boldsymbol\epsilon=\sum_{a=1}^m\epsilon_a\,\mathbf n_a\otimes\mathbf n_a,
$$
$$
\boldsymbol\epsilon_\pm
=\sum_{a=1}^m\langle\epsilon_a\rangle_\pm\mathbf n_a\otimes\mathbf n_a,
\qquad
\langle x\rangle_+=\max(x,0),
\qquad
\langle x\rangle_-=\min(x,0).
$$
Miehe 拉伸与压缩能量密度为$$
\psi_0^+(\boldsymbol\epsilon)
=\frac{\lambda}{2}\langle\operatorname{tr}\boldsymbol\epsilon\rangle_+^2
+\mu\,\boldsymbol\epsilon_+:\boldsymbol\epsilon_+,
$$
$$
\psi_0^-(\boldsymbol\epsilon)
=\frac{\lambda}{2}\langle\operatorname{tr}\boldsymbol\epsilon\rangle_-^2
+\mu\,\boldsymbol\epsilon_-:\boldsymbol\epsilon_-.
$$
相应的拉伸与压缩应力为$$
\boldsymbol\sigma_0^+
=\lambda\langle\operatorname{tr}\boldsymbol\epsilon\rangle_+\mathbf I
+2\mu\boldsymbol\epsilon_+,
\qquad
\boldsymbol\sigma_0^-
=\lambda\langle\operatorname{tr}\boldsymbol\epsilon\rangle_-\mathbf I
+2\mu\boldsymbol\epsilon_-.
$$

取残余刚度退化函数$$
g_\kappa(d)=(1-\kappa)(1-d)^2+\kappa,
\qquad
g_\kappa'(d)=-2(1-\kappa)(1-d),
\qquad
0<\kappa\ll1.
$$
内部能量为$$
\mathcal H(\mathbf u,d)
=\int_\Omega
\left[g_\kappa(d)\psi_0^+(\boldsymbol\epsilon(\mathbf u))
+\psi_0^-(\boldsymbol\epsilon(\mathbf u))\right]dx
+\frac12a_d(d,d),
$$
其中 AT2 双线性形式为$$
a_d(d,\xi)
=\int_\Omega\left(\frac{G_c}{\ell}d\xi+G_c\ell\nabla d\cdot\nabla\xi\right)dx.
$$
外力泛函和总势能分别为$$
\ell_{\mathrm{ext}}(t;\mathbf u)
=\int_\Omega\mathbf b(t)\cdot\mathbf u\,dx
+\int_{\Gamma_N}\bar{\mathbf t}(t)\cdot\mathbf u\,ds,
\qquad
\mathcal E(t;\mathbf u,d)=\mathcal H(\mathbf u,d)-\ell_{\mathrm{ext}}(t;\mathbf u).
$$

与 BDF1 笔记相同，当前模型没有显式施加 $d_t\ge0$、$d^{n+1}\ge d^n$ 或 $0\le d\le1$，因此仍是允许卸载愈合的黏性损伤模型。

## 2. 不变的能量分裂与连续 RLM

定义未损伤弹性双线性形式$$
a_u(\mathbf u,\mathbf v)
=\int_\Omega
\boldsymbol\epsilon(\mathbf v):\mathbb C_0:\boldsymbol\epsilon(\mathbf u)\,dx.
$$
内部能量具有精确分裂$$
\mathcal H(U)
=\frac12a_u(\mathbf u,\mathbf u)
+\frac12a_d(d,d)
+\mathcal E_1(U),
\qquad
U=(\mathbf u,d),
$$
其中非线性能量为$$
\mathcal E_1(\mathbf u,d)
=\int_\Omega
[g_\kappa(d)-1]\psi_0^+(\boldsymbol\epsilon(\mathbf u))\,dx.
$$
来自同一个 $\mathcal E_1$ 的非线性变分力为$$
\mathcal N_u(U;\mathbf v)
=\int_\Omega
[g_\kappa(d)-1]\boldsymbol\sigma_0^+(\mathbf u):
\boldsymbol\epsilon(\mathbf v)\,dx,
$$
$$
\mathcal N_d(U;\xi)
=\int_\Omega
g_\kappa'(d)\psi_0^+(\boldsymbol\epsilon(\mathbf u))\xi\,dx.
$$

引入无量纲标量乘子 $q(t)$ 和具有能量量纲的参数 $\alpha>0$。连续 RLM 系统为$$
a_u(\mathbf u,\mathbf v)+q\,\mathcal N_u(U;\mathbf v)
=\ell_{\mathrm{ext}}(t;\mathbf v),
\qquad
\forall\mathbf v\in V_0,
$$
$$
\eta(d_t,\xi)+a_d(d,\xi)+q\,\mathcal N_d(U;\xi)=0,
\qquad
\forall\xi\in W,
$$
$$
\frac{d}{dt}\mathcal E_1(U)
+\alpha\frac{d}{dt}q^2
=q\left[\mathcal N_u(U;\mathbf u_t)+\mathcal N_d(U;d_t)\right].
$$
连续一致分支为 $q\equiv1$。相应的松弛内部能量为$$
\widetilde{\mathcal H}(U,q)
=\mathcal H(U)+\alpha(q^2-1).
$$

## 3. 固定步长 BDF2 离散

令 $t_n=n\Delta t$，并对任意离散序列 $x^n$ 定义$$
\delta_Bx^{n+1}=3x^{n+1}-4x^n+x^{n-1},
\qquad
D_2x^{n+1}=\frac{\delta_Bx^{n+1}}{2\Delta t},
$$
以及二阶显式外推$$
\overline x^{\,n+1}=2x^n-x^{n-1},
\qquad
\overline U^{\,n+1}=2U^n-U^{n-1}.
$$
在外推状态冻结非线性力：$$
\overline{\mathcal N}_u^{\,n+1}(\mathbf v)
=\mathcal N_u(\overline U^{\,n+1};\mathbf v),
\qquad
\overline{\mathcal N}_d^{\,n+1}(\xi)
=\mathcal N_d(\overline U^{\,n+1};\xi).
$$

对 $n\ge1$，BDF2 场方程为$$
a_u(\mathbf u^{n+1},\mathbf v)
+q^{n+1}\overline{\mathcal N}_u^{\,n+1}(\mathbf v)
=\ell_{\mathrm{ext}}(t_{n+1};\mathbf v),
\qquad
\forall\mathbf v\in V_0,
$$
$$
\frac{\eta}{2\Delta t}(\delta_Bd^{n+1},\xi)
+a_d(d^{n+1},\xi)
+q^{n+1}\overline{\mathcal N}_d^{\,n+1}(\xi)=0,
\qquad
\forall\xi\in W.
$$
对应的 BDF2 标量闭合先写成$$
\begin{aligned}
&3\mathcal E_1(U^{n+1})-4\mathcal E_1(U^n)+\mathcal E_1(U^{n-1})\\
&\quad+\alpha\left[3(q^{n+1})^2-4(q^n)^2+(q^{n-1})^2\right]\\
&=q^{n+1}\left[
\overline{\mathcal N}_u^{\,n+1}(\delta_B\mathbf u^{n+1})
+\overline{\mathcal N}_d^{\,n+1}(\delta_Bd^{n+1})
\right].
\end{aligned}
$$
对于 Miehe 谱能量，直接使用 $\mathcal E_1(U^{n+1})$ 会产生非多项式或分段光滑的一维标量方程，但未知量仍然只有一个。第 6 节直接对该方程求根；第 5 节保留方向二次上界作为计算更便宜的备选闭合。

## 4. 四个常系数线性分支

写成关于 $q^{n+1}$ 的仿射分解$$
\mathbf u^{n+1}
=\mathbf u_a^{n+1}+q^{n+1}\mathbf u_b^{n+1},
\qquad
d^{n+1}
=d_a^{n+1}+q^{n+1}d_b^{n+1}.
$$
位移基准分支满足$$
a_u(\mathbf u_a^{n+1},\mathbf v)
=\ell_{\mathrm{ext}}(t_{n+1};\mathbf v),
\qquad
\mathbf u_a^{n+1}=\mathbf u_D(t_{n+1})\text{ on }\Gamma_D,
$$
位移响应分支满足$$
a_u(\mathbf u_b^{n+1},\mathbf v)
=-\overline{\mathcal N}_u^{\,n+1}(\mathbf v),
\qquad
\mathbf u_b^{n+1}=0\text{ on }\Gamma_D.
$$
相场基准分支满足$$
\frac{3\eta}{2\Delta t}(d_a^{n+1},\xi)
+a_d(d_a^{n+1},\xi)
=\frac{\eta}{2\Delta t}(4d^n-d^{n-1},\xi),
$$
相场响应分支满足$$
\frac{3\eta}{2\Delta t}(d_b^{n+1},\xi)
+a_d(d_b^{n+1},\xi)
=-\overline{\mathcal N}_d^{\,n+1}(\xi).
$$

固定 $\Delta t$ 时，位移矩阵仍对应 $a_u$，相场矩阵对应$$
K_d^{\mathrm{BDF2}}
\longleftrightarrow
\frac{3\eta}{2\Delta t}(\cdot,\cdot)+a_d(\cdot,\cdot).
$$
两者均可预分解并在全部 BDF2 时间步复用。与 BDF1 相比，只需更换相场质量项系数与基准分支历史右端；四分支结构不变。

## 5. 共同方向量与备选 QM 二次闭合

先定义两条标量闭合路线共用的方向能量：$$
\phi_n(q)
=\mathcal E_1(U_a^{n+1}+qU_b^{n+1}),
\qquad
U_a^{n+1}=(\mathbf u_a^{n+1},d_a^{n+1}),
\qquad
U_b^{n+1}=(\mathbf u_b^{n+1},d_b^{n+1}).
$$
以连续一致点 $q=1$ 为中心，记$$
\phi_n^*=\phi_n(1),
\qquad
g_n^*=\phi_n'(1)
=D\mathcal E_1(U_a^{n+1}+U_b^{n+1})[U_b^{n+1}].
$$
令 $U_*^{n+1}=U_a^{n+1}+U_b^{n+1}$，则同一个方向导数也可由变分力计算：$$
g_n^*
=\mathcal N_u(U_*^{n+1};\mathbf u_b^{n+1})
+\mathcal N_d(U_*^{n+1};d_b^{n+1}).
$$
给定包含候选根的搜索区间 $I_n$，取 $S_n\ge0$ 并定义方向二次代理$$
\Pi_n(q)
=\phi_n^*+g_n^*(q-1)+\frac{S_n}{2}(q-1)^2.
$$
QM 验收条件是$$
\phi_n(q)\le\Pi_n(q),
\qquad
q\in I_n.
$$
实际实现至少要在最终候选根处直接检查该不等式。若希望得到区间保证，则需让 $S_n$ 控制 $\phi_n'$ 在 $I_n$ 上的 Lipschitz 常数。Miehe 谱分解在主应变过零或重根附近只有有限光滑性，因此方向差分估计只能用于给出 $S_n$ 初值，不能代替最终上界检查。

用 $\Pi_n(q^{n+1})$ 替换标量闭合中的 $\mathcal E_1(U^{n+1})$，得到$$
\begin{aligned}
&3\Pi_n(q^{n+1})-4\mathcal E_1(U^n)+\mathcal E_1(U^{n-1})\\
&\quad+\alpha\left[3(q^{n+1})^2-4(q^n)^2+(q^{n-1})^2\right]\\
&=q^{n+1}\left[
\overline{\mathcal N}_u^{\,n+1}(\delta_B\mathbf u^{n+1})
+\overline{\mathcal N}_d^{\,n+1}(\delta_Bd^{n+1})
\right].
\end{aligned}
$$
定义两种闭合路线共用的两个已知标量：$$
\begin{aligned}
c_0^n={}&
\overline{\mathcal N}_u^{\,n+1}
(3\mathbf u_a^{n+1}-4\mathbf u^n+\mathbf u^{n-1})\\
&+\overline{\mathcal N}_d^{\,n+1}
(3d_a^{n+1}-4d^n+d^{n-1}),
\end{aligned}
$$
$$
c_1^n
=\overline{\mathcal N}_u^{\,n+1}(\mathbf u_b^{n+1})
+\overline{\mathcal N}_d^{\,n+1}(d_b^{n+1}).
$$
于是右端为 $q^{n+1}[c_0^n+3q^{n+1}c_1^n]$，而 $q^{n+1}$ 满足一元二次方程$$
A_n(q^{n+1})^2+B_nq^{n+1}+C_n=0,
$$
其中$$
A_n
=3\left(\alpha+\frac{S_n}{2}-c_1^n\right),
\qquad
B_n=3(g_n^*-S_n)-c_0^n,
$$
$$
\begin{aligned}
C_n={}&
3\left(\phi_n^*-g_n^*+\frac{S_n}{2}\right)
-4\mathcal E_1(U^n)+\mathcal E_1(U^{n-1})\\
&-4\alpha(q^n)^2+\alpha(q^{n-1})^2.
\end{aligned}
$$

由两个响应分支可得$$
\begin{aligned}
c_1^n={}&
-a_u(\mathbf u_b^{n+1},\mathbf u_b^{n+1})
-\frac{3\eta}{2\Delta t}\|d_b^{n+1}\|_{L^2(\Omega)}^2\\
&-a_d(d_b^{n+1},d_b^{n+1})
\le0,
\end{aligned}
$$
因此 $A_n\ge3(\alpha+S_n/2)>0$。若$$
\beta_n=4(q^n)^2-(q^{n-1})^2>0
$$
并且$$
\alpha>
\max\left\{
0,\,
\frac{
3(\phi_n^*-g_n^*+S_n/2)
-4\mathcal E_1(U^n)+\mathcal E_1(U^{n-1})
}{\beta_n}
\right\},
$$
则 $C_n<0$，从而二次方程恰有一个正根。这只是便于检查的充分条件；当 $\beta_n\le0$ 或条件不满足时，仍可直接检查判别式和全部候选根。

令 $D_n=B_n^2-4A_nC_n$。为避免二次公式的消去误差，可计算$$
\widehat q
=-\frac12\left[
B_n+\operatorname{copysign}(\sqrt{D_n},B_n)
\right],
\qquad
q_1=\frac{\widehat q}{A_n},
\qquad
q_2=\frac{C_n}{\widehat q}.
$$
只在 $D_n\in[-\varepsilon_D,0)$ 时才可把判别式截断为零。每个候选根都必须通过有限性、二次多项式残差、$q\in I_n$ 和 QM 上界检查。若有多个可接受根，应连续跟踪 $q\equiv1$ 的物理解支；可用到二阶预测值 $2q^n-q^{n-1}$ 的距离作为选择准则，同时持续监测 $|q^{n+1}-1|$。

## 6. 首选方案：原始非线性能量的一维直接闭合

如果不使用 $\Pi_n$，则直接保留$$
\phi_n(q)=\mathcal E_1(U_a^{n+1}+qU_b^{n+1}).
$$
沿用第 5 节定义的 $c_0^n$ 和 $c_1^n$，精确标量闭合可写成$$
F_n(q^{n+1})=0,
$$
其中一维非线性残差为$$
\begin{aligned}
F_n(q)={}&
3\phi_n(q)
+3(\alpha-c_1^n)q^2
-c_0^nq\\
&-4\mathcal E_1(U^n)+\mathcal E_1(U^{n-1})\\
&-4\alpha(q^n)^2+\alpha(q^{n-1})^2.
\end{aligned}
$$
该式只是把仿射表示 $U^{n+1}(q)=U_a^{n+1}+qU_b^{n+1}$ 代回原始 BDF2 标量关系，没有对 $\mathcal E_1$ 作上界、插值或多项式拟合。由于 $c_1^n\le0$，显式二次部分满足 $\alpha-c_1^n>0$；但 $\phi_n(q)$ 仍可能非凸，所以这并不自动保证全局唯一根。

若需要导数，方向导数为$$
\phi_n'(q)
=\mathcal N_u(U_a^{n+1}+qU_b^{n+1};\mathbf u_b^{n+1})
+\mathcal N_d(U_a^{n+1}+qU_b^{n+1};d_b^{n+1}),
$$
从而$$
F_n'(q)
=3\phi_n'(q)+6(\alpha-c_1^n)q-c_0^n.
$$
Miehe 能量通常是 $C^1$ 而不是全局 $C^2$，所以 $F_n'$ 可用于带保护的 Newton 或 Newton–Brent 混合迭代，但不宜依赖未经处理的二阶导数。纯 Brent 法只需函数值，通常更稳健。

### 6.1 Brent 括区与根分支选择

取二阶预测中心$$
q_{\mathrm{pred}}^{n+1}=2q^n-q^{n-1}.
$$
一个可直接实现的括区流程如下：

1. 先给定物理解搜索区间 $I_n=[q_{\min},q_{\max}]$，通常取 $q_{\min}>0$，并把 $q_{\mathrm{pred}}^{n+1}$ 与 $1$ 都纳入候选扫描点。先计算 $F_n(1)$；若其归一化残差已满足容差，可直接接受 $q^{n+1}=1$。
2. 否则把 $q_{\mathrm{pred}}^{n+1}$ 投影到 $I_n$，以该点为中心设置小半径 $r_0$，计算局部区间端点及中心处的 $F_n$。
3. 若相邻点出现异号，则得到 Brent 所需的括区；否则按固定倍率扩大半径，直到找到异号区间或到达 $I_n$ 边界。
4. 为避免漏掉同一区间内的多个根，可在 $I_n$ 上先作很便宜的一维粗扫描，对每个异号子区间分别调用 Brent。
5. 对全部收敛根检查 $|F_n(q)|$、括区宽度、有限性和 $q\in I_n$。若有多个可接受根，选择最接近 $q_{\mathrm{pred}}^{n+1}$ 且与 $q\equiv1$ 连续的正根。

偶重根不产生符号变化，单纯 Brent 扫描可能漏掉它。若粗扫描发现 $|F_n|$ 有接近零的局部极小值但没有异号区间，可从该点启动带区间保护的 Newton 或割线法，并仍以最终残差作为验收条件。若搜索区间内没有合格根，应减小 $\Delta t$ 后重做外推和四个线性分支；不能用“最小化 $|F_n|$”得到的非零残差点代替根。

### 6.2 每次函数计算的实际代价

一次 $F_n(q)$ 计算只需要形成$$
d(q)=d_a^{n+1}+qd_b^{n+1},
\qquad
\boldsymbol\epsilon(q)
=\boldsymbol\epsilon(\mathbf u_a^{n+1})
+q\boldsymbol\epsilon(\mathbf u_b^{n+1}),
$$
并在积分点计算$$
\phi_n(q)
=\int_\Omega
[g_\kappa(d(q))-1]\psi_0^+(\boldsymbol\epsilon(q))\,dx.
$$
因此每次标量迭代只进行一次全域能量求积和积分点谱分解，不需要重新组装或求解任何全局线性系统。实现时可缓存 $d_a^{n+1}$、$d_b^{n+1}$、$\boldsymbol\epsilon(\mathbf u_a^{n+1})$、$\boldsymbol\epsilon(\mathbf u_b^{n+1})$ 以及求积权重，使一次函数计算退化为积分点循环和一次全局求和。

场方程中的 $\overline{\mathcal N}^{\,n+1}$、历史能量 $\mathcal E_1(U^n)$、$\mathcal E_1(U^{n-1})$ 与标量迭代中的 $\phi_n(q)$ 必须采用完全相同的网格、积分点和求积规则。否则即使 Brent 把代码中的 $F_n$ 解到机器精度，离散功—能恒等式仍会出现系统误差。

### 6.3 推荐的验收与回退

标量残差宜按能量尺度归一化，例如要求$$
\frac{|F_n(q^{n+1})|}
{E_{\mathrm{scale}}^n}
\le\tau_F,
$$
其中 $E_{\mathrm{scale}}^n$ 至少覆盖 $|\phi_n(q^{n+1})|$、$|\mathcal E_1(U^n)|$、$|\mathcal E_1(U^{n-1})|$ 和 $\alpha$ 相关项。还应同时检查 $|q^{n+1}-1|$、与预测根的距离以及下一节的精确 $G$-能量平衡残差。

若频繁无法括住根，依次检查搜索区间是否过窄、时间步是否过大、积分是否一致以及 $\alpha$ 是否过小。$\alpha$ 应在一次完整计算中保持固定；若要改变 $\alpha$，应从初始时刻重新计算，不能在中途直接替换而继续沿用旧的 $q$ 历史和能量历史。

## 7. BDF2 的 $G$-能量与功—能关系

定义线性能量范数$$
\|U\|_{\mathcal L}^2
=a_u(\mathbf u,\mathbf u)+a_d(d,d).
$$
BDF2 极化恒等式为$$
\begin{aligned}
a(x^{n+1},\delta_Bx^{n+1})
=\frac12\big[
&\|x^{n+1}\|_a^2
+\|2x^{n+1}-x^n\|_a^2\\
&-\|x^n\|_a^2
-\|2x^n-x^{n-1}\|_a^2\\
&+\|x^{n+1}-2x^n+x^{n-1}\|_a^2
\big].
\end{aligned}
$$
据此定义两步松弛 $G$-能量$$
\begin{aligned}
\mathscr G_{\mathrm{RLM}}^n={}&
\frac14\left[
\|U^n\|_{\mathcal L}^2
+\|2U^n-U^{n-1}\|_{\mathcal L}^2
\right]\\
&+\frac32\mathcal E_1(U^n)
-\frac12\mathcal E_1(U^{n-1})\\
&+\alpha\left[
\frac32(q^n)^2-\frac12(q^{n-1})^2-1
\right].
\end{aligned}
$$
该量包含原始非线性能量的精确历史值，但它是 BDF2 的两步修正能量，不等于单时刻物理内部能量 $\mathcal H(U^n)$。

它与单时刻松弛内部能量的差为$$
\begin{aligned}
\mathscr G_{\mathrm{RLM}}^n-\widetilde{\mathcal H}(U^n,q^n)
&=\frac14\left[
\|2U^n-U^{n-1}\|_{\mathcal L}^2-\|U^n\|_{\mathcal L}^2
\right]\\
&+\frac12\left[
\mathcal E_1(U^n)-\mathcal E_1(U^{n-1})
\right]\\
&+\frac{\alpha}{2}\left[
(q^n)^2-(q^{n-1})^2
\right].
\end{aligned}
$$
因此在有限时间步下，$G$-能量单调并不等价于逐点物理内部能量单调；二者的差也应随时间步加密而趋于零。

令$$
\Delta_B\mathbf u^{n+1}
=3\mathbf u^{n+1}-4\mathbf u^n+\mathbf u^{n-1},
\qquad
\Delta_Bd^{n+1}
=3d^{n+1}-4d^n+d^{n-1},
$$
$$
\Delta^2U^{n+1}
=U^{n+1}-2U^n+U^{n-1}.
$$
取 $\mathbf z_B^{n+1}$ 为边界量 $3\mathbf u_D^{n+1}-4\mathbf u_D^n+\mathbf u_D^{n-1}$ 的任意有限元提升，并定义一致 Dirichlet 反力功$$
\mathcal R_D^{n+1}(\mathbf z_B^{n+1})
=a_u(\mathbf u^{n+1},\mathbf z_B^{n+1})
+q^{n+1}\overline{\mathcal N}_u^{\,n+1}(\mathbf z_B^{n+1})
-\ell_{\mathrm{ext}}(t_{n+1};\mathbf z_B^{n+1}).
$$
本步的代数 BDF2 外功为$$
W_{\mathrm{ext},B}^{n+1}
=\frac12\left[
\ell_{\mathrm{ext}}(t_{n+1};\Delta_B\mathbf u^{n+1})
+\mathcal R_D^{n+1}(\mathbf z_B^{n+1})
\right].
$$

### 7.1 原始能量直接闭合

若 $q^{n+1}$ 是第 6 节一维原始能量方程 $F_n(q)=0$ 的合格根，则离散场方程和精确标量关系给出$$
\begin{aligned}
\mathscr G_{\mathrm{RLM}}^{n+1}
-\mathscr G_{\mathrm{RLM}}^n
&+\frac14
\|\Delta^2U^{n+1}\|_{\mathcal L}^2\\
&+\frac{\eta}{4\Delta t}
\|\Delta_Bd^{n+1}\|_{L^2(\Omega)}^2\\
&=W_{\mathrm{ext},B}^{n+1}.
\end{aligned}
$$
这里没有 QM 缺口，也没有方向代理能量；$\mathcal E_1(U^{n+1})$、$\mathcal E_1(U^n)$ 和 $\mathcal E_1(U^{n-1})$ 都是原始 Miehe 非线性能量在相应离散状态上的直接求积值。

### 7.2 方向二次上界闭合

若采用第 5 节的 QM 二次方程，则定义最终根处的上界缺口$$
\varepsilon_{\mathrm{QM}}^{n+1}
=\Pi_n(q^{n+1})-\mathcal E_1(U^{n+1})\ge0.
$$
在线性分支、标量方程和外功离散均求解一致时，可得到离散恒等式$$
\begin{aligned}
\mathscr G_{\mathrm{RLM}}^{n+1}
-\mathscr G_{\mathrm{RLM}}^n
&+\frac14
\|\Delta^2U^{n+1}\|_{\mathcal L}^2\\
&+\frac{\eta}{4\Delta t}
\|\Delta_Bd^{n+1}\|_{L^2(\Omega)}^2
+\frac32\varepsilon_{\mathrm{QM}}^{n+1}\\
&=W_{\mathrm{ext},B}^{n+1}.
\end{aligned}
$$
因此只要 QM 上界成立，就有功—能不等式$$
\begin{aligned}
\mathscr G_{\mathrm{RLM}}^{n+1}
-\mathscr G_{\mathrm{RLM}}^n
&\le
W_{\mathrm{ext},B}^{n+1}
-\frac14\|\Delta^2U^{n+1}\|_{\mathcal L}^2\\
&-\frac{\eta}{4\Delta t}
\|\Delta_Bd^{n+1}\|_{L^2(\Omega)}^2.
\end{aligned}
$$
其中黏性耗散也可写成$$
\frac{\eta}{4\Delta t}\|\Delta_Bd^{n+1}\|^2
=\eta\Delta t\|D_2d^{n+1}\|^2.
$$
当 $W_{\mathrm{ext},B}^{n+1}=0$ 时，原始能量直接闭合使 $\mathscr G_{\mathrm{RLM}}^n$ 按精确等式耗散；QM 闭合则额外包含非负缺口 $3\varepsilon_{\mathrm{QM}}^{n+1}/2$。这里的无条件稳定是关于 $\Delta t$ 的能量估计，并不保证任意大时间步都能解析起裂和快速损伤，也不保证每一步自动存在可接受标量根。

## 8. 启动、逐步算法与失败回退

BDF2 需要 $(U^0,q^0)$ 和 $(U^1,q^1)$ 两层数据。取 $q^0=1$，首步可以沿用 [[RLM-BDF1_Mode_I_Miehe]] 的四个 BDF1 线性分支，但把 QM 标量方程换成原始能量直接闭合。定义$$
\phi_0(q)=\mathcal E_1(U_a^1+qU_b^1),
$$
$$
c_0^0
=\mathcal N_u^0(\mathbf u_a^1-\mathbf u^0)
+\mathcal N_d^0(d_a^1-d^0),
\qquad
c_1^0
=\mathcal N_u^0(\mathbf u_b^1)
+\mathcal N_d^0(d_b^1),
$$
则首步标量方程为$$
F_0(q)
=\phi_0(q)+(\alpha-c_1^0)q^2-c_0^0q
-\mathcal E_1(U^0)-\alpha(q^0)^2=0.
$$
它同样可用第 6 节的括区和 Brent 流程求解，因此整段计算都不需要 QM 代理。对足够光滑的解，单个 BDF1 启动步产生 $O(\Delta t^2)$ 的首步误差，通常足以维持后续 BDF2 的全局二阶精度；若问题刚性很强或启动阶段立即起裂，可改用二阶 Crank–Nicolson、SDIRK2 或经验证的两个半步启动。

每个 $n\ge1$ 的时间步按以下顺序执行：

1. 计算 $\overline U^{\,n+1}=2U^n-U^{n-1}$，并组装 $\overline{\mathcal N}_u^{\,n+1}$ 与 $\overline{\mathcal N}_d^{\,n+1}$。
2. 用预分解矩阵求解 $\mathbf u_a^{n+1}$、$\mathbf u_b^{n+1}$、$d_a^{n+1}$ 和 $d_b^{n+1}$。
3. 计算 $c_0^n$、$c_1^n$ 和历史常数，建立只需能量求积的 $F_n(q)$。
4. 首选原始能量路线：扫描并扩展括区，用 Brent 求解全部候选异号区间，验收残差后选择连续物理解支。
5. 若选用 QM 路线：计算 $\phi_n^*$、$g_n^*$ 和 $S_n$，组装并求解 $A_nq^2+B_nq+C_n=0$，再检查方向上界。
6. 接受根后更新 $U^{n+1}=U_a^{n+1}+q^{n+1}U_b^{n+1}$，记录原始能量、外功、耗散和对应路线的 $G$-能量残差。
7. 若搜索区间内没有合格根，减小 $\Delta t$，重新生成外推状态并重做本步。

无论采用哪条路线，都不应接受非零标量残差、搜索区间外的根或与 $q\equiv1$ 不连续的跳根。QM 路线还不能忽略负判别式和上界失败；直接路线则不能用 $|F_n|$ 的非零局部极小点冒充方程根。

## 9. 二阶一致性与建议诊断量

场方程的 BDF2 差分与二阶外推在光滑时间区间内均为二阶。原始能量直接闭合没有 QM 截断误差，但仍要求启动值二阶准确、所选根分支连续、标量求根误差显著小于时间离散误差，并要求解在当前时间区间内足够光滑。

若采用 QM，仍需额外控制方向代理误差。由于标量方程相当于把时间导数乘以 $2\Delta t$ 后离散，一个便于检查的充分尺度是$$
\varepsilon_{\mathrm{QM}}^{n+1}=O(\Delta t^3).
$$
若 $q^{n+1}-1=O(\Delta t^2)$，且 $\Pi_n$ 在 $q=1$ 与 $\phi_n$ 一阶相切，则通常有 $\varepsilon_{\mathrm{QM}}^{n+1}=O(\Delta t^4)$，不会降低二阶精度；但这一性质应由网格和时间步收敛试验确认，不能只由能量稳定性推出。

每步至少记录以下量：

- $|q^{n+1}-1|$ 与 $|q^{n+1}-(2q^n-q^{n-1})|$；
- 原始能量路线的 $F_n(q^{n+1})$ 归一化残差，或 QM 路线的二次多项式归一化残差；
- Brent 最终括区宽度、函数计算次数和本步找到的合格根数；
- QM 路线额外记录 $\varepsilon_{\mathrm{QM}}^{n+1}$ 及其与 $\Delta t^3$ 的比值；
- 位移与相场线性方程残差；
- $G$-能量平衡残差，其中直接路线取 $\chi^{n+1}=0$，QM 路线取 $\chi^{n+1}=3\varepsilon_{\mathrm{QM}}^{n+1}/2$：$$
\begin{aligned}
r_G^{n+1}={}&
\mathscr G_{\mathrm{RLM}}^{n+1}-\mathscr G_{\mathrm{RLM}}^n
-W_{\mathrm{ext},B}^{n+1}\\
&+\frac14\|\Delta^2U^{n+1}\|_{\mathcal L}^2
+\frac{\eta}{4\Delta t}\|\Delta_Bd^{n+1}\|^2
+\chi^{n+1}.
\end{aligned}
$$

在离散、求积和线性求解完全一致时，$r_G^{n+1}$ 应接近标量求根、能量求积和线性求解容差的组合。若只检查 $\mathscr G_{\mathrm{RLM}}^{n+1}\le\mathscr G_{\mathrm{RLM}}^n$，则无法区分真正的稳定耗散与实现误差。

## 10. 三种离散闭合的直接对应

| 项目 | BDF1–QM | BDF2 原始能量直接闭合 | BDF2–QM |
|---|---|---|---|
| 相场时间差分 | $(d^{n+1}-d^n)/\Delta t$ | $(3d^{n+1}-4d^n+d^{n-1})/(2\Delta t)$ | 同左 |
| 非线性力状态 | $U^n$ | $2U^n-U^{n-1}$ | 同左 |
| 相场矩阵质量系数 | $\eta/\Delta t$ | $3\eta/(2\Delta t)$ | 同左 |
| 场分支数 | 四个 | 四个 | 四个 |
| 标量闭合 | 一元二次方程 | 一维非线性方程 $F_n(q)=0$ | 一元二次方程 |
| 标量求解 | 解析稳定二次公式 | Brent 或带保护的 Newton/割线法 | 解析稳定二次公式 |
| 每次标量迭代 | 不需重复能量积分 | 一次 $\mathcal E_1(U_a+qU_b)$ 求积 | 不需重复能量积分 |
| 非线性能量 | 方向二次上界 | 原始 $\mathcal E_1$ | 方向二次上界 |
| 能量关系 | 一步不等式 | 两步精确 $G$-能量等式 | 带 QM 缺口的两步等式 |
| 启动要求 | 只需初值 | 需要二阶准确的第一步 | 需要二阶准确的第一步 |

原始能量直接闭合的主要代价是若干次积分点能量计算，而不是额外的全局线性或非线性场求解。对于大规模有限元问题，只要全局线性求解远贵于一次能量归约，这通常是值得优先尝试的方案。

## 11. 适用范围与下一步验证

- $\eta$ 仍是黏性材料或正则化参数，$\Delta t$ 仍是物理时间分辨率；二者不能互换。
- 固定步长推导不能直接用于变步长 BDF2。变步长时差分系数、外推系数、仿射分支、标量方程和 $G$-稳定矩阵都必须随步长比重新推导。
- 当前格式仍不含不可逆约束、区间约束和惯性项。事后裁剪 $d$ 会破坏变分关系与功—能恒等式。
- Miehe 谱能量在主应变换号或重根处的有限光滑性可能导致外推误差放大、Brent 残差曲线出现折点、QM 的 $S_n$ 回溯频繁或时间阶暂时下降。
- 建议先用光滑制造解验证 $L^2$ 时间二阶，再做单边缺口 Mode-I 的升载—保持试验，并分别检查反力—位移曲线、起裂时间、峰值反力、$q$ 偏差、标量迭代次数、原始能量闭合和 $G$-能量闭合。
- 与 BDF1 对照时，应在同一空间网格、同一求积规则和同一物理加载历史下做 $\Delta t$ 加密，避免把空间误差、谱分解误差或载荷插值误差误判为时间离散误差。

## 12. 参考文献

1. C. Miehe, F. Welschinger, M. Hofacker, “Thermodynamically consistent phase-field models of fracture: Variational principles and multi-field FE implementations,” *International Journal for Numerical Methods in Engineering*, 83 (2010), 1273–1311. [DOI](https://doi.org/10.1002/nme.2861)
2. Q. Cheng, C. Liu, J. Shen, “A new Lagrange multiplier approach for gradient flows,” *Computer Methods in Applied Mechanics and Engineering*, 367 (2020), 113070. [DOI](https://doi.org/10.1016/j.cma.2020.113070)
3. X. Jing, J. Zhao, “Relaxed Lagrange Multiplier (RLM) Schemes for Phase Field Models Preserving the Relaxed Original Energy Dissipation Law,” arXiv:2607.00355 (2026). [arXiv](https://arxiv.org/abs/2607.00355)

> [!note]
> 本文给出的是“准静态位移平衡 + 黏性相场梯度流 + Miehe 谱分解”的固定步长 BDF2 构造，并同时给出原始非线性能量直接闭合与方向二次上界闭合。完整的二阶收敛定理、标量根的全局存在唯一性、不可逆约束下的稳定性以及变步长 $G$-稳定性仍需分别证明。
