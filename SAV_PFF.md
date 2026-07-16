# 准静态相场断裂模型的标量辅助变量（SAV）完全线性化数值方案知识库

本知识库提炼了基于标量辅助变量（SAV, Scalar Auxiliary Variable）方法构建的准静态相场断裂模型完全线性化、无迭代交替求解（Staggered）数值算法的核心理论与实现步骤。

---

## 1. 物理模型与连续体SAV公式化

### 1.1 连续体控制方程
在引入粘性耗散势 $\Phi_{\text{visc}}(\dot{d}) = \frac{\eta}{2} \dot{d}^2$ 后，动态相场断裂系统的控制方程为：
*   **位移的静态平衡方程（Momentum Balance）**：
    $$\nabla \cdot \sigma + \mathbf{b} = 0 \quad \text{in } \Omega$$
*   **裂纹的相场演化方程（Phase-Field Evolution）**：
    $$\eta\dot{d} = -\frac{\partial E}{\partial d} = -g'(d)\mathcal{H} - G_c \left( \frac{1}{\ell_c} d - \ell_c\nabla^2 d \right) \quad \text{in } \Omega$$

其中，$E$ 为总能量泛函，$\mathbf{u}$ 为位移矢量，$\mathbf{b}$ 为外力。$d \in [0, 1]$ 为标量损伤相场，$d=1$ 代表完全损伤，$d=0$ 代表完好，$\eta \ge 0$ 为动力学粘性系数，$G_c$ 为临界能量释放率，$\ell_c$ 为特征长度尺度。$\sigma = \frac{\partial \psi}{\partial \varepsilon}$ 为柯西应力张量，$\mathcal{H} = \max\{\psi^+\}$ 为历史应变能场，用于满足裂纹不可逆约束。

### 1.2 引入标量辅助变量 (SAV)
为了处理体弹性应变能 $E_{\text{bulk}}(\mathbf{u}, d)$ 与相场 $d$ 之间的强非线性耦合，引入一个全局的时间相关标量辅助变量 $r(t)$：
$$r(t) = \sqrt{E_{\text{bulk}}(\mathbf{u}, d) + S_0} = \sqrt{\int_{\Omega} \left[ g(d)\psi_0^+(\varepsilon) + \psi_0^-(\varepsilon) \right] \text{d}\Omega + S_0}$$

其中 $S_0 > 0$ 为正稳定化常数，用于保证根号内的项严格大于零。通过链式法则，其连续时间导数 $\dot{r}(t)$ 展开为：
$$\dot{r}(t) = \frac{1}{2\sqrt{E_{\text{bulk}} + S_0}} \int_{\Omega} \left( \frac{\partial \psi(\varepsilon, d)}{\partial \varepsilon} : \dot{\varepsilon} + \frac{\partial \psi(\varepsilon, d)}{\partial d} \dot{d} \right) \text{d}\Omega$$

---

## 2. 完全线性化数值方案与能量分裂

由于受损的柯西应力张量高度非线性且非凸，在隐式求解时需要高成本的非线性迭代。为此，引入**线性稳定化技术**，将总体变弹性应变能 $E_{\text{bulk}}(\mathbf{u}, d)$ 分裂为常数线性参考弹性应能 $E_0(\mathbf{u})$（作为隐式稳定器）和非线性/受损偏差能 $E_{\text{nl}}(\mathbf{u}, d)$（进行显式处理）：

$$E_{\text{bulk}}(\mathbf{u}, d) = E_0(\mathbf{u}) + E_{\text{nl}}(\mathbf{u}, d)$$

其中：
*   **参考无损能量**：
    $$E_0(\mathbf{u}) = \frac{1}{2} \int_{\Omega} \varepsilon(\mathbf{u}) : \mathbb{C}_0 : \varepsilon(\mathbf{u}) \text{d}\Omega$$
*   **非线性偏差能量**：
    $$E_{\text{nl}}(\mathbf{u}, d) = \int_{\Omega} \left[ \psi(\varepsilon(\mathbf{u}), d) - \frac{1}{2} \varepsilon(\mathbf{u}) : \mathbb{C}_0 : \varepsilon(\mathbf{u}) \right] \text{d}\Omega$$

这里，$\mathbb{C}_0 = \mathbb{C}_{\text{pristine}}$ 为材料无损状态下的常数第四阶弹性张量（对于各向同性材料，$\mathbb{C}_0 : \varepsilon = \lambda \text{tr}(\varepsilon)\mathbf{I} + 2\mu\varepsilon$）。

此时，重新定义标量辅助变量 $r(t)$，使其仅缩放非线性偏差部分：
$$r(t) = \sqrt{E_{\text{nl}}(\mathbf{u}, d) + S_0}$$

---

## 3. 时间离散多步更新系统

设当前时间步 $t^n$ 的状态变量 $(\mathbf{u}^n, d^n, r^n)$ 及历史应变能场 $\mathcal{H}^n$ 已知，更新至 $t^{n+1}$ 的完全线性化 staggered 格式定义为以下多行耦合系统：

$$\mathcal{H}^n = \max \left( \mathcal{H}^{n-1}, \psi_0^+(\varepsilon(\mathbf{u}^n)) \right), \quad \xi^n = \frac{\mathcal{H}^n}{\sqrt{E_{\text{nl}}^n + S_0}}$$

$$\eta \frac{d^{n+1} - d^n}{\Delta t} - G_c \ell_c \nabla^2 d^{n+1} + \frac{G_c}{\ell_c} d^{n+1} = 2 r^{n+1}(1 - d^n)\xi^n$$

$$- \nabla \cdot \left( \mathbb{C}_0 : \varepsilon(\mathbf{u}^{n+1}) \right) = \nabla \cdot \left[ \frac{r^{n+1}}{\sqrt{E_{\text{nl}}^n + S_0}} \left( \sigma^n - \mathbb{C}_0 : \varepsilon(\mathbf{u}^n) \right) \right] + \mathbf{b}^{n+1}$$

$$r^{n+1} = r^n + \frac{1}{2\sqrt{E_{\text{nl}}^n + S_0}} \int_{\Omega} \left[ \tilde{\sigma}^n : \left( \varepsilon(\mathbf{u}^{n+1}) - \varepsilon(\mathbf{u}^n) \right) + \delta_d\psi^n \left( d^{n+1} - d^n \right) \right] \text{d}\Omega$$

其中，在 $t^n$ 步计算的各本构项为：
*   $\tilde{\sigma}^n = \sigma^n - \mathbb{C}_0 : \varepsilon(\mathbf{u}^n)$
*   $\sigma^n = g(d^n) \frac{\partial \psi_0^+}{\partial \varepsilon}(\varepsilon(\mathbf{u}^n)) + \frac{\partial \psi_0^-}{\partial \varepsilon}(\varepsilon(\mathbf{u}^n))$
*   $\delta_d\psi^n = g'(d^n)\psi_0^+(\varepsilon(\mathbf{u}^n)) = -2(1-d^n)\psi_0^+(\varepsilon(\mathbf{u}^n))$

---

## 4. 基于线性叠加（Linear Superposition）的非迭代求解算法

尽管系统中的 $r^{n+1}$、$\mathbf{u}^{n+1}$ 和 $d^{n+1}$ 在全局上相互耦合，但该耦合是**严格线性**的。通过类 Sherman-Morrison 分割技术，可将更新场表示为基准空间解的线性组合，从而完全规避非线性迭代。

### 4.1 步骤 1：求解相场基准解
由于相场方程关于未知数 $r^{n+1}$ 呈线性关系，设 $d^{n+1} = d_1 + r^{n+1} d_2$。代入离散相场方程，可解耦为两个线性椭圆型偏微分方程：
$$\left( \frac{\eta}{\Delta t} + \frac{G_c}{\ell_c} - G_c \ell_c \nabla^2 \right) d_1 = \frac{\eta}{\Delta t} d^n$$
$$\left( \frac{\eta}{\Delta t} + \frac{G_c}{\ell_c} - G_c \ell_c \nabla^2 \right) d_2 = 2(1 - d^n)\xi^n$$

> **计算优势**：两方程共享同一个对称正定（SPD）的左端线性算子 $\mathbf{K}_d$。该算子在时间步循环外仅需在 $t=0$ 时进行一次 Cholesky 分解。在每个时间步，求解 $d_1$ 和 $d_2$ 仅需进行两次快速的回代消元（Back-substitution）。

### 4.2 步骤 2：求解位移基准解
同理，设 $\mathbf{u}^{n+1} = \mathbf{u}_1 + r^{n+1} \mathbf{u}_2$。代入动量平衡方程，解耦为两个线性弹性力学方程：
$$- \nabla \cdot \left( \mathbb{C}_0 : \varepsilon(\mathbf{u}_1) \right) = \mathbf{b}^{n+1}$$
$$- \nabla \cdot \left( \mathbb{C}_0 : \varepsilon(\mathbf{u}_2) \right) = \frac{1}{\sqrt{E_{\text{nl}}^n + S_0}} \nabla \cdot \left[ \sigma^n - \mathbb{C}_0 : \varepsilon(\mathbf{u}^n) \right]$$

> **计算优势**：两偏微分方程同样共享恒定的线弹性刚度算子 $\mathbf{K}_u = \mathbf{K}_0$。该算子在整个模拟周期内为常数，只需在初始时刻进行一次因式分解，每步求解只需两次回代消元。

### 4.3 步骤 3：求解全局标量 $r^{n+1}$ 
将基准解关系式 $d^{n+1} = d_1 + r^{n+1} d_2$ 和 $\mathbf{u}^{n+1} = \mathbf{u}_1 + r^{n+1} \mathbf{u}_2$ 代入 $r^{n+1}$ 的更新方程，得到关于全局标量 $r^{n+1}$ 的一元线性代数方程：
$$r^{n+1} = \frac{B}{A}$$

其中，标量常数 $A$ 和 $B$ 通过全域数值积分计算：
$$A = 1 - \frac{1}{2\sqrt{E_{\text{nl}}^n + S_0}} \int_{\Omega} \left[ \tilde{\sigma}^n : \varepsilon(\mathbf{u}_2) + \delta_d\psi^n d_2 \right] \text{d}\Omega$$
$$B = r^n + \frac{1}{2\sqrt{E_{\text{nl}}^n + S_0}} \int_{\Omega} \left[ \tilde{\sigma}^n : \left( \varepsilon(\mathbf{u}_1) - \varepsilon(\mathbf{u}^n) \right) + \delta_d\psi^n \left( d_1 - d^n \right) \right] \text{d}\Omega$$

> **计算优势**：在离散层面上，这些积分直接退化为简单的向量-向量内积（Vector-vector inner products），计算复杂度为极低的 $\mathcal{O}(N)$。

### 4.4 步骤 4：物理场与速度重建
在显式确定 $r^{n+1}$ 的值后，通过线性叠加重建真实的 $t^{n+1}$ 步相场和位移场，并显式更新速度场：
*   $$d^{n+1} = d_1 + r^{n+1} d_2$$
*   $$\mathbf{u}^{n+1} = \mathbf{u}_1 + r^{n+1} \mathbf{u}_2$$

---

## 5. 算法综合流程图 (逻辑参考)

```
[已知第 n 步状态变量] -> 更新历史应变能场 H^n, 缩放因子 ξ^n
                             |
                             v
           [Cholesky 回代消元 (利用初始化时分解的 Kd, Ku)]
           ├── 求解相场基准解 d1, d2
           └── 求解位移基准解 u1, u2
                             |
                             v
           [组装积分 A, B (向量内积 O(N) 复杂度)]
           └── 计算全局标量 r^(n+1) = B / A
                             |
                             v
           [物理场重建与更新]
           ├── d^(n+1) = d1 + r^(n+1) * d2
           └── u^(n+1) = u1 + r^(n+1) * u2
```

该方案使得每一个时间步的更新仅涉及**四次回代消元**和**两次向量内积**，在保证完全线性化和不需任何非线性迭代的前提下，实现了极高的计算效率。