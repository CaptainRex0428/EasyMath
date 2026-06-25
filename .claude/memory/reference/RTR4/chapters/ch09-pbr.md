# Ch9 Physically Based Shading 基于物理的着色（**核心参考**）

> **与 EasyMath 关联**：本章节是 PBR 渲染的物理基础。当 `Color.h` 扩展 PBR 相关 API 时，公式参考本章节。
> **完成度**：P0（核心参考）— 完整抽取，2026-06-25

**核心问题**：如何用物理原理建模光与材质的交互？

**章节定位**：实时渲染的"基础物理"——所有现代 PBR 引擎（Unreal、Unity HDRP、Filament、Disney BRDF）的理论基础。

---

## 9.1 光的物理学

### 9.1.1 光波基础

- 光是**电磁横波**：电场 E 和磁场 H 在垂直于传播方向的平面上振荡
- E 和 H 相互垂直，幅度比固定
- 单色光：单一波长 $\lambda$（400-700 nm 可见光）
- 多色光：多个波长的组合
- 偏振：线偏振 vs 非偏振
- 相位速度（真空中光速）$c \approx 3 \times 10^8$ m/s

**能量与振幅**：能量流（irradiance）$\propto$ 振幅²

**干涉**（关键现象）：
- **相长干涉**：$n$ 个相同波叠加，irradiance 是单波的 $n^2$ 倍
- **相消干涉**：$n$ 个相反相位的波，irradiance 为 0
- 不违反能量守恒——空间分布平均后守恒

### 9.1.2 粒子散射

- **瑞利散射**（Rayleigh）：粒子直径 < 波长/10
  - 蓝色天空
  - 散射 $\propto 1/\lambda^4$（蓝光散射强 4 倍）
  - 方向：前向/后向对称

- **米氏散射**（Mie）：粒子直径 ≈ 或 > 波长
  - 云、雾
  - 强前向散射
  - 波长依赖性弱

- **Tyndall 散射**：胶体中的特殊散射

- **聚集团簇效应**：$n$ 个分子聚集成团簇，散射增强 $n$ 倍（相同密度下）

### 9.1.3 介质（Medium）

- **均匀介质**：所有分子相同且均匀分布
  - 不同分子散射波相消干涉 → **无散射**
  - 只改变相位速度（折射率 $n$）和振幅（吸收）
- **复折射率**：$\hat{n} = n + i\kappa$（含衰减率 $\kappa$）
- **非均匀介质**：均匀介质 + 散射粒子

**非均匀介质的视觉表现**：
- 散射 → 灰暗、不透明
- 吸收 → 颜色变化
- 高散射 + 低吸收 → 白色

### 9.1.4 表面

**关键事实**：从光学看，表面 = 不同折射率介质的分界面。

**两侧折射率**：$n_1$（外）、$n_2$（内）

**表面光学性质**：
- 表面要求平行于表面的电场分量连续
- 散射波只向两个方向：反射 + 透射
- 频率不变（无荧光）

**反射角 = 入射角**（Snell 反射定律）
**透射角**：$\eta_1 \sin\theta_1 = \eta_2 \sin\theta_2$（Snell 折射定律）

**金属 vs 电介质**：
- 金属：自由电子 → 高吸收 + 高反射
- 电介质：透明 → 折射进入内部

**几何光学 vs 波动光学**：
- 几何光学：忽略干涉、衍射（光 = 射线）
- 波动光学：1-100 波长范围的细节（章节 9.11）

**微观几何**：
- 表面 < 1 像素 → 微表面
- 多种微表面法线 → 不同方向反射
- 统计建模 → NDF

### 9.1.5 次表面散射

**问题**：光进入表面后多次散射再出射（皮肤、蜡、玛瑙）。

**关键尺度**：
- 进出距离 < 着色尺度 → 局部次表面散射（漫反射）
- 进出距离 > 着色尺度 → 全局次表面散射（章节 14.6）

---

## 9.2 相机模型

### 9.2.1 针孔相机（理想化）

- 没有镜头，极小光圈（数学点）
- 每个传感器收集一束光线（圆锥）
- 适合渲染近似

### 9.2.2 透镜相机

- 大光圈 → 收集更多光
- 引入景深（DoF）
- 镜头带来：模糊、眩光、畸变

**辐射度量 vs 像素着色**：
- 着色样本 → 1 条光线
- AA → 重建连续信号

---

## 9.3 The BRDF（双向反射分布函数）

### 9.3.1 BRDF 定义

$$
f_r(\mathbf{l}, \mathbf{v}) = \frac{dL_o(\mathbf{v})}{dE(\mathbf{l})}
$$

**物理意义**：从方向 $\mathbf{l}$ 入射的辐照度 $E$ → 方向 $\mathbf{v}$ 出射的辐射度 $L$ 的"比例"。

**入射方向/出射方向**：
- 入射光方向 $\mathbf{l}$（指向表面）
- 观察方向 $\mathbf{v}$（指向相机）

### 9.3.2 BRDF 的参数化

**球坐标**：
- $\theta$：相对于法线的仰角
- $\phi$：绕法线的方位角

**各向同性 BRDF**：仅 3 个标量（$\theta_l, \theta_v, \phi_l - \phi_v$）。

**各向异性 BRDF**：4 个标量（$\theta_l, \theta_v, \phi_l, \phi_v$）。

### 9.3.3 反射方程（Reflectance Equation）

$$
L_o(\mathbf{v}) = \int_{\Omega^+} f_r(\mathbf{l}, \mathbf{v}) L_i(\mathbf{l}) (\mathbf{n} \cdot \mathbf{l}) \, d\omega_l
$$

**几何项** $(\mathbf{n} \cdot \mathbf{l})$：Lambert 余弦定律（与法线的点积）。

**球坐标展开**：
$$
L_o(\mathbf{v}) = \int_{0}^{2\pi} \int_{0}^{\pi/2} f_r(\mathbf{l}, \mathbf{v}) L_i(\mathbf{l}) \cos\theta_l \sin\theta_l \, d\theta_l \, d\phi_l
$$

**半向量表示**（重要）：常用 $\mathbf{h} = (\mathbf{l} + \mathbf{v}) / ||\mathbf{l} + \mathbf{v}||$。

### 9.3.4 BRDF 方向的范围

BRDF 只在 $\mathbf{n} \cdot \mathbf{l} > 0$ 且 $\mathbf{n} \cdot \mathbf{v} > 0$ 时有定义。

**插值法线的副作用**：实际可能产生 $\mathbf{n} \cdot \mathbf{l} < 0$（如法线贴图边缘）。

**实用处理**：
- `max(0, dot)`（最常用，简单）
- `abs(dot) + 0.0001`（寒霜引擎做法，避免除以 0）
- **soft clamp**：超过 $\pi/2$ 时逐渐趋于 0

### 9.3.5 BRDF 的两个物理约束

#### Helmholtz 互易性

$$
f_r(\mathbf{l}, \mathbf{v}) = f_r(\mathbf{v}, \mathbf{l})
$$

入射和出射可以互换。

**实际**：实时 BRDF 常违反互易性，无明显瑕疵（但 path tracing 严格需要）。

#### 能量守恒

$$
\forall \mathbf{l}: \quad R(\mathbf{l}) = \int_{\Omega^+} f_r(\mathbf{l}, \mathbf{v}) (\mathbf{n} \cdot \mathbf{v}) \, d\omega_v \le 1
$$

**定向半球反射率** $R(\mathbf{l})$ 必须在 $[0, 1]$ 内。

**注意**：BRDF 本身可以 $> 1$（高光中心是高度集中的分布），但积分后必须 $\le 1$。

### 9.3.6 各向同性参数化示例

球坐标 $(\theta, \phi)$：

$$
L_o(\mathbf{v}) = \int_0^{2\pi} \int_0^{\pi/2} f_r(\theta_l, \phi_l, \theta_v, \phi_v) L_i(\theta_l, \phi_l) \cos\theta_l \sin\theta_l \, d\theta_l \, d\phi_l
$$

### 9.3.7 Lambertian BRDF

最简单的 BRDF：恒定值。

$$
f_r(\mathbf{l}, \mathbf{v}) = \frac{\rho_{ss}}{\pi}
$$

其中 $\rho_{ss}$ 是**次表面反照率**（subsurface albedo）。

**反射率**：
$$
R_{\text{lambert}} = \rho_{ss}
$$

**术语**：
- $\rho_{ss}$ / albedo / diffuse color
- 注意：albedo 是天文术语（行星反射能力），不是反射率

### 9.3.8 BRDF 的可视化

**输入方向固定，绘制所有输出方向上的 BRDF 值**：
- 漫反射 = 半球
- 镜面 = 椭球（波瓣）
- 镜面波瓣厚度 → 反射模糊程度

**互易性** = 反向可视化（不同入射对单个出射的贡献）。

---

## 9.4 光照

### 9.4.1 方向光（Directional Light）

**极限情况**：面光源缩为无穷小，光源颜色 $L_l$ 固定。

反射方程简化为：
$$
L_o(\mathbf{v}) = f_r(\mathbf{l}, \mathbf{v}) L_l (\mathbf{n} \cdot \mathbf{l})^+ = f_r(\mathbf{l}, \mathbf{v}) L_l \cos\theta_l
$$

### 9.4.2 精确光源

**与方向光类似**，但加入距离衰减：
$$
L_i \propto \frac{1}{r^2}
$$

**注意**：真实物理中��真正的"点光源"，但作为面光源的极限近似是合理的。

### 9.4.3 多光源求和

$$
L_o = \sum_{i=1}^{n} f_r(\mathbf{l}_i, \mathbf{v}) L_i (\mathbf{n} \cdot \mathbf{l}_i)^+
$$

**BRDF 中除以 $\pi$ 的常被反射方程中的 $\pi$ 抵消**——简化代码。

**注意**：用论文 BRDF 时注意是否带 $\pi$ 系数。

---

## 9.5 菲涅尔反射（**重要**）

### 9.5.1 完整菲涅尔方程

菲涅尔方程描述了**入射光在表面上的反射/折射比例**。

**两种偏振**（s 偏振垂直、p 偏振平行）：

$$
F_s = \left| \frac{n_1 \cos\theta_i - n_2 \cos\theta_t}{n_1 \cos\theta_i + n_2 \cos\theta_t} \right|^2
$$

$$
F_p = \left| \frac{n_2 \cos\theta_i - n_1 \cos\theta_t}{n_2 \cos\theta_i + n_1 \cos\theta_t} \right|^2
$$

**非偏振光的反射率**：
$$
F = \frac{F_s + F_p}{2}
$$

**使用半向量**（BRDF 中）：$F(\mathbf{h}, \mathbf{l})$，用 $\mathbf{h}$ 代替法线。

### 9.5.2 Schlick 近似（**实时标准**）

$$
F(\mathbf{h}, \mathbf{l}) \approx F_0 + (1 - F_0) (1 - \mathbf{h} \cdot \mathbf{l})^5
$$

**关键性质**：
- 当 $\mathbf{h} \cdot \mathbf{l} = 1$（法向入射）→ $F = F_0$
- 当 $\mathbf{h} \cdot \mathbf{l} = 0$（掠射）→ $F = 1$（白色）
- 简单、廉价、足够准确

**$F_0$ 计算**（电介质，$n_1 = 1$）：
$$
F_0 = \left( \frac{n_2 - 1}{n_2 + 1} \right)^2
$$

**默认 $F_0$**（非金属）：$0.04$（很多电介质接近此值）

### 9.5.3 修正 Schlick（更精确）

$$
F(\mathbf{h}, \mathbf{l}) \approx F_0 + (1 - F_0) (1 - \mathbf{h} \cdot \mathbf{l})^p
$$

**变化 $p$**：可调过渡的"速度"。Gulbrandsen 给出更精确的金属近似。

### 9.5.4 典型 $F_0$ 值

**电介质**（$F_0 \approx 0.04$ 或更低，无色）：
- 玻璃：0.04
- 水：0.02
- 塑料：0.04-0.05
- 木材：~0.04
- 皮肤：~0.028

**金属**（$F_0 > 0.5$，有色）：

| 金属 | F0 (R, G, B) 线性 |
|------|---------------------|
| 金 (Gold) | (1.000, 0.766, 0.336) |
| 银 (Silver) | (0.972, 0.960, 0.915) |
| 铝 (Aluminum) | (0.913, 0.926, 0.924) |
| 铜 (Copper) | (0.955, 0.638, 0.538) |
| 铁 (Iron) | (0.560, 0.570, 0.580) |
| 镍 (Nickel) | (0.660, 0.609, 0.526) |

### 9.5.5 半导体

性质介于电介质和金属之间（如硅、锗）。
- $F_0$ 在 0.1-0.5 之间
- 波长依赖性强

### 9.5.6 内反射

$n_1 > n_2$ 时，光从内部射出到外部。

**临界角**（全内反射开始）：
$$
\theta_c = \arcsin(n_2 / n_1)
$$

- 临界角以内 → 正常折射
- 临界角以外 → 全内反射（$F = 1$）

**内反射曲线**是外反射曲线的"压缩"版本（临界角处变白）。

**简化计算**：用透射角 $\theta_t$ 代替入射角 $\theta_i$。

**应用**：水下气泡看起来像银球（高反射率）。

**金属不发生内反射**（会迅速吸收内部光线）。

---

## 9.6 微观几何（Microgeometry）

### 9.6.1 关键概念

- 不可见的不规则性（< 1 像素）→ 微观几何
- 微观几何对每个像素可见多个微表面
- 不同微表面有不同法线 → 不同方向反射
- 宏观法线 = 微表面法线的平均

**核心影响**：
- 表面越粗糙 → 反射越模糊
- 反射锥变宽变暗
- 高光变"宽"

### 9.6.2 各向同性 vs 各向异性

- **各向同性**：旋转对称（统一粗糙度）
- **各向异性**：有方向偏好（拉丝金属、头发、木纹）

### 9.6.3 shadowing / masking

- **shadowing**：微表面遮挡入射光
- **masking**：微表面遮挡反射光
- 与微表面高度分布相关

### 9.6.4 高度相关性

- 微表面高度和法线相关 → shadowing 和 masking 有效
- 凹陷部分被遮挡 → 表面看起来"光滑"

### 9.6.5 逆反射（Retroreflection）

- 微观不规则 > 次表面散射距离时
- 粗糙表面优先反射回入射方向
- 例子：路面、积雪

### 9.6.6 相互反射

- 微表面之间多次反射
- 电介质：每次被菲涅尔衰减 → 微弱
- 金属：唯一可见漫反射来源
- 有色金属：多次反射颜色更深

### 9.6.7 几何与着色尺度的关系

- 表面不规则 < 波长 → 忽略（不影响光）
- 表面不规则 > 100× 波长 → 等同平面
- 表面不规则 1-100× 波长 → 衍射（章节 9.11）

---

## 9.7 微表面理论（Microfacet Theory）

### 9.7.1 基本思想

- 微观几何 = 微表面集合
- 每个微表面平坦，有微表面法线 $\mathbf{m}$
- 微表面作为**完美菲涅尔镜面**（micro-BRDF = 菲涅尔镜面）
- 总 BRDF = 所有微表面贡献的和

### 9.7.2 BRDF 通用公式

$$
f_r(\mathbf{l}, \mathbf{v}) = \int_{\Omega_{\mathbf{m}}^+} f_{\mu}(\mathbf{l}, \mathbf{v}, \mathbf{m}) \frac{D(\mathbf{m}) G(\mathbf{l}, \mathbf{v}, \mathbf{m})}{4 (\mathbf{n} \cdot \mathbf{l})(\mathbf{n} \cdot \mathbf{v})} (\mathbf{n} \cdot \mathbf{m}) \, d\omega_{\mathbf{m}}
$$

其中：
- $f_{\mu}$：微表面 BRDF（菲涅尔）
- $D(\mathbf{m})$：NDF（法线分布）
- $G(\mathbf{l}, \mathbf{v}, \mathbf{m})$：masking-shadowing
- 分母 $4 (\mathbf{n} \cdot \mathbf{l})(\mathbf{n} \cdot \mathbf{v})$：雅可比变换

### 9.7.3 NDF 的归一化

**面积归一化**（整个球面积分）：
$$
\int_{\Omega_{\mathbf{m}}} D(\mathbf{m}) \, d\omega_{\mathbf{m}} = 1
$$

**投影面积归一化**（投影到宏观表面）：
$$
\int_{\Omega_{\mathbf{m}}^+} D(\mathbf{m}) (\mathbf{n} \cdot \mathbf{m}) \, d\omega_{\mathbf{m}} = 1
$$

后者更常用（投影面积之和 = 宏观面积）。

### 9.7.4 Masking 函数

**Masking 函数** $G_1(\mathbf{v}, \mathbf{m})$：从方向 $\mathbf{v}$ 看，法线为 $\mathbf{m}$ 的微表面可见的比例。

**联合 masking-shadowing**：$G(\mathbf{l}, \mathbf{v}, \mathbf{m}) = G_1(\mathbf{v}, \mathbf{m}) G_1(\mathbf{l}, \mathbf{m})$（可分离形式，简化但有误差）。

**Smith 函数**：Heitz 2014 证明 Smith masking 是最准确的物理合法函数。

**Smith Lambda 函数**：
$$
G_1(\mathbf{v}, \mathbf{m}) = \frac{\chi^+(\mathbf{v} \cdot \mathbf{m})}{1 + \Lambda(\mathbf{v})}
$$
其中 $\chi^+$ 是指示函数（$\mathbf{v} \cdot \mathbf{m} > 0$ 时为 1）。

**Heitz 2014**：Smith 是唯一满足法线-masking 独立性的合法函数。

---

## 9.8 表面反射 BRDF 模型

### 9.8.1 法线分布函数 NDF

#### Beckmann NDF（经典，光学界用得多）

$$
D(\mathbf{m}) = \frac{1}{\pi \alpha^2 \cos^4\theta_m} \exp\left(-\frac{\tan^2\theta_m}{\alpha^2}\right)
$$

- $\theta_m$：宏观法线 $\mathbf{n}$ 与微表面法线 $\mathbf{m}$ 的夹角
- $\alpha$：与微表面均方根斜率成正比
- 形状不变（shape-invariant）→ 简化 $\Lambda$ 推导

**Smith $\Lambda$（Beckmann）**：
$$
\Lambda(\mathbf{v}) = \frac{\text{erf}(a) - 1}{2} + \frac{1}{2a\sqrt{\pi}} e^{-a^2}
$$
其中 $a = 1/(\alpha \tan\theta_v)$，$\text{erf}$ 是误差函数。

**近似**（避免 erf）：
$$
\Lambda(\mathbf{v}) \approx \begin{cases} \frac{1 - 1.259 a + 0.396 a^2}{3.535 a + 2.181 a^2}, & a < 1.6 \\ 0, & a \ge 1.6 \end{cases}
$$

#### Blinn-Phong NDF（移动端友好）

$$
D(\mathbf{m}) = \frac{n+1}{2\pi} (\mathbf{n} \cdot \mathbf{m})^n
$$

- $n$：粗糙度参数（越大越光滑）
- 简单、廉价
- $\alpha$ 与 $n$ 关系（等价变换）：$\alpha^2 = 2/(n+2)$（Walter 2007）

**非形状不变** → $\Lambda$ 无解析解，用 Beckmann 等价。

#### GGX / Trowbridge-Reitz NDF（**行业标准**）

$$
D(\mathbf{m}) = \frac{\alpha^2}{\pi [(\mathbf{n} \cdot \mathbf{m})^2 (\alpha^2 - 1) + 1]^2}
$$

- **Disney 映射**：$\alpha = \text{roughness}^2$（让滑块线性对应视觉）
- 形状不变
- Smith $\Lambda$ 简单：

$$
\Lambda(\mathbf{v}) = \frac{-1 + \sqrt{1 + \frac{1}{\alpha^2 \tan^2\theta_v}}}{2}
$$

**变量只以 $1/\tan^2\theta_v$ 形式出现** → 避免开方。

**视觉特征**：长尾 → 高光中心锐利 + 周围微光（haze）。

#### 广义 Trowbridge-Reitz (GTR) NDF

$$
D(\mathbf{m}) = \frac{c(\gamma)}{\pi [(\mathbf{n} \cdot \mathbf{m})^2 (\alpha^2 - 1) + 1]^\gamma}
$$

- $\gamma$：控制尾部长度
- $\gamma = 2$：GGX
- $\gamma$ 小：长尾（更接近真实材质）
- $\gamma$ 大：短尾（接近 Beckmann）
- 归一化系数 $c(\gamma)$ 复杂

**非形状不变**。

#### Karis 近似（GGX + Smith）

Karis 2013 给出 Smith $\Lambda$ 的近似：

$$
\Lambda(\mathbf{v}) \approx \frac{0.5 \cdot \text{sign}(\cos\theta_v)}{|\cos\theta_v| + \alpha \tan\theta_v}
$$

**优化**：Hammon 2017 进一步用此简化 GGX + Smith masking 的组合（消掉一些项）。

### 9.8.2 各向异性 NDF

**核心思想**：沿两个主轴有不同粗糙度 $\alpha_x, \alpha_y$。

**Ward 各向异性**（简单）：
$$
D(\mathbf{m}) = \frac{1}{\alpha_x \alpha_y \sqrt{(\mathbf{n} \cdot \mathbf{m})^4}} \exp\left( -\frac{(\mathbf{m} \cdot \mathbf{x}/\alpha_x)^2 + (\mathbf{m} \cdot \mathbf{y}/\alpha_y)^2}{(\mathbf{n} \cdot \mathbf{m})^2} \right)
$$

**GGX 各向异性**：
$$
D(\mathbf{m}) = \frac{1}{\pi \alpha_x \alpha_y \left[ \frac{(m_x/\alpha_x)^2 + (m_y/\alpha_y)^2 + m_z^2}{} \right]^2}
$$

**应用**：拉丝金属、头发、毛发。

### 9.8.3 多次反弹（Kulla-Conty）

**问题**：单次反射 BRDF 会过暗（忽略微表面相互反射）。

**Kulla-Conty 2017**：用额外的 BRDF 项补偿。

$$
L_o \approx L_{\text{direct}} + L_{\text{energyCompensation}}
$$

---

## 9.9 次表面散射 BRDF

### 9.9.1 次表面反照率

**关键参数**：
- $\sigma_s$：散射系数
- $\sigma_a$：吸收系数
- 平均自由程 $1/(\sigma_a + \sigma_s)$

**次表面反照率**：
$$
\rho_{ss} = \frac{\sigma_s}{\sigma_a + \sigma_s}
$$

**简化模型**：漫反射 $f_r = \rho_{ss} / \pi$。

### 9.9.2 散射距离与着色尺度

**关键判断**：
- 散射距离 < 着色尺度 → 局部模型（漫反射 + Lambertian）
- 散射距离 > 着色尺度 → 全局模型（章节 14.6）

### 9.9.3 光滑表面次表面模型

**Wrap lighting**（半 Lambert）：
$$
\text{wrapped} = \max(0, \mathbf{n} \cdot \mathbf{l} + 0.5) / 1.5
$$

软化光照过渡。

### 9.9.4 粗糙表面次表面模型

**法线模糊**：在像素着色中，对法线做小范围模糊。
- 改善粗糙皮肤的次表面效果

**预积分皮肤**：
- 烘焙"blurred normal → translucency"查找表
- 运行时只查表

**Christophe Schlick 经验模型**：
$$
f_r = \frac{\rho}{\pi} \cdot \frac{1}{1 + F(1-\rho) \cdot (1-\mathbf{n} \cdot \mathbf{l})^2}
$$

### 9.9.5 Hanrahan-Krueger 模型

基于辐射传输方程的精确模型：
$$
T(r) = \frac{\alpha'}{4\pi} \left[ (\sigma_{tr} d + 1) \frac{e^{-\sigma_{tr} d}}{\sigma_{tr} d^2} + \frac{e^{-\sigma_{tr} d}}{2d} \right]
$$
其中 $d$ 是厚度，$\sigma_{tr}$ 是有效消光系数。

---

## 9.10 布料 BRDF

布料是**非传统材质**——纤维结构，不是连续表面。

### 9.10.1 经验布料模型

**核心**：增加 sheen 项（模拟绒毛）。

**Ashikhmin 模型**：在主 BRDF 上加 "sheen" lobe。

### 9.10.2 微表面布料模型

**用微面反射**（纤维截面作为微表面）。

### 9.10.3 微圆柱体布料模型

**用细长圆柱体表示纤维**。

**关键参数**：纤维密度、纤维朝向、纤维反射率。

### 9.10.4 布料 NDF

**亮点**：圆柱体 NDF，不是微面 NDF。
- 例：所有朝向随机 → 光线在多个圆柱上散射

---

## 9.11 波动光学 BRDF

**适用**：薄膜干涉（油膜、肥皂泡）、衍射（CD 光栅）。

### 9.11.1 衍射模型

**原理**：微观周期性结构 → 衍射光栅。
- 反射率随波长和方向变化（虹彩效果）
- 例：CD、蝴蝶翅膀

### 9.11.2 薄膜干涉模型

**原理**：上下表面反射光干涉。
- 厚度变化 → 颜色变化

**条件**（相长干涉）：
$$
2 n d \cos\theta_t = m \lambda
$$
其中 $d$ 是膜厚，$n$ 是折射率，$m$ 是整数阶数。

**应用**：油膜、肥皂泡、珍珠。

---

## 9.12 分层材质

**适用**：漆面（透明漆 + 不透明金属）、皮肤（表皮 + 真皮 + 血液）。

**模型**：每层用 BRDF/BSSRDF 表达，层间用 Beer-Lambert 衰减。

**Weidlich-Wilkie 模型**：层间积分。

**简化**：预计算为查找表。

---

## 9.13 混合和过滤材质

### 9.13.1 混合问题

当像素跨两种材质时（边界），直接混合可能产生混叠。

### 9.13.2 过滤法线与法线分布

- **过滤法线**：将粗糙度高的微面求平均
- **过滤 NDF**：NDF 自身具有可过滤性
- **Heitz 2018**：详细推导

### 9.13.3 去相关

**关键问题**：不同材质的法线要"去相关"——避免混合时出现错误高光。

---

## 9.14 BRDF 速查表

| BRDF 模型 | 复杂度 | 真实感 | 适用 |
|-----------|--------|--------|------|
| Lambert | 低 | 低 | 简单着色 |
| Phong | 低 | 低 | 旧游戏 |
| Blinn-Phong | 低 | 中 | 移动端 |
| Cook-Torrance (Beckmann) | 中 | 高 | 经典 PBR |
| Cook-Torrance (GGX) | 中 | 高 | 现代 PBR 标准 |
| GTR | 中 | 极高 | 工业级 |
| Sheen + 标准 | 高 | 极高 | 布料 |

---

## 与 EasyMath 项目关联

### 当前实现
- `Color<T>` 已支持颜色空间，可作为 BRDF 的输入/输出
- 暂无 BRDF 工具函数

### 未来 PBR API 设计的参考

```cpp
// 概念性 API（EasyMath 未来可加）
namespace EasyMath {
    // BRDF 输入
    struct BRDFInputs {
        Vector3 albedo;
        float roughness;
        float metalness;
        Color F0;  // 基色反射率
    };

    // 微表面 BRDF 三件套
    template<typename T>
    T ggxNDF(const Vector<T, 3>& n, const Vector<T, 3>& h, T alpha);

    template<typename T>
    T smithLambda(T cosTheta, T alpha);

    template<typename T>
    T smithG1(T cosThetaV, T cosThetaL, T alpha);

    template<typename T>
    Color<T> schlickFresnel(const Color<T>& F0, T hDotL);

    // 完整 Cook-Torrance
    template<typename T>
    Color<T> evaluateBRDF(const BRDFInputs& mat,
                          const Vector<T, 3>& n,
                          const Vector<T, 3>& v,
                          const Vector<T, 3>& l);

    // 几何项
    template<typename T>
    T geometryTermGGX(T nDotV, T nDotL, T alpha);

    // 各向异性版本（未来）
    template<typename T>
    T ggxNDFAniso(const Vector<T, 3>& n, const Vector<T, 3>& h,
                   T alphaX, T alphaY);
}
```

### 公式参考（已实现公式的对应）

- **Rec.709 亮度系数**（Ch8）：EasyMath `Color::luminance()` 用 $0.2126, 0.7152, 0.0722$
- **sRGB ↔ Linear 公式**（Ch8）：`Color::toLinear()` / `Color::toSRGB()` 用分段函数
- **WCAG 对比度**（Ch8）：`Color::contrastRatio()` 用 $(L_1 + 0.05) / (L_2 + 0.05)$

## 关键术语索引

| 术语 | 释义 |
|------|------|
| BRDF | 双向反射分布函数 |
| 微表面 | 微观尺度的随机倾斜小面 |
| NDF | 法线分布函数 |
| 菲涅尔 | 入射角相关的反射率 |
| Schlick 近似 | 菲涅尔的快速近似 |
| GGX | 流行的 NDF 形式 |
| Trowbridge-Reitz | GGX 的原名 |
| GTR | 广义 Trowbridge-Reitz |
| 几何遮蔽 | masking + shadowing |
| 互易性 | BRDF 的对称性 |
| 能量守恒 | BRDF 半球积分 $\le 1$ |
| $F_0$ | 法向入射反射率 |
| 多次反弹 | 微面间的互相反射 |
| Kulla-Conty | 能量补偿近似 |
| 微圆柱体 | 布料的纤维表示 |
| sheen | 布料绒毛的额外 lobe |
| 薄膜干涉 | 上下表面光干涉 |
| 衍射 | 周期性结构的虹彩效果 |
| $V(\mathbf{v})$ | Smith masking 函数 |
| $\Lambda$ | Smith Lambda 函数 |
| 形状不变 | NDF 缩放可统一（GGX 是） |
| $\alpha$ | 粗糙度参数 |

## 关键公式

### 微表面 BRDF（**Cook-Torrance**）

$$
f_r(\mathbf{l}, \mathbf{v}) = \frac{F(\mathbf{h}, \mathbf{l}) \cdot G(\mathbf{l}, \mathbf{v}) \cdot D(\mathbf{h})}{4 (\mathbf{n} \cdot \mathbf{l})(\mathbf{n} \cdot \mathbf{v})}
$$

### Schlick 菲涅尔

$$
F(\mathbf{h}, \mathbf{l}) = F_0 + (1 - F_0) (1 - \mathbf{h} \cdot \mathbf{l})^5
$$

### GGX NDF

$$
D(\mathbf{m}) = \frac{\alpha^2}{\pi [(\mathbf{n} \cdot \mathbf{m})^2 (\alpha^2 - 1) + 1]^2}
$$

### Smith $\Lambda$（GGX）

$$
\Lambda(\mathbf{v}) = \frac{-1 + \sqrt{1 + \frac{1}{\alpha^2 \tan^2\theta_v}}}{2}
$$

### Smith Masking

$$
G_1(\mathbf{v}, \mathbf{m}) = \frac{1}{1 + \Lambda(\mathbf{v})}
$$

### Karis 近似

$$
\Lambda(\mathbf{v}) \approx \frac{0.5}{|\cos\theta_v| + \alpha \tan\theta_v}
$$

## 抽取历史

- 2026-06-25：基于 `RTR4/_raw/ch009.txt`（168KB）抽取
- 来源章节号：Ch9
- 跳过内容：Ch9.6-9.7 部分示例代码（保留公式即可）
- 置信度：高
