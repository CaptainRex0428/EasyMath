# Ch9 Physically Based Shading 基于物理的着色（**核心参考**）

> **与 EasyMath 关联**：本章节是 PBR 渲染的物理基础。当 `Color.h` 扩展 PBR 相关 API 时，公式参考本章节。

**核心问题**：如何用物理原理建模光与材质的交互？

---

## 9.1 光的物理学

### 9.1.1 粒子（Particles）

**爱因斯坦光量子假说**：光 = 离散的光子（packet of energy）。
- 光子能量 $E = h\nu$（$h$ = 普朗克常数，$\nu$ = 频率）
- **波动光学**在宏观上等价于**几何光学**

**应用**：高光、彩虹、薄膜干涉需波动光学；漫反射、镜面用几何光学足够。

### 9.1.2 介质（Medium）

- **真空**：无吸收
- **参与介质**（空气、水、烟）：有吸收和散射
- **吸收系数** $\sigma_a$、**散射系数** $\sigma_s$、**消光系数** $\sigma_t = \sigma_a + \sigma_s$

### 9.1.3 表面（Surfaces）

- 宏观表面（三角形）→ 微观表面（凹凸）
- 表面属性：粗糙度、折射率、金属度

### 9.1.4 次表面散射（Subsurface Scattering）

**问题**：光进入表面后多次散射再出射 → 皮肤、蜡、大理石效果。
**简化**：用 BRDF 近似（用经验参数），或用 BSSRDF（双向散射表面反射分布函数）。

---

## 9.2 相机（Camera）

**理想相机模型**：
- **针孔相机**：忽略透镜，简单但无景深
- **薄透镜模型**：有景深

**关键参数**：
- 焦距 $f$
- 光圈直径 $A$
- F 数 = $f/A$（小光圈大 F 数 → 大景深）
- 曝光时间 $\Delta t$
- ISO 感光度

**曝光值**：
$$EV = \log_2\left(\frac{f^2}{\Delta t \cdot ISO}\right)$$

---

## 9.3 The BRDF

**定义**：双向反射分布函数（Bidirectional Reflectance Distribution Function）。

$$
f_r(\mathbf{l}, \mathbf{v}) = \frac{dL_o(\mathbf{v})}{dE(\mathbf{l})}
$$

**物理意义**：从方向 $\mathbf{l}$ 入射的辐照度，对方向 $\mathbf{v}$ 出射的辐射度贡献。

**关键性质**：
- 亥姆霍兹互易性：$f_r(\mathbf{l}, \mathbf{v}) = f_r(\mathbf{v}, \mathbf{l})$
- 能量守恒：$\int f_r(\mathbf{l}, \mathbf{v}) (\mathbf{n} \cdot \mathbf{l}) d\mathbf{l} \le 1$
- 非负

**渲染方程**（辐射度形式）：
$$
L_o(\mathbf{p}, \mathbf{v}) = \int_{\Omega^+} f_r(\mathbf{l}, \mathbf{v}) L_i(\mathbf{p}, \mathbf{l}) (\mathbf{n} \cdot \mathbf{l}) d\mathbf{l}
$$

---

## 9.4 光照（Illumination）

**点光源**（无穷小、距离 $r$）：
$$E = \frac{I}{r^2}$$ （平方反比）

**面光源**（有限大小）：需要对光源面积积分（Ch10）。

---

## 9.5 菲涅尔反射（Fresnel）

**核心**：反射率随入射角变化——正视时反射弱，掠射时反射强（"水看起来底部更暗，远处更亮"）。

**完整公式**（不同偏振）：
- s 偏振：$r_s = \frac{n_1\cos\theta_i - n_2\cos\theta_t}{n_1\cos\theta_i + n_2\cos\theta_t}$
- p 偏振：$r_p = \frac{n_2\cos\theta_i - n_1\cos\theta_t}{n_2\cos\theta_i + n_1\cos\theta_t}$
- 反射率 $F = (r_s^2 + r_p^2)/2$

**Schlick 近似**（实时渲染标准）：
$$
F(\mathbf{h}, \mathbf{l}) \approx F_0 + (1 - F_0) (1 - \mathbf{h} \cdot \mathbf{l})^5
$$

**$F_0$ 值**：
- **电介质**（非金属）：$F_0 \approx 0.04$
- **金属**：$F_0$ = 金属颜色（金、银、铜等各有 RGB）

### 9.5.1 外反射
- 光从空气进入表面 → 反射

### 9.5.2 典型菲涅尔值
- 各种材质的 $F_0$ 常数（电介质 0.04，水 0.02，玻璃 0.04）
- 金属：金 (1.0, 0.71, 0.29), 银 (0.97, 0.96, 0.91), 铝 (0.91, 0.92, 0.92), 铜 (0.95, 0.64, 0.54)

### 9.5.3 内反射
- 光从表面内部反射回内部（用于透明材质）

---

## 9.6 微观几何（Microgeometry）

**核心**：宏观平面在微观尺度上有大量随机倾斜的小面。

- 法线 $\mathbf{n}$（宏观）= 所有微观法线的平均
- 不同微观面方向 → 不同的反射方向
- **NDF** $D(\mathbf{h})$：微观法线的分布密度

**关键**：NDF 的形状决定高光形状。

---

## 9.7 微表面理论（Microfacet Theory）

**核心假设**：每个微观面都是**完美镜面**，只有朝向不同。

$$
f_r(\mathbf{l}, \mathbf{v}) = \frac{F(\mathbf{h}, \mathbf{l}) \cdot G(\mathbf{l}, \mathbf{v}, \mathbf{n}) \cdot D(\mathbf{h}, \alpha)}{4 (\mathbf{n} \cdot \mathbf{l})(\mathbf{n} \cdot \mathbf{v})}
$$

**Cook-Torrance 三件套**：

| 符号 | 名称 | 物理意义 | 直观 |
|------|------|---------|------|
| $D(\mathbf{h})$ | 法线分布函数 NDF | 微观法线方向分布 | 高光形状 |
| $F(\mathbf{h}, \mathbf{l})$ | 菲涅尔 | 不同角度的反射率 | 边缘更亮 |
| $G(\mathbf{l}, \mathbf{v}, \mathbf{n})$ | 几何遮蔽 | 微面互相遮挡 | 自阴影 |

**为什么分母有 $4 (\mathbf{n} \cdot \mathbf{l})(\mathbf{n} \cdot \mathbf{v})$**？——雅可比变换（将微观空间映射到宏观空间）。

---

## 9.8 表面反射 BRDF 模型

### 9.8.1 法线分布函数 NDF

**各向同性 NDF**（沿 $\mathbf{n}$ 轴对称）：

**Beckmann**：
$$
D(\mathbf{h}) = \frac{1}{\pi \alpha^2 \cos^4\theta_h} \exp\left(-\frac{\tan^2\theta_h}{\alpha^2}\right)
$$
其中 $\theta_h = \arccos(\mathbf{n} \cdot \mathbf{h})$。

**GGX / Trowbridge-Reitz**（更现代，长尾）：
$$
D(\mathbf{h}) = \frac{\alpha^2}{\pi [(\mathbf{n} \cdot \mathbf{h})^2 (\alpha^2 - 1) + 1]^2}
$$

其中 $\alpha = \text{roughness}^2$（Disney 重映射，让滑块"线性"对应视觉效果）。

**各向异性 NDF**（沿两个主轴粗糙度不同）：

$$
D(\mathbf{h}) = \frac{1}{\pi \alpha_x \alpha_y \left( \frac{h_x^2}{\alpha_x^2} + \frac{h_y^2}{\alpha_y^2} + h_z^2 \right)^2}
$$

用于拉丝金属、头发等。

### 9.8.2 多次反弹的表面反射

单次微面反射只考虑第一次反射——实际中微观面间会**互相反射**（能量守恒修正）。

**Kulla-Conty 近似**：用一个额外的 BRDF 项补偿多次反弹。

---

## 9.9 次表面散射 BRDF

**与表面 BRDF 的区别**：光从入射点进入，在表面下散射后从不同点出射。

**简化模型**（仍按 BRDF 形式表达）：
- 漫反射 albedo
- 表面透射率
- 法线分布（与表面 BRDF 类似）

**物理参数**（PBRT 风格）：
- 散射系数 $\sigma_s$、吸收系数 $\sigma_a$
- 各向异性因子 $g$
- 平均自由程 $1/(\sigma_a + \sigma_s)$

**简化**（实时）：wrap lighting（半 Lambert）、预积分皮肤着色（Ch14）。

---

## 9.10 布料 BRDF

布料是**非传统材质**——纤维结构，不是连续表面。

**三种模型**：
1. **经验布料模型**：增加 sheen 项，模拟绒毛
2. **微表面布料**：用微面反射（纤维截面）
3. **微圆柱体布料**：用细长圆柱体表示纤维

**关键参数**：纤维密度、纤维朝向、纤维反射率。

---

## 9.11 波动光学 BRDF

**适用场景**：薄膜干涉（油膜、肥皂泡）、衍射（CD 光栅）。

### 9.11.1 衍射模型
- 微观周期性结构 → 衍射光栅
- 反射率随波长和方向变化（虹彩效果）

### 9.11.2 薄膜干涉
- 上下表面反射光干涉
- 厚度变化 → 颜色变化
- **公式**：$2 n d \cos\theta_t = m\lambda$（相长干涉）

---

## 9.12 分层材质

**适用**：漆面（透明漆 + 不透明金属）、皮肤（表皮 + 真皮 + 血液）。

**模型**：每层用 BRDF/BSSRDF 表达，层间用 Beer-Lambert 衰减。

**简化**：预计算为查找表。

---

## 9.13 混合和过滤材质

**问题**：混合不同粗糙度/颜色的材质时（如像素跨两种材质），避免混叠。

**方案**：
- **过滤法线**：将粗糙度高的微面求平均
- **过滤 NDF**：NDF 自身具有可过滤性
- **去相关**：不同材质不要"对齐"高光

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
        Vector3 F0;
    };

    // 漫反射 + 镜面 BRDF（Cook-Torrance）
    template<typename T>
    Color<T> evaluateBRDF(const BRDFInputs& mat,
                          const Vector<T, 3>& n,
                          const Vector<T, 3>& v,
                          const Vector<T, 3>& l);

    // GGX NDF
    template<typename T>
    T ggxNDF(const Vector<T, 3>& n, const Vector<T, 3>& h, T alpha);

    // Schlick 菲涅尔
    template<typename T>
    Color<T> schlickFresnel(const Color<T>& F0, T hDotL);
}
```

## 关键术语索引

| 术语 | 释义 |
|------|------|
| BRDF | 双向反射分布函数 |
| 微表面 | 微观尺度的随机倾斜小面 |
| NDF | 法线分布函数 |
| 菲涅尔 | 入射角相关的反射率 |
| Schlick 近似 | 菲涅尔的快速近似 |
| GGX | 流行的 NDF 形式 |
| 几何遮蔽 | 微面互相遮挡 |
| 互易性 | BRDF 的对称性 |
| 能量守恒 | BRDF 积分 $\le 1$ |
| $F_0$ | 法向入射反射率 |
| 多次反弹 | 微面间的互相反射 |
| Kulla-Conty | 能量补偿近似 |

## 关键公式（见 `formulas.md` §E）

- **Cook-Torrance BRDF**
- **Schlick 菲涅尔**
- **GGX NDF**
- **微表面能量守恒条件**

## 与 RTR4 其他章节关联

- **Ch5 Shading**：传统经验着色模型
- **Ch8 Light**：辐射度量学基础
- **Ch10 局部光照**：BRDF 的实际应用
- **Ch11 GI**：渲染方程
- **Ch14 体积**：次表面散射
- **Ch23 硬件**：BRDF 求值在 GPU 上
