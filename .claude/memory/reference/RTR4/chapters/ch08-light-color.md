# Ch8 Light and Color 光与颜色（**核心参考**）

> **与 EasyMath 关联**：直接对应 `include/Color.h`。本章的 sRGB、Rec709、WCAG 对比度、辐照度等概念已被 `Color<T>` 实现。
> **完成度**：P0（核心参考）— 完整抽取，2026-06-25

**核心问题**：光是什么？颜色是什么？如何在计算机中量化和表示？

---

## 8.1 光量（Light Quantities）

### 8.1.1 辐射度量学（Radiometry）

**物理基础**——测量光的**能量**。

| 符号 | 名称 | 含义 | 单位 |
|------|------|------|------|
| $Q$ | 辐射能（Radiant energy） | 总能量 | J |
| $\Phi$ | 辐射通量（Radiant flux） | 单位时间能量 | W |
| $E$ | 辐照度（Irradiance） | 单位面积入射通量 | W/m² |
| $L$ | 辐射度（Radiance） | 单位面积单位立体角 | W/(m²·sr) |
| $I$ | 辐射强度（Radiant intensity） | 单位立体角 | W/sr |

**立体角**（Solid angle）：
$$
\Omega = \frac{A}{r^2} \text{（立体弧度 sr）}
$$

**球面度（steradian）**：球的总立体角 = $4\pi$ sr。

**关键关系**：
- 辐照度 $E$ 是辐射度 $L$ 在半球上的积分：
$$
E = \int_{\Omega^+} L \cos\theta \, d\omega
$$
- 辐照度按 $\cos\theta$ 衰减（**Lambert 余弦定律**）

**微分立体角**：
- 球坐标：$d\omega = \sin\theta \, d\theta \, d\phi$
- 余弦参数：$d\omega = d(\cos\theta) d\phi$

### 8.1.2 光度学（Photometry）

**人眼感知**的辐射度量——加权 $V(\lambda)$ 函数（眼睛对绿光最敏感）。

| 物理量 | 光度量 | 单位 |
|--------|--------|------|
| 辐射通量 $\Phi$ | 光通量 | 流明 lm |
| 辐照度 $E$ | 照度 | 勒克斯 lx = lm/m² |
| 辐射度 $L$ | 亮度 | 尼特 cd/m² |

**$V(\lambda)$**：CIE 标准观察者函数，峰值在 **555nm 绿光**。

**关键直觉**：辐射度量是"物理真实"，光度学是"人眼感知"。渲染可以任选其一（一般用物理量，最后做色调映射）。

**常见照度参考**：
- 阳光直射：~100,000 lx
- 室内照明：~500 lx
- 月光：~0.1 lx

### 8.1.3 色度学（Colorimetry）

**三色学说**：人眼有三种锥体细胞（S/M/L），任意颜色 = 三种基色的加权和。

**CIE XYZ 色度图**：把可见光波长映射到 (x, y) 色度坐标。

**关键概念**：
- **色温**（Color temperature）：黑体辐射温度
- **白点**（White point）：D65（6500K）是最常见的标准白
- **色域**（Color gamut）：设备能显示的颜色范围

**CIE 标准光源**：

| 光源 | 描述 | 色温 (K) |
|------|------|----------|
| A | 钨丝灯 | 2856 |
| D50 | 日光（早晨/傍晚） | 5000 |
| D55 | 日光（中午） | 5500 |
| **D65** | 日光（标准白） | **6500** |
| F 系列 | 荧光灯 | 2700-6500 |

**显示器白点**：sRGB 用 D65。

### 8.1.4 使用 RGB 颜色进行渲染

**问题**：三色学说 + 显示器三基色（红绿蓝磷光体）→ RGB 颜色空间。

**常用 RGB 空间**：

| 空间 | 白点 | 伽马 | 用途 |
|------|------|------|------|
| **sRGB** | D65 | ~2.2 | 显示器/网络（最常用） |
| Rec.709 | D65 | ~2.2 | HDTV |
| Rec.2020 | D65 | ~2.2 | UHDTV |
| Adobe RGB | D65 | ~2.2 | 印刷 |
| ProPhoto RGB | D50 | ~1.8 | 宽色域 |

**为什么物理计算用 Linear，存储用 sRGB？**
- 物理量（光强）线性叠加正确
- 显示器有伽马响应（CRT 物理 + 视觉感知联合导致），需要编码

**sRGB ↔ Linear 公式**：

**sRGB → Linear**（解码）：
$$
C_{\text{lin}} = \begin{cases} C_{\text{srgb}} / 12.92, & C_{\text{srgb}} \le 0.04045 \\ \left(\frac{C_{\text{srgb}} + 0.055}{1.055}\right)^{2.4}, & C_{\text{srgb}} > 0.04045 \end{cases}
$$

**Linear → sRGB**（编码）：
$$
C_{\text{srgb}} = \begin{cases} 12.92 \cdot C_{\text{lin}}, & C_{\text{lin}} \le 0.0031308 \\ 1.055 \cdot C_{\text{lin}}^{1/2.4} - 0.055, & C_{\text{lin}} > 0.0031308 \end{cases}
$$

**C++ 片段**（EasyMath `Color::toLinear` / `toSRGB`）：
```cpp
LinearColor<T> lin = srgb.toLinear();
sRGBColor<T> srgb = lin.toSRGB();
```

### 8.1.5 亮度（Luminance）

**物理定义**：单位立体角内的光通量。

**Rec.709 亮度**（HDTV）：
$$
Y_{709} = 0.2126 R + 0.7152 G + 0.0722 B
$$

**Rec.601 亮度**（NTSC/老电视）：
$$
Y_{601} = 0.299 R + 0.587 G + 0.114 B
$$

**Rec.2020 亮度**（UHDTV）：
$$
Y_{2020} = 0.2627 R + 0.678 G + 0.0593 B
$$

> EasyMath `Color::luminance(Rec709)` 等即此实现。

### 8.1.6 CIE XYZ 与 LMS

- **LMS**：长波、中波、短波锥体响应
- **XYZ**：CIE 1931 标准三色匹配函数

**转换矩阵**（CIE 1931 → linear sRGB）：
$$
\begin{pmatrix} R \\ G \\ B \end{pmatrix} = M_{\text{XYZ→sRGB}} \begin{pmatrix} X \\ Y \\ Z \end{pmatrix}
$$

具体矩阵略（标准值可在 CIE 文献查）。

---

## 8.2 从场景到屏幕

### 8.2.1 HDR 显示编码

**LDR**（Low Dynamic Range）：8 位/通道，$[0, 1]$。
**HDR**（High Dynamic Range）：10-16 位/通道，可超过 1（真实光强）。

**HDR 优势**：
- 自然光照下动态范围可达 **100,000:1**（远超 LDR）
- 物理正确的光照计算
- 后处理灵活（曝光调整、色调映射）

**HDR 显示标准**：
- **PQ 曲线**（Perceptual Quantizer，SMPTE ST 2084）：HDR10、Dolby Vision
- **HLG**（Hybrid Log-Gamma，BBC/NHK）：广播 HDR

### 8.2.2 色调映射（Tone Mapping）

**核心问题**：HDR 范围 → LDR 显示范围的转换。

**关键概念**：
- **色调映射**（Tone Mapping）：改变分布
- **色调再现**（Tone Reproduction）：模拟胶片/相机响应
- **曝光**（Exposure）：用光圈/快门/ISO 控制入射光量

**曝光公式**（手动控制）：
$$
L_{\text{out}} = f(E_{\text{scale}} \cdot L_{\text{in}})
$$
其中 $E_{\text{scale}}$ 是曝光因子，$f$ 是色调映射函数。

**常见色调映射算子**：

| 算子 | 公式/思路 | 特点 |
|------|----------|------|
| **Clamping** | $L_{\text{out}} = \min(L, 1)$ | 简单，会过曝 |
| **Reinhard** | $L_{\text{out}} = L / (1 + L)$ | 简单、自然 |
| **Reinhard 扩展** | $L_{\text{out}} = L \cdot (1 + L/L_{white}^2) / (1 + L)$ | 可调白点 |
| **Uncharted 2** | 复杂分段 | 电影感 |
| **ACES** | Academy Color Encoding | 业界标准 |
| **曝光调节** | $L_{in} \cdot k$ | 模拟相机光圈 |

**简化版 Reinhard**：
```cpp
Color<T> toneMap(const LinearColor<T>& hdr, T exposure = 1.0) {
    Color<T> c = hdr * exposure;
    return Color<T>{
        c[0] / (T(1) + c[0]),
        c[1] / (T(1) + c[1]),
        c[2] / (T(1) + c[2]),
        c[3]
    };
}
```

### 8.2.3 颜色分级（Color Grading）

**目标**：调整最终图像的色彩倾向（电影感、暖色/冷色、特定氛围）。

**常见操作**：
- 升降曝光
- 调整对比度（S 曲线）
- 调整饱和度
- 色调偏移（lift/gamma/gain）
- 查找表（LUT）

---

## 8.3 关键概念总结

### 8.3.1 辐照度 vs 辐射度

- **辐照度 $E$**：单位面积入射通量（与方向无关）
- **辐射度 $L$**：单位面积单位立体角的通量（与方向有关）

**BRDF 关系**（核心）：
$$
L_o(\mathbf{v}) = \int_{\Omega^+} f_r(\mathbf{l}, \mathbf{v}) L_i(\mathbf{l}) (\mathbf{n} \cdot \mathbf{l}) d\mathbf{l}
$$
输出辐射度 = 输入辐射度 × BRDF × 几何项的半球积分。

### 8.3.2 朗伯余弦定律

漫反射表面接收的光量按 $\cos\theta$ 衰减：
$$
E = L \cdot \cos\theta
$$

**重要**：这是 BRDF 中 $(\mathbf{n} \cdot \mathbf{l})$ 项的来源。

### 8.3.3 双向表面反射分布函数（BRDF）

定义（详见 Ch9）：
$$
f_r(\mathbf{l}, \mathbf{v}) = \frac{dL_o(\mathbf{v})}{dE(\mathbf{l})}
$$

### 8.3.4 渲染方程

$$
L_o(\mathbf{p}, \mathbf{v}) = L_e(\mathbf{p}, \mathbf{v}) + \int_{\Omega} f_r(\mathbf{l}, \mathbf{v}) L_i(\mathbf{p}, \mathbf{l}) (\mathbf{n} \cdot \mathbf{l}) d\mathbf{l}
$$

（详见 Ch9、Ch11）

---

## 8.4 渲染中的色彩管理

**正确流程**：
1. **线性空间**计算（光强、BRDF、积分）
2. **sRGB 编码**到帧缓冲
3. 显示器自动**sRGB 解码**显示

**常见错误**：在 sRGB 空间做光照计算 → ��理不正确（亮度过强）。

**EasyMath 处理**：
- `LinearColor` 用于物理计算
- `sRGBColor` 用于显示
- `luminance()` 自动转换到 linear 再计算
- `perceivedLuminance()` 直接在当前空间计算

---

## 8.5 辐射度量学完整公式

### 8.5.1 微分立体角

球坐标：
$$
d\omega = \sin\theta \, d\theta \, d\phi
$$

余弦参数 $\mu = \cos\theta$：
$$
d\omega = d\mu \, d\phi
$$

### 8.5.2 辐照度半球积分

$E$ 是 $L$ 在半球上的余弦加权积分：
$$
E = \int_{\Omega^+} L(\mathbf{l}) \cos\theta \, d\omega = \int_0^{2\pi} \int_0^{\pi/2} L(\theta, \phi) \cos\theta \sin\theta \, d\theta \, d\phi
$$

### 8.5.3 各向同性点光源辐照度

距离点光源 $r$ 处的辐照度：
$$
E = \frac{I}{r^2}
$$
其中 $I$ 是辐射强度。

---

## 8.6 颜色感知与差异

### 8.6.1 ΔE（色差）

**CIELAB ΔE76**（简单）：
$$
\Delta E = \sqrt{(\Delta L^*)^2 + (\Delta a^*)^2 + (\Delta b^*)^2}
$$
**CIEDE2000**（更准确）：复杂公式，考虑亮度、色调、饱和度。

### 8.6.2 视觉敏锐度

- 暗光下人眼主要用杆体细胞（无色）
- 亮光下用锥体细胞（彩色）
- 中间过渡区（中间视觉）

---

## 8.7 WCAG（无障碍）

**W3C WCAG 2.1 标准**：

| 等级 | 最小对比度（正文） | 大字体 |
|------|-------------------|--------|
| AA | 4.5:1 | 3:1 |
| AAA | 7:1 | 4.5:1 |

**公式**（Rec.709 相对亮度 + 0.05 偏移）：
$$
\text{contrast}(C_1, C_2) = \frac{L_1 + 0.05}{L_2 + 0.05},\quad L_1 \ge L_2
$$

> EasyMath `Color::contrastRatio()` 即此实现。

---

## 与 EasyMath 项目关联

### 已实现 ✓
- `sRGBColor<T>` / `LinearColor<T>` / `HSV` / `HSL` / `HSI` 类型
- sRGB ↔ Linear 转换（`toLinear` / `toSRGB`）
- RGB ↔ HSx 转换（`toHSV` / `toHSL` / `toHSI`）
- 亮度计算（`luminance` / `perceivedLuminance` / `relativeLuminance`）
- Rec.601/709/2020 标准
- WCAG 对比度（`contrastRatio`）

### 公式参考

- **Rec.709 亮度**：`Color::luminance()` 实现 $0.2126, 0.7152, 0.0722$
- **WCAG 对比度**：`Color::contrastRatio()` 实现
- **sRGB 伽马**：`Color::toLinear()` / `toSRGB()` 用分段函数

### 未来 PBR 渲染时需要的概念
- 辐射度 $L$ → 出射光
- 辐照度 $E$ → 入射光累加
- 立体角积分
- 色调映射（Reinhard / ACES）
- HDR / 16 位浮点帧缓冲
- 多个标准光源（A/D50/D65/F）
- 色彩空间转换（XYZ、ACEScg 等）

## 关键术语索引

| 术语 | 释义 |
|------|------|
| 辐射度量学 | 物理光测量（能量单位） |
| 光度学 | 人眼感知的光测量（流明） |
| 色度学 | 颜色匹配科学 |
| CIE XYZ | 标准色度空间 |
| 立体角 | 3D 角度（sr 球面度） |
| 辐照度 E | 单位面积入射功率 |
| 辐射度 L | 单位面积单位立体角功率 |
| sRGB | 标准 RGB（伽马编码） |
| Linear RGB | 物理线性 RGB |
| 色温 | 黑体辐射温度（K） |
| 色调映射 | HDR → LDR 转换 |
| WCAG | Web 内容无障碍指南 |
| ΔE | 颜色差异度量 |
| Rec.709/601/2020 | 不同的 RGB 编码标准 |
| HDR | 高动态范围 |
| LDR | 低动态范围 |
| 立体弧度 | sr，球面度 |
| $V(\lambda)$ | 人眼光谱响应函数 |
| 锥体细胞 | L/M/S 三种颜色感知细胞 |
| D65 | 标准白光（6500K） |

## 关键公式

### Rec.709 亮度
$$Y = 0.2126 R + 0.7152 G + 0.0722 B$$

### Rec.2020 亮度
$$Y = 0.2627 R + 0.678 G + 0.0593 B$$

### WCAG 对比度
$$CR = \frac{L_1 + 0.05}{L_2 + 0.05},\quad L_1 \ge L_2$$
（$L$ 是在 linear 空间用 Rec.709 计算的亮度）

### 朗伯余弦
$$E_{\text{surface}} = L_{\text{incident}} \cdot \max(0, \mathbf{n} \cdot \mathbf{l})$$

### 辐照度半球积分
$$E = \int_{\Omega^+} L(\mathbf{l}) (\mathbf{n} \cdot \mathbf{l}) d\omega$$

### sRGB → Linear
$$C_{\text{lin}} = \begin{cases} C_{\text{srgb}} / 12.92, & C \le 0.04045 \\ ((C + 0.055)/1.055)^{2.4}, & C > 0.04045 \end{cases}$$

### 渲染方程
$$L_o(\mathbf{p}, \mathbf{v}) = L_e + \int_{\Omega^+} f_r(\mathbf{l}, \mathbf{v}) L_i(\mathbf{l}) (\mathbf{n} \cdot \mathbf{l}) d\mathbf{l}$$

## 抽取历史

- 2026-06-25：基于 `RTR4/_raw/ch008.txt`（55KB）抽取
- 来源章节号：Ch8
- 跳过内容：Ch8.1 一些具体光源表
- 置信度：高
