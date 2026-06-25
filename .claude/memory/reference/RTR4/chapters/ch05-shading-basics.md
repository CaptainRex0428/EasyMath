# Ch5 Shading Basics 着色基础

**核心问题**：光线如何与表面交互产生最终颜色？

---

## 5.1 着色模型（Shading Model）

**从简单到复杂**：

| 模型 | 公式 | 特点 |
|------|------|------|
| Lambert | $I_d = k_d (\mathbf{n} \cdot \mathbf{l})$ | 漫反射 |
| Phong | $I_s = k_s (\mathbf{r} \cdot \mathbf{v})^\alpha$ | 高光 |
| Blinn-Phong | $I_s = k_s (\mathbf{n} \cdot \mathbf{h})^\alpha$ | 更高效 |
| Cook-Torrance | 见 Ch9 | 物理正确 |

**Phong 模型完整**：
$$
I = k_a I_a + \sum_{i} \left( k_d (\mathbf{n} \cdot \mathbf{l}_i) + k_s (\mathbf{r}_i \cdot \mathbf{v})^{\alpha} \right) I_i
$$

**C++ 片段**（EasyMath 待实现）：
```cpp
Color<T> phongShading(const Color<T>& lightColor,
                      const Vector<T, 3>& normal,
                      const Vector<T, 3>& lightDir,
                      const Vector<T, 3>& viewDir,
                      T kd, T ks, T shininess) {
    T diff = std::max(T(0), dot(normal, lightDir));
    Vector<T, 3> reflectDir = reflect(-lightDir, normal);
    T spec = std::pow(std::max(T(0), dot(reflectDir, viewDir)), shininess);
    return lightColor * (kd * diff + ks * spec);
}
```

---

## 5.2 光源

### 5.2.1 方向光
无穷远 → 光线方向恒定（适合太阳光）。

### 5.2.2 精确光源（Punctual Lights）
- **点光源**（泛光灯）：向各方向辐射，距离衰减
- **聚光灯**：圆锥限制
- **其他**：矩形光、圆盘光（属于面光源，Ch10）

**衰减函数**：
- 物理：$1/r^2$
- 经验：$1/(a + br + cr^2)$
- 修正：避免光源处无穷大

### 5.2.3 其他光源类型
- **面光源**（Ch10）
- **环境光**（IBL, Ch10）
- **HDRI 环境贴图**（Ch10）

---

## 5.3 实现着色模型

### 5.3.1 计算频率
- **顶点级**：法线、UV、光照方向插值
- **像素级**：最终颜色、纹理
- **频率太高会闪烁**（高光、走样）

### 5.3.2 实现示例
- BRDF 公式 → 着色器代码
- 纹理采样 → albedo / normal / roughness

### 5.3.3 材质系统
- 材质 = 着色器 + 参数 + 纹理
- 材质实例 vs 材质类

---

## 5.4 锯齿和抗锯齿

### 5.4.1 采样和滤波理论
**Nyquist 定理**：采样率 ≥ 2× 信号最高频率。

**重建**：理想低通滤波（sinc 函数）→ 无走样，但无限核。
**实用**：截断 sinc（mitchell、gaussian）→ 模糊 + 振铃。

**走样**：高频信号被欠采样 → 阶梯、闪烁、锯齿。

### 5.4.2 基于屏幕的抗锯齿

| 技术 | 原理 | 优缺点 |
|------|------|--------|
| **SSAA** | 超采样 → 降采样 | 慢、效果好 |
| **MSAA** | 几何边缘超采样 | 解决几何锯齿，**不解决着色锯齿** |
| **FXAA** | 后处理边缘检测 | 廉价、模糊 |
| **TAA** | 时域累积 + 抖动 | 接近电影级，需运动向量 |
| **SMAA** | 形态学 | 比 FXAA 质量好 |
| **DLSS / FSR** | AI 升频 | 现代技术，需特定硬件 |

---

## 5.5 透明度、Alpha、合成

### 5.5.1 混合顺序
- **画家算法**：从远到近绘制
- **Z-sort**：透明物体先排序
- **顺序无关透明度**（OIT）：复杂但解决任意顺序

### 5.5.2 顺序无关透明度算法
- **深度剥离**（Depth Peeling）
- **A-buffer**
- **Moment-based OIT**

### 5.5.3 Alpha 预乘
- 预乘 alpha：`color *= alpha`
- 合成时只需加法（避免乘法）

---

## 5.6 显示编码

最终输出从 Linear → sRGB。

**关键**：不要在 sRGB 空间做光照计算。

---

## 与 EasyMath 关联

**中**——EasyMath 暂无光照模型，但 Color.h 可作为材质颜色输入。

**未来**：
- Blinn-Phong 简单光照（入门示例）
- PBR BRDF（参考 Ch9）

## 关键术语

| 术语 | 释义 |
|------|------|
| Lambert | 漫反射 |
| Phong / Blinn-Phong | 高光模型 |
| MSAA | 多重采样抗锯齿 |
| TAA | 时域抗锯齿 |
| OIT | 顺序无关透明度 |
| Alpha 预乘 | 预乘 alpha 的合成 |
| 显示编码 | Linear → sRGB |
