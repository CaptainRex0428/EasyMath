# Ch7 Shadows 阴影

**核心问题**：如何让场景中的阴影渲染得又快又好？

---

## 7.1 平面阴影

### 7.1.1 投影阴影
将物体投影到平面上（用投影矩阵）。

**限制**：仅适用平面接收器。

### 7.1.2 软阴影
通过调整投影矩阵的边界 → 软化阴影边缘。

---

## 7.2 曲面上的阴影

- 球面阴影：在曲面上渲染
- **阴影纹理**：用阴影纹理模拟柔和过渡

---

## 7.3 阴影体算法（Shadow Volumes）

**核心**：从光源沿轮廓挤出（extrude）→ 形成"阴影体"。

**优点**：硬阴影精确
**缺点**：大量几何，依赖模板缓冲

**关键思想**：观察者在阴影体内 → 阴影。

---

## 7.4 阴影贴图（Shadow Map）

**核心**：从光源视角渲染深度到 shadow map，主相机比较深度。

**步骤**：
1. 以光源为相机，渲染深度图
2. 主相机渲染时，对每个片元做 lookup
3. 比较深度：片元深度 > shadow map 深度 → 阴影

**优点**：通用、GPU 友好
**缺点**：锯齿、走样、自阴影（shadow acne）、彼得潘（Pete-Panning）

### 7.4.1 分辨率增强
- **PCF**（Ch7.5）
- **VSM**（Variance Shadow Map）：用方差
- **ESM**（Exponential Shadow Map）：用指数
- **MSM**（Mipmap Shadow Map）：逐 mipmap 滤波

---

## 7.5 PCF（Percentage-Closer Filtering）

**核心**：对 shadow map 做区域采样（不是单点），按"通过/失败"比例决定阴影强度。

**标准实现**：
- $3 \times 3$ 或 $5 \times 5$ 采样
- 取"通过"比例 → 阴影强度
- 廉价软阴影

**缺点**：均匀软 → 物理上不对。

---

## 7.6 PCSS（Percentage-Closer Soft Shadows）

**核心**：根据**遮挡物到接收点距离**动态调整 PCF 半径。

**步骤**：
1. 找"平均遮挡者"距离
2. 距离小 → 半径小（锐利阴影）
3. 距离大 → 半径大（柔和阴影）

**优点**：真实感软阴影
**缺点**：多次采样，较贵

---

## 7.7 过滤阴影贴图

VSM、ESM、MSM 等方法，将"通过"测试转化为统计量。

---

## 7.8 体积阴影技术

- **不规则 Z-buffer**（Ch7.9）
- **阴影体积**（Ch7.3）
- **半影体积**（soft shadow volume）

---

## 7.9 不规则 Z-buffer

**核心**：每个像素存深度范围（min, max）→ 软阴影。

---

## 7.10 其他应用

- 阴影贴图用于其他效果：SSAO、反射、间接光照
- 软粒子

---

## 与 EasyMath 关联

**低**——EasyMath 是数学库，不涉及渲染管线。

**但是**：
- 阴影算法的数学（如 PCF 的滤波核、PCSS 的距离衰减）可能作为 Color.h 的工具函数
- 视锥体测试（Ch19）可能用于 light frustum

## 关键术语

| 术语 | 释义 |
|------|------|
| Shadow Map | 阴影贴图 |
| Shadow Volume | 阴影体 |
| PCF | Percentage-Closer Filtering |
| PCSS | PCF + 半影 |
| VSM | Variance Shadow Map |
| ESM | Exponential Shadow Map |
| Shadow Acne | 阴影表面上的自阴影 |
| Pete-Panning | 阴影与物体分离 |
| 自阴影 | 物体在自身表面上的阴影 |
| 接触硬化阴影 | PCSS 的视觉效果 |

## 关键算法对比

| 算法 | 质量 | 速度 | 适用 |
|------|------|------|------|
| Shadow Map | 低 | 快 | 简单场景 |
| PCF | 中 | 快 | 游戏 |
| PCSS | 高 | 慢 | 主机游戏 |
| Shadow Volume | 高（硬） | 中 | 工业级 |
| VSM/ESM | 中 | 快 | 移动 |
