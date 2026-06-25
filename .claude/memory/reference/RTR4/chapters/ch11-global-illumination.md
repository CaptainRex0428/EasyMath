# Ch11 Global Illumination 全局光照

**核心问题**：如何让光在场景中"反射多次"产生真实感？

---

## 11.1 渲染方程

**完整形式**（包含自发光、多次反射、阴影）：
$$
L_o(\mathbf{p}, \mathbf{v}) = L_e(\mathbf{p}, \mathbf{v}) + \int_{\Omega} f_r(\mathbf{l}, \mathbf{v}) L_i(\mathbf{p}, \mathbf{l}) (\mathbf{n} \cdot \mathbf{l}) d\mathbf{l}
$$

**包含环境光遮蔽修正**：
$$
L_o(\mathbf{p}, \mathbf{v}) = L_e(\mathbf{p}, \mathbf{v}) + \int_{\Omega} \text{vis}(\mathbf{p}, \mathbf{l}) f_r(\mathbf{l}, \mathbf{v}) L_i(\mathbf{l}) (\mathbf{n} \cdot \mathbf{l}) d\mathbf{l}
$$

其中 $\text{vis}$ 是可见性函数（0 或 1）。

---

## 11.2 通用 GI

### 11.2.1 辐射度（Radiosity）

**核心**：将表面离散为 patch，假设漫反射 → 求解线性方程组。

**优点**：漫反射场景完美
**缺点**：O(n²) 内存，O(n³) 计算

**加速**：层次辐射度、多重网格

### 11.2.2 光线追踪

**核心**：从相机发射射线，递归追踪反射/折射。

**路径追踪**：Monte Carlo 估计渲染方程。
- **BRDF 采样**（重要性采样）
- **光源采样**（直接光）
- **MIS**（Multiple Importance Sampling）：组合两种采样

**降噪**：必须（每像素 1 样本会噪声）。

---

## 11.3 环境光遮蔽（AO）

**核心**：在角落/缝隙处减少环境光。

### 11.3.1 AO 理论
$$
\text{AO}(\mathbf{p}) = \frac{1}{\pi} \int_{\Omega^+} \text{vis}(\mathbf{p}, \mathbf{l}) (\mathbf{n} \cdot \mathbf{l}) d\mathbf{l}
$$
半球上未被遮挡方向的余弦积分。

### 11.3.2 可见性和 Obscurance
- **可见性**：是否被几何遮挡
- **Obscurance**：与距离相关的衰减

### 11.3.3 考虑相互反射
加入"反射 AO" → 多次反弹修正。

### 11.3.4 预计算 AO
- **AO 贴图**：bake 到顶点或纹理
- 静态场景适用

### 11.3.5 动态 AO
- **SSAO**（屏幕空间 AO）
- **HBAO**（Horizon-Based AO）
- **GTAO**（Ground Truth AO）
- **RTX AO**（光线追踪 AO）

### 11.3.6 屏幕空间方法
- 从深度/法线缓冲重建几何
- 在屏幕空间做半球采样

### 11.3.7 用 AO 着色
$$
L_{ao} = L_{env} \cdot \text{AO}
$$

---

## 11.4 定向遮蔽（Directional Occlusion）

**核心**：AO 是**标量**，但遮蔽应该是**方向性**的（如对天光的遮挡 vs 对地面的遮挡）。

**SSDO**（Screen Space Directional Occlusion）

---

## 11.5 漫反射 GI

### 11.5.1 表面预照明
预计算每个表面点的环境光遮蔽。

### 11.5.2 定向表面预照明
存储多个方向的环境光（用 SH）。

### 11.5.3 预计算传输
存储两个表面点之间的"传输"（光线传播）。

### 11.5.4 存储方法
- **Lightmap**：UV 空间贴图
- **顶点烘焙**：顶点属性
- **SH probe**：体素探针

### 11.5.5 动态漫反射 GI
- **LPV**（Light Propagation Volumes）：体素传播
- **RSM**（Reflective Shadow Maps）：从光源反射
- **Voxel Cone Tracing**：体素锥追踪
- **SSGI**（Screen Space GI）
- **DDGI**（Dynamic Diffuse GI）：Unity 的探针方案

### 11.5.6 光照传播体积
- 把场景体素化
- 在体素网格中迭代传播光

### 11.5.7 基于体素的方法
- 稀疏体素八叉树
- 每帧重建 + cone tracing

### 11.5.8 屏幕空间方法
- **SSGI**：从屏幕深度/法线重建 GI
- 限制：屏幕外无信息

### 11.5.9 其他
- **DDGI**（Dynamic Diffuse GI）
- **Probe-based**：场景中放置探针

---

## 11.6 镜面 GI

### 11.6.1 局部环境贴图
- **Reflection Capture**（UE）
- **Reflection Probe**（Unity）
- 在关键位置放置 cubemap

### 11.6.2 环境贴图动态更新
- 每帧更新关键探针
- 适合动态光源/物体

### 11.6.3 体素方法
- 用体素反射探针

### 11.6.4 平面反射
- 第二个相机渲染镜面
- 性能高、效果好
- 限制：仅平面

### 11.6.5 屏幕空间方法（SSR）
- 从屏幕深度/法线做 ray marching
- **接触硬化反射**
- 限制：屏幕外无反射

---

## 11.7 统一方法

**VXGI / GI-1.0**：把体素和 cone tracing 统一处理漫反射+镜面+多次反射。

---

## 与 EasyMath 关联

**低**——GI 是渲染器层面。

**但是**：
- SH 函数可能成为未来 EasyMath 的工具
- AO 公式中的余弦半球积分

## 关键术语

| 术语 | 释义 |
|------|------|
| GI | Global Illumination |
| 渲染方程 | 完整光照方程 |
| AO | 环境光遮蔽 |
| SSAO | 屏幕空间 AO |
| HBAO | 基于地平线 AO |
| GTAO | Ground Truth AO |
| RSM | 反射阴影图 |
| LPV | 光照传播体积 |
| SH | 球谐函数 |
| Lightmap | 光照贴图 |
| SSR | 屏幕空间反射 |
| Probe | 环境探针 |
| DDGI | 动态漫反射 GI |
