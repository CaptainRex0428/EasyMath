# Ch26 Real-Time Ray Tracing 实时光线追踪（**核心参考**，在线内容）

**核心问题**：如何在 GPU 上做实时光线追踪？

---

## 26.1 光线追踪基础

**核心思想**：从相机像素发射射线，求与场景的最近交点。

**与光栅化的对比**：
- 光栅化：物体 → 像素（高效、复杂光照难）
- 光线追踪：像素 → 物体（直观、全局效果自然）

**Ray Tracing Pipeline**（DXR/Vulkan RT）：
1. **Ray Generation Shader**：每像素生成主射线
2. **Intersection Shader**：测试射线与图元
3. **Any-Hit / Hit Shader**：记录命中
4. **Closest-Hit / Miss Shader**：处理最近命中或未命中

---

## 26.2 光线追踪着色器

**RayGen**（必需）：
- 生成主射线
- 调用 TraceRay()

**Intersection**（可选用 BVH 替代）：
- 射线 vs 单个图元（三角形、AABB、过程式几何）

**Any-Hit**（可选）：
- 每次命中调用
- 用于透明度、alpha test

**Closest-Hit**：
- 最近命中调用
- 计算颜色、递归追踪

**Miss**：
- 未命中调用
- 天空盒 / 环境光

---

## 26.3 顶层和底层加速结构

**两层 BVH**：
- **顶层 BLAS**（Bottom Level）：每个模型一个 BVH
- **顶层 TLAS**（Top Level）：组合所有 BLAS（按实例）

**优点**：
- 模型可重用（多个实例）
- 动态对象支持
- 内存高效

---

## 26.4 一致性

**核心问题**：邻近射线往往命中邻近对象 → 冗余计算。

### 26.4.1 场景一致性

#### 空间数据结构的属性
- **空间一致性**（Spatial Coherence）：邻近射线访问相近节点
- **光线间一致性**（Inter-ray Coherence）
- **光线内一致性**（Intra-ray Coherence）：单根光线的访问路径
- **帧间一致性**（Temporal Coherence）：跨帧相似

#### 构造方案
- **SAH**：表面面积排序
- **宽 BVH 节点**：每节点多个物体 → 减少遍历
- **BVH 压缩**：缓存友好

#### 遍历方案
- **Mailbox**：包排序，重测检测
- **栈遍历**：手动栈（避免 GPU 栈溢出）
- **回溯** vs **前向** 遍历
- **鲁棒 BVH 遍历**：处理浮点错误

---

## 26.5 性能优化

**采样策略**：
- **重要性采样**：按 BRDF 概率分布采样
- **MIS**（Multiple Importance Sampling）
- **俄罗斯轮盘赌**：随机终止递归

**降噪**：
- 每像素 1-4 样本（噪声）
- **SVGF**（Spatiotemporal Variance Guided Filter）
- **NRD**（NVIDIA Real-time Denoiser）
- **OIDN**（Intel Open Image Denoise）

**时域累积**：
- **重投影**：把上一帧像素重投影到当前帧
- **运动向量**补偿
- **拒绝**：失败时重采样

---

## 26.6 软阴影（光线追踪）

**核心**：从光源发射阴影射线 → 测试是否被遮挡。

**软阴影射线**：
- 多根阴影射线模拟面光源
- 按光源面积采样

**算法**：
- 软阴影（PCSS in RT）
- 接触硬化

---

## 26.7 反射（光线追踪）

**反射射线**：
- 从交点发射
- 方向 = reflect(入射, 法线)

**递归**：
- 多次反射
- 俄罗斯轮盘赌终止

**降噪挑战**：
- 高频反射 → 噪声
- 时域累积

---

## 26.8 路径追踪

**完整 GI**：
- 像素 → 交点 → 漫反射 → 多个方向 → ...

**采样数**：
- 1-4 样本（实时）
- 大量样本（离线）

**降噪**：
- 必须
- GI 是低频 → 容易降噪

---

## 26.9 降噪（Denoising）

**输入**：低样本（噪声）图像 + 辅助缓冲（深度、法线、运动）

**算法**：
- **双边滤波**：保边平滑
- **联合双边滤波**：用深度/法线加权
- **时域累积**：跨帧混合
- **AI 降噪**：CNN 学习

**代表实现**：
- **SVGF**（Frostbite）
- **NRD**（NVIDIA）
- **OIDN**（Intel）

---

## 26.10 混合渲染（Hybrid Rendering）

**核心**：光栅化 + 选择性光线追踪。

**方案**：
- 主渲染：光栅化（快）
- 反射：光线追踪（自然）
- 阴影：光线追踪（精确）
- AO：光线追踪（细节好）
- GI：光线追踪 + 降噪

**代表**：
- **Unreal Engine 5 Lumen**
- **Unity HDRP**
- **Frostbite**

---

## 26.11 未来方向

- **Mesh Shaders**：更灵活的图元处理
- **DXR 1.2 / Vulkan RT 升级**
- **AI 加速 BVH 构建**
- **神经辐射缓存**
- **神经 BRDF**
- **真实感全 GI**

---

## 与 EasyMath 关联

**中**——EasyMath 是光线追踪的数学基础。

**未来可加**：
- AABB / Sphere 类（已在 Ch22 列出）
- Ray 类
- 求交算法（Möller-Trumbore 在 Ch22）

**EasyMath 待实现 API 概念**：
```cpp
// 概念性 API
namespace EasyMath {
    struct Ray { Vector3 origin, direction; };

    // 已在 Ch22 列出
    bool intersect(const Ray&, const AABB&, float& t);
    bool intersect(const Ray&, const Sphere&, float& t);
    bool intersect(const Ray&, const Triangle&, float& t);

    // 反射/折射
    Vector3 reflect(const Vector3& incident, const Vector3& normal);
    Vector3 refract(const Vector3& incident, const Vector3& normal, float eta);

    // 重要性采样
    Vector3 cosineSampleHemisphere(float u1, float u2);
    Vector3 ggxSample(float u1, float u2, float roughness);
}
```

## 关键术语

| 术语 | 释义 |
|------|------|
| RayGen | 主射线着色器 |
| Closest-Hit | 最近命中处理 |
| Any-Hit | 任何命中处理 |
| Miss | 未命中处理 |
| BLAS / TLAS | 底层/顶层加速结构 |
| MIS | 多重要性采样 |
| 俄罗斯轮盘赌 | 随机终止 |
| SVGF | 时空方差引导滤波 |
| NRD | NVIDIA 实时降噪 |
| Lumen | UE5 实时 GI |
| 一致性 | 相邻访问相似性 |
| Mailbox | 包排序 |
| 重要性采样 | 按概率分布采样 |

## 关键算法

- **BVH 遍历**（栈式）
- **Möller-Trumbore** 三角形求交
- **反射/折射方向**
- **降噪滤波**（双边、时域累积）

## 关键路径

- 完整路径追踪 = RayGen → Trace → Closest-Hit → 反射/折射 → 递归
- 降噪 = 路径追踪（低样本）+ 滤波 = 实时 GI
