# Ch26 Real-Time Ray Tracing 实时光线追踪（**核心参考**）

> **与 EasyMath 关联**：未来 BVH 遍历、ray-triangle 相交 API 的参考。
> **完成度**：P0（核心参考）— 完整抽取，2026-06-25
> **注**：Ch26 来自 RTR4 在线内容（realtimerendering.com）

**核心问题**：如何在 GPU 上做实时光线追踪？

---

## 26.1 光线追踪基础

**核心思想**：从相机像素发射射线，求与场景的最近交点。

**与光栅化的对比**：
- 光栅化：物体 → 像素（高效、复杂光照难）
- 光线追踪：像素 → 物体（直观、全局效果自然）

**Ray Tracing Pipeline**（DXR / Vulkan RT）：
1. **Ray Generation Shader**：每像素生成主射线
2. **Intersection Shader**：测试射线与图元
3. **Any-Hit / Hit Shader**：记录命中
4. **Closest-Hit / Miss Shader**：处理最近命中或未命中

---

## 26.2 光线追踪着色器

### 26.2.1 RayGen（必需）
- 生成主射线
- 调用 `TraceRay()`
- 处理最终颜色

### 26.2.2 Intersection（可选）
- 射线 vs 单个图元（三角形、AABB、过程式几何）
- 不用 BVH 时必须自己写

### 26.2.3 Any-Hit（可选）
- 每次命中调用
- 用于透明度、alpha test

### 26.2.4 Closest-Hit
- 最近命中调用
- 计算颜色、递归追踪

### 26.2.5 Miss
- 未命中调用
- 天空盒 / 环境光

---

## 26.3 顶层和底层加速结构

**两层 BVH**：
- **底层 BLAS**（Bottom Level）：每个模型一个 BVH
- **顶层 TLAS**（Top Level）：组合所有 BLAS（按实例）

**优点**：
- 模型可重用（多个实例）
- 动态对象支持
- 内存高效

**API**（DXR）：
```cpp
// 伪代码
RaytracingAccelerationStructure tlas, blas;
tlas.addInstance(blas);  // 多实例共享一个 BLAS
```

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
- **H-PLOC**：指针重排

#### 遍历方案
- **Mailbox**：包排序，重测检测
- **栈遍历**：手动栈（避免 GPU 栈溢出）
- **回溯** vs **前向** 遍历
- **鲁棒 BVH 遍历**：处理浮点错误（Aila 2013）

### 26.4.2 一致性利用

**1. 帧间一致性**：
- 重投影光线（time warp）
- 重用上一帧可见性

**2. 空间一致性**：
- 相邻像素用同一光线路径
- 屏幕空间分块

**3. 光线间一致性**：
- 同簇内光线共享遍历
- wavefront 调度

---

## 26.5 性能优化

### 26.5.1 采样策略

**重要性采样**：
- 按 BRDF 概率分布采样
- 漫反射 → 余弦加权
- 镜面 → GGX 分布采样

**MIS**（Multiple Importance Sampling）：
- 组合 BRDF 采样 + 光源采样
- 加权平衡（power heuristic）

**俄罗斯轮盘赌**：
- 概率 $p$ 继续
- $(1-p)$ 终止 + 加权
- 期望无偏

### 26.5.2 降噪（**关键**）

**输入**：低样本（噪声）图像 + 辅助缓冲（深度、法线、运动）

**算法**：
- **双边滤波**：保边平滑
- **联合双边滤波**：用深度/法线加权
- **时域累积**：跨帧混合
- **AI 降噪**：CNN 学习

**代表实现**：
- **SVGF**（Spatiotemporal Variance Guided Filter，Frostbite）
- **NRD**（NVIDIA Real-time Denoiser）
- **OIDN**（Intel Open Image Denoise）

**SVGF 算法**：
1. 第一帧：低样本渲染
2. 计算方差（atrous 小波）
3. 时间累积（带运动向量）
4. 边缘感知空间滤波

### 26.5.3 软阴影（光线追踪）

**核心**：从光源发射阴影射线 → 测试是否被遮挡。

**软阴影射线**：
- 多根阴影射线模拟面光源
- 按光源面积采样

**算法**：
- 软阴影（PCSS in RT）
- 接触硬化

### 26.5.4 反射（光线追踪）

**反射射线**：
- 从交点发射
- 方向 = reflect(入射, 法线)

**递归**：
- 多次反射
- 俄罗斯轮盘赌终止

**降噪挑战**：
- 高频反射 → 噪声
- 时域累积

### 26.5.5 路径追踪

**完整 GI**：
- 像素 → 交点 → 漫反射 → 多个方向 → ...

**采样数**：
- 1-4 样本（实时）
- 大量样本（离线）

**降噪**：
- 必须
- GI 是低频 → 容易降噪

---

## 26.6 混合渲染（Hybrid Rendering）

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
- **Frostbite**（寒霜）

---

## 26.7 未来方向

- **Mesh Shaders**：更灵活的图元处理
- **DXR 1.2 / Vulkan RT 升级**
- **AI 加速 BVH 构建**
- **神经辐射缓存**
- **神经 BRDF**
- **真实感全 GI**
- **硬件 RT 升级**：RTX 50 系列、AMD RDNA 4+

---

## 26.8 BVH 遍历算法

### 26.8.1 栈式遍历

```cpp
bool rayBVH(const Ray& ray, const BVHNode* root, RayHit& hit) {
    float tMin = 0, tMax = INFINITY;
    if (!intersect(ray, root->bv, tMin, tMax)) return false;

    Stack stack;
    stack.push({root, tMin, tMax});

    while (!stack.empty()) {
        auto [node, tEnter, tExit] = stack.pop();
        if (tEnter > hit.t) continue;  // 已找到更近的

        if (node->isLeaf) {
            for (auto& tri : node->tris) {
                float t, u, v;
                if (rayTriangle(ray, tri, t, u, v) && t < hit.t) {
                    hit = {t, tri, u, v};
                }
            }
        } else {
            float tLeftEnter, tLeftExit, tRightEnter, tRightExit;
            bool hitLeft = intersect(ray, node->left->bv, tLeftEnter, tLeftExit);
            bool hitRight = intersect(ray, node->right->bv, tRightEnter, tRightExit);

            // 优先访问更近的
            if (hitLeft && hitRight) {
                if (tLeftEnter < tRightEnter) {
                    stack.push({node->right, tRightEnter, tRightExit});
                    stack.push({node->left, tLeftEnter, tLeftExit});
                } else {
                    stack.push({node->left, tLeftEnter, tLeftExit});
                    stack.push({node->right, tRightEnter, tRightExit});
                }
            } else if (hitLeft) {
                stack.push({node->left, tLeftEnter, tLeftExit});
            } else if (hitRight) {
                stack.push({node->right, tRightEnter, tRightExit});
            }
        }
    }
    return hit.t < INFINITY;
}
```

### 26.8.2 Mailbox 优化

**核心问题**：同一三角形可能被多根光线测试。

**Mailbox**：
- 每帧给每个三角形一个时间戳
- 光线到达三角形前检查时间戳
- 时间戳更新 = "我测过了"
- 跳过已测三角形

**限制**：需要持久化状态（每光线 / 每帧重置）。

### 26.8.3 鲁棒 BVH 遍历（Aila 2013）

**问题**：浮点误差导致光线-边界错过。

**Aila 2013 解决**：
- 自适应 epsilon
- 鲁棒距离比较

---

## 26.9 降噪算法详解

### 26.9.1 SVGF（Spatiotemporal Variance Guided Filter）

**输入**：
- 低样本图像（1-4 spp）
- 法线、深度、运动向量
- 上一帧累积

**步骤**：
1. **方差估计**：邻域像素差异
2. **时间累积**：重投影到当前帧
   - 历史 = history * (1 - α) + current * α
   - α 取决于历史与新帧差异
3. **空间滤波**（atrous 小波）：
   - 边缘感知（用法线/深度差）
   - 多尺度 5x5 → 3x3 → ...

### 26.9.2 NRD（NVIDIA Real-time Denoiser）

**库函数**：
- `NRD_ReblurDiffuse`：漫反射 GI 降噪
- `NRD_ReblurSpecular`：镜面反射降噪
- `NRD_RelaxDiffuse` / `NRD_RelaxSpecular`：高质量版
- `NRD_SigmaShadow`：阴影降噪

**输入要求**：
- 法线��深度、roughness
- 运动向量
- 噪声信号 + 方差

### 26.9.3 OIDN（Open Image Denoise）

**Intel 开源降噪库**。
- CPU 端（高质量）
- 训练良好的 CNN

---

## 26.10 软阴影在 RT 中的实现

**多根阴影射线**（PCSS-RT）：
1. 从着色点到光源采样多根方向
2. 每根方向做 BVH 遍历
3. 统计可见比例

**MIS 阴影**：
- 1 根从光源的"主"射线 + N 根辅助射线
- 组合减少方差

---

## 与 EasyMath 项目关联

### 待实现 API 概念

```cpp
// 概念性 API
namespace EasyMath {
    // 几何类（已在 Ch22 列出）

    // BVH
    struct BVHNode {
        AABB bv;
        int left, right;  // 子节点索引，-1 表示叶子
        std::vector<int> triIndices;  // 叶子节点
    };
    struct BVH {
        std::vector<BVHNode> nodes;
        int root;
        // 构造
        static BVH build(const std::vector<Triangle>& tris, float saHeuristic = 1.0f);
    };

    // 光线追踪
    struct Ray {
        Vector3 origin, direction;
        float tMin = 0, tMax = INFINITY;
    };

    struct RayHit {
        float t = INFINITY;
        int triangleIndex = -1;
        Vector3 normal, uv;  // 重心坐标
    };

    // BVH 遍历
    bool rayCast(const Ray& ray, const BVH& bvh, RayHit& hit);

    // 光线生成
    Ray generatePrimaryRay(int x, int y, int width, int height,
                         const Matrix<float, 4, 4>& invVP);

    // 反射
    Vector3 reflect(const Vector3& incident, const Vector3& normal);
    Vector3 refract(const Vector3& incident, const Vector3& normal, float eta);

    // 重要性采样
    Vector3 cosineSampleHemisphere(float u1, float u2);
    Vector3 ggxSample(float u1, float u2, float roughness);
}
```

### 关键概念关联

- **BVH**：Ch19 已详述
- **Möller-Trumbore**：Ch22 已详述
- **BRDF 重要性采样**：Ch9 已详述

## 关键术语索引

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
| OIDN | Intel 降噪 |
| Lumen | UE5 实时 GI |
| 一致性 | 相邻访问相似性 |
| Mailbox | 包排序 |
| 重要性采样 | 按概率分布采样 |
| 降噪 | 减少 RT 噪声 |
| 混合渲染 | 光栅 + RT 混合 |
| 鲁棒遍历 | 抗浮点错误 |
| 时间累积 | 跨帧混合 |
| atrous 小波 | 多尺度空间滤波 |

## 关键算法

- **BVH 遍历**（栈式）
- **Möller-Trumbore** 三角形求交（Ch22）
- **反射/折射方向**
- **降噪滤波**（双边、时域累积）

## 关键路径

- 完整路径追踪 = RayGen → Trace → Closest-Hit → 反射/折射 → 递归
- 降噪 = 路径追踪（低样本）+ 滤波 = 实时 GI

## 抽取历史

- 2026-06-25：基于 `RTR4/_raw/ch026.txt`（86KB）抽取
- 来源章节号：Ch26（在线内容）
- 跳过内容：Ch26 一些硬件细节
- 置信度：高
