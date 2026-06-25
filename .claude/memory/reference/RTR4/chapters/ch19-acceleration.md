# Ch19 Acceleration Algorithms 加速算法（**核心参考**）

> **与 EasyMath 关联**：与未来 `AABB` / `OBB` / `Sphere` / `Frustum` 类 + 加速算法 API 相关。
> **完成度**：P0（核心参考）— 完整抽取，2026-06-25

**核心问题**：如何让渲染系统处理大型场景？

**性能目标**：
- 更高帧率（60-90 FPS）
- 更高分辨率（4K-8K）
- 更真实材质
- **更多几何复杂度**（无上限）—— 必然需要加速

---

## 19.1 空间数据结构

### 19.1.1 空间数据结构概览

**核心思想**：层次化组织 → 查询从 $O(n)$ 降到 $O(\log n)$。

**主要类型**：

| 类型 | 空间细分 | 适用 | 特点 |
|------|---------|------|------|
| **BVH**（层次包围体） | 否 | 通用、动态 | 灵活、易构造 |
| **BSP 树** | 是 | 静态、刚体 | 严格前后排序 |
| **k-D 树** | 是 | 静态、点云 | 轴对齐 BSP |
| **八叉树** | 是 | 均匀分布 | 规则、3D 网格 |
| **场景图** | 否 | 层级关系 | 关注模型关系 |

**空间细分 vs 物体细分**：
- **空间细分**：BSP、八叉树——细分空间
- **物体细分**：BVH——细分物体

### 19.1.2 层次包围体（BVH）

**核心思想**：每个内部节点是子树的包围体。

**树结构**：
- 根节点 = 整个场景的 BV
- 内部节点 = 子节点 + 自己的 BV
- 叶子节点 = 实际几何

**构造**：
- 自下而上：合并相邻物体的 BV
- 自顶而下：递归切分

**优势**：
- 灵活（任意维度、任意形状）
- 动态更新容易
- 缓存友好（线性 BVH）

### 19.1.3 BVH 节点度数

- 二叉树：$k = 2$（最简单）
- $k = 4, 8$：性能更好（高度低）
- $k$ 越大 → 树越矮 → 遍历步骤少 → 但每节点工作多

**推荐**：
- 二叉树：$k = 2$
- 高扇出：$k = 4, 8$
- 构造简化：$k = 2$ 沿最长轴；$k = 4$ 沿两个最长轴；$k = 8$ 沿所有轴

### 19.1.4 BVH 遍历（光线求交）

```
rayBVH(ray, node, tMin, tMax):
    if (ray not hit node.bv) return NO_HIT
    if (node is leaf):
        return test ray vs node.geometry
    closest = +INF
    for each child in node.children:
        t = rayBVH(ray, child, tMin, tMax)
        if (t < closest):
            closest = t
    return closest
```

**最近距离剪枝**：跟踪当前最近命中 $t_{best}$，子节点命中 $> t_{best}$ → 跳过。

### 19.1.5 BVH 动态更新

**问题**：物体移动 → BV 可能不再包含子物体 → BVH 失效。

**解决方案**：
1. **移除 + 重新插入**：失效节点从树中删除，重新插入根
2. **递归扩展**：扩展父节点 BV 包含移动后的子节点
3. **时域包围体（TBV）**：基于运动极限
4. **自下而上重拟合**（refit）：不重建，只更新 BV
5. **部分重建**：选择最差的部分重建

### 19.1.6 BSP 树（二叉空间划分）

**两种 BSP**：
- **轴对齐 BSP（k-D 树）**：用轴对齐平面划分
- **多边形对齐 BSP**：用多边形所在平面划分

**轴对齐 BSP 构造**：
- AABB 包围场景
- 选一个轴 + 一个位置 → 划分平面
- 物体分类（与平面相交的小物体存到上层或切分）
- 递归直到停止条件

**优势**：粗略的前后排序（用于画家算法）。

**多边形对齐 BSP 构造**：
- 选一个多边形 → 分割平面
- 其他多边形按平面分类（相交的被切分）
- 递归

**优势**：**严格**的前后排序（画家算法用）。

**应用**：Doom（90 年代经典）、碰撞检测、静态场景。

### 19.1.7 八叉树

**核心**：每层把 box 同时沿 3 个轴在中心切分 → 8 个子 box。

**结构**：
- 根：场景包围盒
- 每层：8 倍
- 叶子：足够小或物体足够少

**优势**：
- 高度规则（内存布局简单）
- 邻居查询快（24/26 邻居）
- 适合均匀分布

**缺点**：
- 物体分布不均时浪费空间
- 解决方法：**松散八叉树**（loose octree，BV 重叠）

### 19.1.8 缓存无关与缓存感知

**核心问题**：BVH 树结构的内存布局 → 缓存命中率。

**线性 BVH（LBVH）**：
- 节点连续存储（DFS 顺序）
- 父节点 + 2 子节点连续 → 缓存友好

**Morton 曲线布局**：
- 用 Morton（Z-order）曲线排序节点
- 局部性极好

**H-PLOC**（Heuristic Pointer-rearranging LOC）：
- 运行时调整指针

### 19.1.9 场景图

**核心**：节点有几何 + 变换 + 子节点。

**作用**：
- 层级变换继承
- 包围体层级（用于剔除）
- 动画控制
- 选择性渲染

**与 BVH 区别**：场景图是**逻辑结构**，BVH 是**空间索引**。

---

## 19.2 剔除技术概览

**5 层剔除**（按开销从小到大）：

1. **对象级剔除**：整体不渲染（最粗）
2. **背面剔除**（Ch19.3）：O(1)/三角形
3. **视锥体剔除**（Ch19.4）：O(1)/对象（BV 测试）
4. **入口剔除**（Ch19.5）：场景分区
5. **细节剔除**（Ch19.6）：小三角形/小物体
6. **遮挡剔除**（Ch19.7）：被其他物体遮挡

**关键原则**：每层剔除只在前一层不通过时考虑；越早剔除越省。

---

## 19.3 背面剔除

**核心**：背向相机的三角形不渲染。

**方法**：
- 比较三角形法线与视线方向点积
- 顺时针/逆时针 vs 观察方向

**EasyMath 支持**：通过 `Vector` 的叉乘计算法线。

```cpp
Vector3 normal = cross(v1 - v0, v2 - v0);
bool isFrontFacing = dot(normal, viewDir) > 0;
```

---

## 19.4 视锥体剔除

**核心**：视锥外对象不渲染。

**算法**：包围体（BV）vs 视锥体（6 平面）。

**平面提取**（Gribb-Hartmann）：
```cpp
// 从 VP 矩阵提取 6 个平面（OpenGL 行主序）
Plane planes[6];
// 左
planes[0] = {vp.row(3) + vp.row(0)};
// 右
planes[1] = {vp.row(3) - vp.row(0)};
// 下、上、近、远 类似
```

**BV vs 视锥**：
- 球：简单距离比较
- AABB：用 p-vertex 思路（法线正方向最远点）
- OBB：投影后测试

```cpp
bool frustumAABB(const Frustum& f, const AABB& box) {
    for (auto& p : f.planes) {
        Vector3 pv;  // p-vertex
        for (int i = 0; i < 3; i++)
            pv[i] = (p.normal[i] >= 0) ? box.max[i] : box.min[i];
        if (dot(p.normal, pv) + p.d < 0) return false;
    }
    return true;
}
```

---

## 19.5 入口剔除（Portal Culling）

**核心**：通过"门"（如房间门口）连接空间区域。

**算法**：
- 找到相机所在区域
- 通过"门"找到相邻区域
- 递归找所有可见区域

**适用**：室内场景、关卡设计。

---

## 19.6 细节剔除和小三角形剔除

**核心**：屏幕覆盖小的对象/三角形不渲染。

**阈值**：< 10 像素 → 跳过。

**好处**：
- 减少 overdraw
- 避免像素着色器的浪费

**应用**：树、草、装饰物。

---

## 19.7 遮挡剔除（Occlusion Culling）

**核心**：被其他对象挡住的不渲染。

### 19.7.1 遮挡查询

**原理**：用上一帧深度或遮挡物 → 查询"我是否被挡住"。

**OpenGL 枚举**：
- `GL_SAMPLES_PASSED`：精确（可能慢）
- `GL_ANY_SAMPLES_PASSED`：粗略（快）

**DirectX 11+**：occlusion query 同理。

**异步性**：GPU 处理查询 → CPU 继续 → 后续读结果。

**谓词渲染**（predicated rendering）：
- 查询 + draw call 一起提交
- GPU 自动跳过被遮挡的 draw call

**优化策略**：
- 优先查询最可能被遮挡的物体
- 跟踪历史：常可见物体少查
- Mattausch 2016：批量查询、抖动采样

### 19.7.2 层次 Z 缓冲（HiZ）

**核心**：把深度图缩成层级。

**构造**：mipmap-like，每个 mip 层级存窗口内最远深度。

**查询**：
1. 物体包围盒投影到屏幕
2. 根据投影大小选 mip 层级
3. 比较包围盒最近深度 vs mip 深度
4. 若最近深度 > mip 深度 → 遮挡

**mip 层级选择**：
$$
\text{level} = \max\left(0, \min\left(\text{maxMipLevel}, \left\lfloor \log_2\left(\frac{\text{maxScreenSize}}{4} \right) \right\rfloor \right)\right)
$$

最多覆盖 4 个深度值（成本可预测）。

**球体最近深度**（观察空间）：
$$
t_{\text{closest}} = \text{clipZ} - r
$$

**双 pass 方法**（Haar & Aaltonen）：
- Pass 1：用上一帧 z-pyramid 剔除
- Pass 2：渲染所有"可能可见"物体，生成新 z-pyramid
- 对被剔除物体再测

**MHDB**（Masked Hierarchical Depth Buffer）：
- 每像素 3 bit（覆盖蒙版 + 2 个深度）
- 更好处理深度不连续

**商业实现**：
- Umbra（被广泛集成）
- Frostbite、寒霜、UE5

---

## 19.8 剔除系统

**多粒度组合**：
- Object 级别（粗）
- Cluster 级别（64-256 三角形/组）
- Triangle 级别（细）

**Cluster 剔除**：物体更小 → 更容易被剔除。

**Triangle 剔除**（GPU）：
- 透视除法后视锥
- 背面
- 退化
- 小三角形
- 遮挡

**间接绘制**（Open-draw indirect / DX execute indirect）：
- GPU 拉取压缩后的三角形列表

**经典实现**：
- Shopf 2008（GPU AI + 剔除）
- Haar & Aaltonen（《刺客信条：大革命》）
- Frostbite / 寒霜
- Engel 2007（可见性缓冲）

---

## 19.9 LOD（Level of Detail）

### 19.9.1 LOD 切换

#### 离散几何 LOD
- 多套不同精度网格
- 按距离切换
- **缺点**：popping（突然跳变）

#### 混合 LOD
- 切换时同时渲染两套
- 淡入淡出（alpha）
- Giegl & Wimmer 2007：先画 LOD1 → LOD2 淡入
- 缺点：渲染成本翻倍

#### Alpha LOD
- 物体随距离增加逐渐变透明
- 完全透明后停止渲染
- 简单、连续
- 缺点：仅完全消失后才省成本

#### CLOD（连续 LOD）
- 边折叠动态简化
- 无 LOD 模型库（运行时生成）
- 缺点：每帧重建 + 难并行

#### Geomorph LOD
- 顶点插值（LOD 之间）
- 与 CLOD 类似但更平滑
- Sander & Mitchell 2007

#### 分数曲面细分
- GPU 曲面细分因子为浮点
- 避免 popping
- 适用 Bezier、位移映射

### 19.9.2 LOD 选择

#### 基于范围
- 简单：距离 ∈ $[d_0, d_1]$ → LOD_n
- 缺点：远近不精确

**延迟切换**（hysteresis）：在边界附近添加缓冲，避免频繁切换。

#### 基于投影面积
- 物体在屏幕上的像素数
- 球体估计：
$$
\text{radius}_\text{projected} = \frac{r}{d_{\text{view}}} \cdot \frac{W/2}{\tan(\text{fov}/2)}
$$
- 像素数 ≈ $\pi r_\text{projected}^2$

#### 其他度量
- 速度（移动越快越粗糙）
- 屏幕像素密度
- 视差

### 19.9.3 限时 LOD 渲染

**核心**：保证帧率，必要时降级。

- GPU 时间预算
- 超过预算 → 切换到更低 LOD
- 雾 + LOD 组合

---

## 19.10 渲染大型场景

### 19.10.1 虚拟纹理
- 大纹理分块
- 按需加载到 GPU
- 适合开放世界

### 19.10.2 纹理转码
- GPU 上解压
- 节省 CPU 时间

### 19.10.3 通用流式传输
- 流式场景、纹理、几何
- LOD + 流式 → 大世界

### 19.10.4 地形渲染
- 块地形（chunk）
- 四叉树 LOD
- 视锥体剔除
- 高度图 + 噪声

---

## 19.11 SAH 公式

**SAH（Surface Area Heuristic）**—— BVH 构建核心：
$$
C_{\text{SAH}}(N, \text{split}) = C_{\text{traverse}} + \frac{S_A}{S_{\text{parent}}} \cdot N_A \cdot C_{\text{intersect}} + \frac{S_B}{S_{\text{parent}}} \cdot N_B \cdot C_{\text{intersect}}
$$

其中：
- $C_{\text{traverse}}$：节点遍历代价（常数）
- $S_A, S_B, S_{\text{parent}}$：A、B 子树、parent 的表面积
- $N_A, N_B$：A、B 子树中的图元数
- $C_{\text{intersect}}$：单个图元求交代价

**关键直觉**：表面积之比 = 该子树被光线/视锥击中的概率。

**优化 BVH 质量**：选择最小化 $\text{SAH}$ 代价的切分。

---

## 与 EasyMath 项目关联

### 待实现（v1.1+）

**几何类**：
- `AABB` / `OBB` / `Sphere`（已在 Ch22 列出）
- `Frustum`
- `Ray`
- `Plane`
- `Triangle`
- `KDTree<T>`（轴对齐 BSP 树）
- `Octree<T>`（八叉树）
- `BVH<T>`（层次包围体）

**剔除 API**：
```cpp
// 视锥体
Frustum extractFrustum(const Matrix<float, 4, 4>& vp);
bool intersect(const Frustum&, const AABB&);

// 遮挡查询接口（GPU）
struct OcclusionQuery {
    void begin();
    void end();
    bool isVisible();  // GPU 结果
};

// LOD
template<typename T>
T computeLODByRange(float distance, const std::vector<float>& thresholds);
T computeLODByProjectedArea(float projectedArea, const std::vector<float>& thresholds);
```

## 关键术语索引

| 术语 | 释义 |
|------|------|
| BVH | 层次包围体 |
| AABB | 轴对齐包围盒 |
| OBB | 定向包围盒 |
| k-D 树 | 轴对齐 BSP 树 |
| 八叉树 | 3D 网格 8 叉划分 |
| 视锥体 | 6 平面定义 |
| 入口剔除 | Portal Culling |
| 遮挡剔除 | Occlusion Culling |
| HiZ | 层次 Z 缓冲 |
| SAH | 表面积启发式 |
| popping | LOD 切换的视觉突变 |
| 混合 LOD | Geomorphing |
| 延迟切换 | Hysteresis |
| CLOD | 连续 LOD |
| 虚拟纹理 | 大纹理分块 |
| 流式传输 | 边下载边渲染 |
| 分数曲面细分 | 浮点细分因子 |
| 谓词渲染 | 查询+draw 合并 |
| Mortion 曲线 | Z-order 空间填充曲线 |

## 关键公式

- **SAH** 见上
- **mip 层级** 见 19.7.2
- **球体投影半径** 见 19.9.2

## 抽取历史

- 2026-06-25：基于 `RTR4/_raw/ch019.txt`（148KB）抽取
- 来源章节号：Ch19
- 跳过内容：Ch19.1 一些图示
- 置信度：高
