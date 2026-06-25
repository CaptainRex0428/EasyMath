# Ch25 Collision Detection 碰撞检测（**核心参考**）

> **与 EasyMath 关联**：未来 GJK、SAT、CCD 等碰撞检测算法的参考。
> **完成度**：P0（核心参考）— 完整抽取，2026-06-25
> **注**：Ch25 来自 RTR4 在线内容（realtimerendering.com）

**核心问题**：如何高效地检测物体之间的碰撞？

---

## 25.1 宽阶段碰撞检测

**目的**：从大量对象中快速排除明显不碰撞的对。

### 25.1.1 扫描剪枝算法（Sweep and Prune, SAP）

**核心**：在某个轴上排序包围盒投影 → 找重叠对。

**步骤**：
1. 计算所有 AABB 在某轴（如 X）的 min/max
2. 按 min 排序
3. 扫描：与当前 min 重叠的所有 max 都是潜在对

**复杂度**：$O(n \log n + k)$，$k$ 是重叠对数。

**优化**：
- 多轴 SAP（X+Y 或 X+Y+Z）：减少漏检
- 插入排序（增量更新）
- 轴选择（最长轴效果最好）

```cpp
// 简化 SAP
struct Endpoint { float value; int id; bool isMin; };
std::vector<Endpoint> ends;
for (auto& box : boxes) {
    ends.push_back({box.min.x, box.id, true});
    ends.push_back({box.max.x, box.id, false});
}
std::sort(ends.begin(), ends.end());
std::set<int> active;
std::vector<std::pair<int, int>> pairs;
for (auto& e : ends) {
    if (e.isMin) {
        for (int id : active) pairs.push_back({id, e.id});
        active.insert(e.id);
    } else {
        active.erase(e.id);
    }
}
```

### 25.1.2 网格（Grid）

**核心**：把空间分格，每格存对象 → 只检查同格/邻格。

**均匀网格**：
- 格大小 = $\sqrt[3]{V/n}$（$V$ 体积，$n$ 物体数）
- 26 邻居（3D）
- 哈希：cell 坐标 → 哈希值

**层次网格**：根据密度调整格子大小。

### 25.1.3 层次包围体（宽阶段）

**核心**：粗 BVH 筛 → 中阶段 BVH 精确测。

---

## 25.2 中阶段碰撞检测

### 25.2.1 BVH 构建

**SAH 启发式**（见 Ch19.1.1）：
$$
C = C_T + \frac{S_A}{S_{total}} N_A C_i + \frac{S_B}{S_{total}} N_B C_i
$$

**构建算法**：
- 自顶向下（递归）—— 最常用
- 自底向上（合并）—— 质量更好但慢
- 增量构建（插入）

**质量**：
- SAH BVH 质量 > 中点切分 2-3 倍

### 25.2.2 BVH 间的碰撞测试

**递归遍历**：
```
testBVH(A, B):
    if (A.bv not intersect B.bv) return
    if (A and B are leaves):
        return narrow test(A.geom, B.geom)
    if (A is leaf):
        for child in B.children: testBVH(A, child)
    elif (B is leaf):
        for child in A.children: testBVH(B, child)
    else:
        // 选择较大 BVH 优先遍历
        testBVH(A.left, B); testBVH(A.right, B)
```

**优化**：
- 缓存（Mailbox）：避免重测同一对
- 启发式选择遍历顺序

### 25.2.3 BVH 成本函数

**SAH 代价**（构建时用）：
$$
C = C_T + \frac{S_A}{S_{total}} N_A C_i + \frac{S_B}{S_{total}} N_B C_i
$$

**运行时代价**（更新时用）：
$$
C = \sum_{i} T_{i,\text{test}} + N_{\text{false-positive}} T_{\text{narrow}}
$$

### 25.2.4 OBB 树

**OBB 选择**：
- 紧致但慢
- 适合旋转物体

**层次结构建立**：
- 自顶向下
- 分裂面（最长轴或 PCA）

**刚体运动处理**：
- 旋转物体 → OBB 需更新
- 层次重建代价高
- 替代：保存世界变换，每次测试时用

---

## 25.3 窄阶段碰撞检测

精确测试一对。

### 25.3.1 图元 vs 图元

- 三角形-三角形（Möller SAT）
- 三角形-AABB
- 球-球
- 球-平面
- 球-多边形

### 25.3.2 距离查询

**最近点**（avoiding actual collision）：
- GJK（Gilbert-Johnson-Keerthi）增量算法
- 用于物理约束、trigger volume

**穿透深度**（separating axis / Minkowski sum）：
- EPA（Expanding Polytope Algorithm）

---

## 25.4 射线碰撞检测

**核心**：找射线最近命中。

**应用**：武器弹道、玩家视线、AI 寻路。

**算法**：
- BVH 遍历
- Möller-Trumbore（每个三角形）
- Mailbox 加速

```cpp
RayHit rayCast(const Ray& ray, const BVH& bvh) {
    RayHit best = {INFINITY, -1};
    stack<BVHNode*> s;
    s.push(bvh.root);
    while (!s.empty()) {
        auto* n = s.top(); s.pop();
        if (!intersect(ray, n->bv, tMin, best.t)) continue;
        if (n->isLeaf) {
            for (int tri : n->triangles) {
                float t, u, v;
                if (rayTriangle(ray, tri, t, u, v) && t < best.t) {
                    best = {t, tri};
                }
            }
        } else {
            s.push(n->left);
            s.push(n->right);
        }
    }
    return best;
}
```

---

## 25.5 使用 BSP 树的动态 CD

**核心**：BSP 树加速射线投射和点查询。

**优势**：
- 静态场景极快
- 适合 FPS 游戏

**劣势**：
- 动态物体更新慢
- 已不流行

---

## 25.6 限时碰撞检测

**核心**：保证 CD 在固定时间内完成。

**方法**：
- **优先级队列**：先检测重要对
- **时间预算**：超时跳过
- **几何时间一致性**：重用上一帧结果
- **增量更新**：BVH 局部更新

---

## 25.7 可变形模型

**核心**：皮肤/布料/绳索的碰撞。

**方法**：
- 简化碰撞网格
- 凸包分解
- BVH 重建

---

## 25.8 连续碰撞检测（CCD）

**核心**：检测运动中的穿透（不只是终点位置）。

**应用**：
- 高速物体（子弹）
- 物理模拟稳定性（避免穿模）

**算法**：
- **保守前进**（Conservative Advancement）
- **Toi 求根**（Time of Impact）
- 找到最早接触时间

---

## 25.9 碰撞响应

**应用力**：

**法向冲量**（分离）：
$$
\mathbf{j}_n = \frac{-(1 + e)(\mathbf{v}_{\text{rel}} \cdot \mathbf{n})}{\frac{1}{m_1} + \frac{1}{m_2}} \mathbf{n}
$$

**摩擦**（切向）：
$$
\mathbf{j}_t = \min(|\mathbf{j}_t^{\text{desired}}|, \mu |\mathbf{j}_n|) \mathbf{t}
$$

**恢复系数**：
$$v_{after} = e \cdot v_{before}$$
其中 $e \in [0, 1]$（$e=0$ 完全非弹性，$e=1$ 完全弹性）。

**核心算法**：
- **冲量法**（Impulse-based）
- **Penalty 力**（弹簧惩罚）
- **约束求解器**（Sequential Impulse）

---

## 25.10 粒子

### 25.10.1 粒子系统
- 大量粒子（爆炸、烟）
- 粒子间碰撞可忽略
- 与场景碰撞（地面、墙）

### 25.10.2 粒子的物理模拟
- 流体（FEM、SPH、PIC/FLIP）
- 布料（mass-spring）
- 刚体（Newton-Euler）

---

## 25.11 动态相交测试

### 25.11.1 球-平面
- 当前距离 + 速度投影
- 球何时撞平面

### 25.11.2 球-球
- 相对速度 + 当前距离
- 求根

### 25.11.3 球-多边形
- 顶点 + 边 + 面
- 最近点
- SAT

### 25.11.4 动态分离轴方法（SAT）
- 考虑运动
- 时间相关测试

---

## 25.12 GJK 算法（关键）

**GJK（Gilbert-Johnson-Keerthi）**：求两凸包之间的最小距离。

**核心思想**：**Minkowski 差**。两凸包 $A - B = \{a - b : a \in A, b \in B\}$。原点是否在差集内 → 相交。

**支持点函数**：
$$
\text{support}(A, \mathbf{d}) = \arg\max_{\mathbf{p} \in A} \mathbf{d} \cdot \mathbf{p}
$$

**Minkowski 差支持点**：
$$
\text{support}_{A-B}(\mathbf{d}) = \text{support}_A(\mathbf{d}) - \text{support}_B(-\mathbf{d})
$$

**GJK 迭代**：构造一个包含原点的单纯形（2D 三角形 / 3D 四面体），不断用支持点替换最远点。当原点被包含 → 相交。

**关键性质**：
- $O(n)$ 单次迭代（$n$ 是顶点数）
- 通常 5-20 次迭代收敛

**EPA**（Expanding Polytope Algorithm）：当 GJK 检测到相交 → 找最近表面 → 穿透深度。

---

## 25.13 SAT（分离轴定理）

**核心定理**：两凸多边形不相交 ⇔ 存在一个分离轴使两者的投影不重叠。

**应用**：凸多边形-凸多边形、AABB-AABB、OBB-OBB 各种组合。

**步骤**：
1. 生成候选分离轴
   - 两多边形各自的法线
   - 两多边形各边两两叉积
2. 对每个分离轴，投影两多边形
3. 如果有任一轴投影不重叠 → 不相交
4. 所有轴都重叠 → 相交

**AABB-AABB SAT**：
- 3 个轴（X/Y/Z）
- 每轴独立比较 min/max

**OBB-OBB SAT**：
- 3 个 A 轴 + 3 个 B 轴 + 9 个 A_i × B_j 叉积
- 15 个分离轴

**三角形-三角形 SAT**：
- 6 个三角形法线
- 1 个叉积（$\mathbf{n}_1 \times \mathbf{n}_2$）
- 9 个边-边叉积
- 13 个分离轴

---

## 25.14 蒙皮/布料/绳索的碰撞

**简化碰撞网格**：
- 视觉网格很细（1000+ 三角形）
- 碰撞网格很粗（10-50 三角形）
- 每个顶点关联到碰撞网格的最近三角形

**凸包分解**：
- 任意网格 → 多个凸包
- 凸包之间做 SAT
- 凸包内部始终是凸的

**Spatial hashing**（布料/粒子）：
- 每帧重建空间哈希
- 邻域查询快

---

## 25.15 关键算法总结

| 算法 | 复杂度 | 适用 |
|------|--------|------|
| SAP | $O(n \log n + k)$ | 宽阶段 |
| 网格 | $O(n + k)$ | 均匀分布 |
| BVH | $O(\log n)$ | 通用 |
| GJK | $O(n)$ 每次迭代 | 凸包最小距离 |
| SAT | $O(n^2)$ 轴数 | 凸多边形 |
| Möller-Trumbore | $O(1)$ | 单三角形-射线 |

---

## 与 EasyMath 项目关联

### 待实现 API 概念

```cpp
// 概念性 API
namespace EasyMath {
    // 几何类（已在 Ch22 列出）

    // GJK 算法
    struct GJKResult {
        bool intersects;
        Vector3 closestPointA;
        Vector3 closestPointB;
        float distance;
    };
    GJKResult gjk(const ConvexHull& a, const ConvexHull& b);

    // SAT
    bool satConvex(const ConvexHull& a, const ConvexHull& b);
    Vector3 satPenetration(const ConvexHull& a, const ConvexHull& b);

    // EPA（穿透深度）
    float epaPenetration(const ConvexHull& a, const ConvexHull& b);

    // 射线-CD
    RayHit rayCast(const Ray& ray, const BVH& bvh);
    RayHit rayCast(const Ray& ray, const std::vector<Triangle>& tris);

    // 距离查询
    float distance(const AABB&, const AABB&);
    float distance(const Sphere&, const Triangle&);
    Vector3 closestPoint(const Triangle&, const Vector3&);

    // CCD
    float timeOfImpact(const Sphere& s, const Vector3& velocity, const Triangle& t);
}
```

### 公式参考

- **Minkowski 差**：`support_{A-B}(d) = support_A(d) - support_B(-d)`
- **GJK 收敛条件**：原点包含在单纯形内
- **恢复系数**：$e \in [0, 1]$

## 关键术语索引

| 术语 | 释义 |
|------|------|
| 扫描剪枝 | SAP，宽阶段 |
| BVH | 层次包围体 |
| SAH | 表面积启发式 |
| GJK | Gilbert-Johnson-Keerthi 距离 |
| EPA | Expanding Polytope Algorithm |
| SAT | 分离轴定理 |
| CCD | 连续碰撞检测 |
| 恢复系数 | 弹性碰撞参数 |
| 冲量 | 碰撞响应 |
| Minkowski 差 | GJK 的核心概念 |
| 支持点 | 凸包上最远点 |
| 单纯形 | 2D 三角形 / 3D 四面体 |
| 穿透深度 | 相交时的最小分离距离 |
| 冲量法 | 碰撞响应算法 |
| 顺序脉冲 | 约束求解器 |
| ToI | 首次接触时间 |
| Mass-spring | 布料模型 |
| FEM | 有限元方法 |
| SPH | 光滑粒子流体动力学 |

## 关键公式

### SAT 判定

对两凸多边形，存在分离轴 $a$ 使投影不重叠 → 不相交。

### GJK 收敛

迭代构造包含原点的单纯形 → 收敛 = 相交。

### Minkowski 差

$$
\text{support}_{A-B}(\mathbf{d}) = \text{support}_A(\mathbf{d}) - \text{support}_B(-\mathbf{d})
$$

### 恢复系数

$$v_{after} = e \cdot v_{before}$$

### 法向冲量

$$
\mathbf{j}_n = \frac{-(1+e)(\mathbf{v}_{rel} \cdot \mathbf{n})}{\frac{1}{m_1} + \frac{1}{m_2}} \mathbf{n}
$$

### 摩擦冲量

$$
|\mathbf{j}_t| \le \mu |\mathbf{j}_n|
$$

## 抽取历史

- 2026-06-25：基于 `RTR4/_raw/ch025.txt`（99KB）抽取
- 来源章节号：Ch25（在线内容）
- 跳过内容：Ch25 一些子小节示例
- 置信度：高
