# Ch22 Intersection Test Methods 相交测试方法（**核心参考**）

> **与 EasyMath 关联**：直接对应未来的 `AABB` / `OBB` / `Sphere` / `Ray` / `Triangle` / `Frustum` 类。
> **完成度**：P0（核心参考）— 完整抽取，2026-06-25

**核心问题**：如何高效判断几何对象是否相交？

**应用场景**：
- 拾取（点击 → 选 3D 对象）
- 碰撞检测（物理 + 游戏）
- 视锥体剔除（Ch19）
- 包围体测试（Ch19、Ch25）
- 光线追踪（Ch26）
- 距离查询

---

## 22.1 GPU 加速的拾取

### 22.1.1 CPU 拾取

- 像素 → 射线 → 场景 BVH 测试
- 框选 → 视锥体 → BVH 测试
- 缺点：位移/细分几何不匹配、alpha 贴图难处理

### 22.1.2 GPU 拾取（ID 渲染）

**核心思想**：
- 每个三角形/网格分配唯一 ID（作为颜色）
- 离屏渲染到小缓冲
- 鼠标点击 → 读像素颜色 → 知道选了哪个物体

**额外信息**：
- 顶点颜色（红绿蓝 = 重心坐标）→ 可算出三角形内相对位置
- 法线、UV 等都可存到屏幕外缓冲

**优化**：
- 静态场景：渲染一次复用
- 动态场景：用微缩窗口 + 离轴相机

---

## 22.2 定义和工具

### 22.2.1 射线

$$
R(t) = \mathbf{o} + t \mathbf{d}
$$

- $\mathbf{o}$：射线原点
- $\mathbf{d}$：单位方向向量（$||\mathbf{d}|| = 1$）
- $t < 0$：射线原点**后面**的点
- $t \ge 0$：射线上的点
- $t$ 是距离（因为 $\mathbf{d}$ 归一化）

**最大距离 $t_{\max}$**：从初始 $t_{\max} = \infty$ 开始，命中时更新为命中距离。射线变成**线段**。

### 22.2.2 隐式 vs 显式表面

**隐式表面**（公式 $f(\mathbf{p}) = 0$）：
- 例：球面 $\mathbf{p} \cdot \mathbf{p} = r^2$
- 优点：点在表面内/外容易判断
- 缺点：渲染需要光线步进或网格化

**显式表面**（参数化）：
- 例：三角网格、Bezier 曲面
- 优点：直接渲染
- 缺点：点在表面内/外需要测试

---

## 22.3 创建包围体（复习 Ch19）

### 22.3.1 AABB

$$
\text{AABB} = \{ \mathbf{p} : b_{min,i} \le p_i \le b_{max,i}, \forall i \in \{x, y, z\} \}
$$

**构建**（从三角形列表）：
```cpp
AABB build(const std::vector<Triangle>& tris) {
    AABB box;
    box.min = box.max = tris[0].v0;
    for (auto& t : tris) {
        box.min = min(box.min, t.v0, t.v1, t.v2);
        box.max = max(box.max, t.v0, t.v1, t.v2);
    }
    return box;
}
```

### 22.3.2 k-DOP（k-Discrete Oriented Polytope）

AABB 的推广：用 $k/2$ 对法线（固定方向）描述。

- $k = 6$：AABB
- $k = 14$：八边形柱体
- $k = 18$、$k = 26$：更紧的包围

### 22.3.3 包围球

```cpp
Sphere minSphere(const std::vector<Vector3>& points) {
    // Ritter 启发式
    Sphere box = computeBoundingBox(points);
    Sphere s;
    s.center = box.center;
    s.radius = 0;
    for (auto& p : points) {
        float d = length(p - s.center);
        if (d > s.radius) {
            s.center = (s.center + p) * 0.5f;
            s.radius = d * 0.5f;
        }
    }
    return s;
}
```

**更紧的球**：Welzl 随机算法（O(n) 期望时间）。

### 22.3.4 凸多面体

- **凸包**：从点集构建（如 QuickHull、Gift Wrapping）
- 紧致但构建成本高

### 22.3.5 OBB

- 中心 = 质心
- 轴 = PCA（主成分分析）
- 半轴长度 = 每轴最远投影

```cpp
OBB build(const std::vector<Vector3>& points) {
    Matrix3 cov = covariance(points);
    // 特征向量分解 → 轴
    auto [eigvals, eigvecs] = jacobi(cov);
    OBB o;
    o.center = mean(points);
    o.axes = eigvecs;  // 每列是一个轴
    for (auto& p : points) {
        Vector3 local = o.axes.transpose() * (p - o.center);
        o.halfExtents = max(o.halfExtents, abs(local));
    }
    return o;
}
```

---

## 22.4 几何概率

**Pick 法则**：随机点、随机三角形、随机光线-三角形相交。

**目的**：通过 Monte Carlo 评估算法效率。

---

## 22.5 经验法则

- **简单几何先测**（AABB 比 OBB 快）
- **利用一致性**（spatial / temporal）
- **BVH 优于 BSP**（更易构造、缓存友好）
- **Knee 平面测试**（视锥 vs AABB，6 个平面逐个测）

---

## 22.6 射线/球体相交

### 22.6.1 数学解法

**代入** $R(t)$ 到球方程：
$$
||\mathbf{o} + t\mathbf{d} - \mathbf{c}||^2 = r^2
$$

**展开**（设 $\mathbf{L} = \mathbf{o} - \mathbf{c}$）：
$$
a t^2 + b t + c = 0
$$
其中：
- $a = \mathbf{d} \cdot \mathbf{d} = 1$（$\mathbf{d}$ 归一化）
- $b = 2 \mathbf{d} \cdot \mathbf{L}$
- $c = \mathbf{L} \cdot \mathbf{L} - r^2$

**判别式** $\Delta = b^2 - 4ac$：
- $\Delta < 0$：无交点
- $\Delta = 0$：相切（一个交点）
- $\Delta > 0$：两个交点

**解**：
$$
t = \frac{-b \pm \sqrt{\Delta}}{2a}
$$
取较小正值（最近交点）。

### 22.6.2 优化解法

**几何解释**：球心到射线的最近距离平方 vs 半径平方。

设 $t_c = -\mathbf{d} \cdot \mathbf{L}$（球心在射线上的最近点的 $t$），
$d^2 = ||\mathbf{o} + t_c \mathbf{d} - \mathbf{c}||^2$。

- $d^2 > r^2$：无交点
- $d^2 = r^2$：相切
- $d^2 \le r^2$：相交，$t = t_c \pm \sqrt{r^2 - d^2}$

**优化**：避免完整平方展开，只算必要项。

```cpp
template<typename T>
bool raySphereIntersect(const Vector<T, 3>& o, const Vector<T, 3>& d,
                       const Vector<T, 3>& c, T r,
                       T& tOut) {
    Vector<T, 3> L = o - c;
    T tca = dot(L, d);
    T d2 = dot(L, L) - tca * tca;
    T r2 = r * r;
    if (d2 > r2) return false;
    T thc = std::sqrt(r2 - d2);
    T t0 = tca - thc;
    T t1 = tca + thc;
    if (t0 < 0) {
        if (t1 < 0) return false;  // 全在背后
        tOut = t1;
    } else {
        tOut = t0;
    }
    return true;
}
```

---

## 22.7 射线/Box 相交

### 22.7.1 平板法（Slab Method）

**核心**：每对轴对齐平面定义一个 "slab"（夹层），求射线与 slab 的相交区间。

**每轴 $i$**：
$$
t_{1,i} = \frac{b_{min,i} - o_i}{d_i},\quad t_{2,i} = \frac{b_{max,i} - o_i}{d_i}
$$

如果 $d_i = 0$：射线与轴平面平行 $\rightarrow$ 如果不在 slab 内则不相交。

**全局区间**：
$$
t_{enter} = \max_i \min(t_{1,i}, t_{2,i})
$$
$$
t_{exit} = \min_i \max(t_{1,i}, t_{2,i})
$$

**命中条件**：$t_{enter} \le t_{exit}$ 且 $t_{exit} \ge 0$（且 $t_{enter} \le t_{MAX}$）

**C++ 实现**（**slab method 经典**）：
```cpp
template<typename T>
bool rayAABBIntersect(const Vector<T, 3>& o, const Vector<T, 3>& d,
                      const Vector<T, 3>& bmin, const Vector<T, 3>& bmax,
                      T& tMin, T& tMax) {
    T tEnter = -std::numeric_limits<T>::infinity();
    T tExit = std::numeric_limits<T>::infinity();
    for (int i = 0; i < 3; i++) {
        if (std::abs(d[i]) < T(1e-8)) {
            // 平行于此轴
            if (o[i] < bmin[i] || o[i] > bmax[i]) return false;
        } else {
            T t1 = (bmin[i] - o[i]) / d[i];
            T t2 = (bmax[i] - o[i]) / d[i];
            if (t1 > t2) std::swap(t1, t2);
            tEnter = std::max(tEnter, t1);
            tExit = std::min(tExit, t2);
            if (tEnter > tExit) return false;
        }
    }
    if (tExit < 0) return false;  // 全在背后
    tMin = tEnter;
    tMax = tExit;
    return true;
}
```

### 22.7.2 射线斜率法

**核心**：避免除法（用倒数）。

**关键技巧**：预处理倒数 $1/d_i$，与 $t_{min,i}$、$t_{max,i}$ 比较。

适用：SIMD/批量求交。

---

## 22.8 射线/三角形相交

### 22.8.1 Möller-Trumbore 算法（**经典**）

**核心**：用克莱姆法则解线性方程组。

设：
- 三角形顶点 $\mathbf{v}_0, \mathbf{v}_1, \mathbf{v}_2$
- $\mathbf{e}_1 = \mathbf{v}_1 - \mathbf{v}_0$
- $\mathbf{e}_2 = \mathbf{v}_2 - \mathbf{v}_0$
- $\mathbf{s} = \mathbf{o} - \mathbf{v}_0$

**求解**（克莱姆法则）：
$$
\begin{pmatrix} t \\ u \\ v \end{pmatrix} = \frac{1}{\mathbf{d} \cdot (\mathbf{e}_1 \times \mathbf{e}_2)} \begin{pmatrix} \mathbf{s} \cdot (\mathbf{e}_1 \times \mathbf{e}_2) \\ \mathbf{d} \cdot (\mathbf{s} \times \mathbf{e}_2) \\ \mathbf{s} \cdot (\mathbf{d} \times \mathbf{e}_1) \end{pmatrix}
$$

**关键点**：
- 分母 $\mathbf{d} \cdot (\mathbf{e}_1 \times \mathbf{e}_2)$ 是**三角形法线点积射线方向**
- 若 $|分母| < \epsilon$：射线与三角形平行
- **命中条件**：$t \ge 0$，$u \ge 0$，$v \ge 0$，$u + v \le 1$

**重心坐标** $(u, v, 1-u-v)$：三角形内点的参数化。

### 22.8.2 C++ 实现

```cpp
template<typename T>
bool rayTriangleIntersect(
    const Vector<T, 3>& o, const Vector<T, 3>& d,
    const Vector<T, 3>& v0, const Vector<T, 3>& v1, const Vector<T, 3>& v2,
    T& t, T& u, T& v) {
    Vector<T, 3> e1 = v1 - v0;
    Vector<T, 3> e2 = v2 - v0;
    Vector<T, 3> pvec = cross(d, e2);
    T det = dot(e1, pvec);

    if (std::abs(det) < T(1e-8)) return false;  // 平行

    T invDet = T(1) / det;
    Vector<T, 3> s = o - v0;
    u = dot(s, pvec) * invDet;
    if (u < T(0) || u > T(1)) return false;

    Vector<T, 3> qvec = cross(s, e1);
    v = dot(d, qvec) * invDet;
    if (v < T(0) || u + v > T(1)) return false;

    t = dot(e2, qvec) * invDet;
    return t >= T(0);
}
```

### 22.8.3 优化版本

**预计算 $\mathbf{e}_1 \times \mathbf{e}_2$**（三角形法线）—— 静态三角形可一次算多次用。

**Watertight Ray-Triangle**（Woop 2013）：对退化三角形鲁棒，IEEE 浮点兼容。

**Plücker 坐标**：用 6D 坐标表示射线和三角形 → 一次点积替代多个叉积。

**Aila 2013**："分支优化的 SIMD 求交"。

---

## 22.9 射线/多边形相交

**策略**：
- 多边形 = 多个三角形
- 对每个三角形用 Möller-Trumbore
- 找最近命中

---

## 22.10 平面/Box 相交

**核心**：将 OBB 变换到 AABB 局部空间。

**步骤**：
1. 提取平面方程（点 + 法线）
2. 将平面变换到 box 局部空间
3. 用 AABB vs 平面测试

**AABB vs 平面**：
- 找"正顶点"（p-vertex，法线正方向最远点）
- 找"负顶点"（n-vertex，法线负方向最远点）
- p-vertex 距离 $< 0$ → 完全在负侧（不相交）
- n-vertex 距离 $> 0$ → 完全在正侧（不相交）
- 否则相交

---

## 22.11 三角形/三角形相交

**应用**：碰撞检测、CSG。

**Möller SAT 三角形算法**（1997）：
- 13 分离轴测试
  - 3 个三角形法线（$\mathbf{n}_1, \mathbf{n}_2, \mathbf{n}_1 \times \mathbf{n}_2$）
  - 9 个边-边叉积方向（$\mathbf{e}_1^i \times \mathbf{e}_2^j$）
- O(9) 比较

---

## 22.12 三角形/Box 相交

**应用**：布娃娃碰撞、毛发碰撞。

**算法**：
- 三角形 3 个顶点 vs box
- 三角形 3 条边 vs box
- 三角形面（平面）vs box
- 总 22 个 SAT 分离轴

---

## 22.13 BV/BV 相交

### 22.13.1 球-球

**最简单**：
$$
||\mathbf{c}_1 - \mathbf{c}_2|| \le r_1 + r_2
$$
**优化**：比较距离**平方**和 $(r_1 + r_2)^2$，避免 sqrt。

```cpp
bool sphereSphere(const Sphere& a, const Sphere& b) {
    T dist2 = dot(a.center - b.center, a.center - b.center);
    T r = a.radius + b.radius;
    return dist2 <= r * r;
}
```

### 22.13.2 球-Box

**算法**：
1. 找 box 上离球心最近的点 $\mathbf{p}$
2. 比较 $||\mathbf{p} - \mathbf{c}||$ 与 $r$

```cpp
Vector<T, 3> closestPointOnAABB(const Vector<T, 3>& p,
                                const Vector<T, 3>& bmin,
                                const Vector<T, 3>& bmax) {
    return clamp(p, bmin, bmax);
}

bool sphereAABB(const Sphere& s, const AABB& box) {
    Vector<T, 3> p = closestPointOnAABB(s.center, box.min, box.max);
    return length(s.center - p) <= s.radius;
}
```

### 22.13.3 AABB-AABB

**最简单**：
```cpp
bool aabbAABB(const AABB& a, const AABB& b) {
    if (a.max.x < b.min.x || a.min.x > b.max.x) return false;
    if (a.max.y < b.min.y || a.min.y > b.max.y) return false;
    if (a.max.z < b.min.z || a.min.z > b.max.z) return false;
    return true;
}
```

**6 次比较，提前退出**。

### 22.13.4 k-DOP/k-DOP

**$k$ 对法线** → $k$ 次 1D 测试（每对法线对应一个轴）。

```cpp
bool kDOPkDOP(const kDOP& a, const kDOP& b) {
    for (auto& axis : kDOP::axes) {
        if (a.max(axis) < b.min(axis)) return false;
        if (b.max(axis) < a.min(axis)) return false;
    }
    return true;
}
```

### 22.13.5 OBB/OBB

**SAT**（分离轴定理）：
- 3 个 A 轴（OBB A 的轴）
- 3 个 B 轴（OBB B 的轴）
- 9 个 A_i × B_j 叉积
- 共 15 分离轴

每轴测试：投影两个 OBB 到该轴，比较区间是否重叠。

---

## 22.14 视锥体相交测试

### 22.14.1 提取视锥体平面

**从 VP 矩阵提取 6 个平面**（Gribb-Hartmann 2001）：

```cpp
Frustum extractFrustum(const Matrix<T, 4, 4>& vp) {
    Frustum f;
    // 左
    f.planes[0] = {vp.row(3) + vp.row(0), vp(3,0) + vp(0,0)};
    // 右
    f.planes[1] = {vp.row(3) - vp.row(0), vp(3,0) - vp(0,0)};
    // 下
    f.planes[2] = {vp.row(3) + vp.row(1), vp(3,1) + vp(1,1)};
    // 上
    f.planes[3] = {vp.row(3) - vp.row(1), vp(3,1) - vp(1,1)};
    // 近
    f.planes[4] = {vp.row(3) + vp.row(2), vp(3,2) + vp(2,2)};
    // 远
    f.planes[5] = {vp.row(3) - vp.row(2), vp(3,2) - vp(2,2)};
    return f;
}
```

**注意**：行主序下用行操作；列主序下用列操作。

### 22.14.2 视锥体/球

```cpp
bool frustumSphere(const Frustum& f, const Sphere& s) {
    for (auto& p : f.planes) {
        T dist = dot(p.normal, s.center) + p.d;
        if (dist < -s.radius) return false;  // 完全在某平面外
    }
    return true;  // 在视锥内或相交
}
```

### 22.14.3 视锥体/Box（AABB）

```cpp
bool frustumAABB(const Frustum& f, const AABB& box) {
    for (auto& p : f.planes) {
        // p-vertex：法线正方向最远点
        Vector<T, 3> pv;
        for (int i = 0; i < 3; i++) {
            pv[i] = (p.normal[i] >= 0) ? box.max[i] : box.min[i];
        }
        if (dot(p.normal, pv) + p.d < 0) return false;  // 完全在外
    }
    return true;
}
```

---

## 22.15 线/线相交

### 22.15.1 二维

**方法 1**：解联立方程
$$
\mathbf{p}_1 + t(\mathbf{p}_2 - \mathbf{p}_1) = \mathbf{p}_3 + s(\mathbf{p}_4 - \mathbf{p}_3)
$$
解 $t, s$。

**方法 2**：用叉积判断共线性 + 求交点

### 22.15.2 三维

三维线段通常不共面（除非平行）→ 找**最短距离线**。

**核心**：6 个方程，求最小距离的两点。

---

## 22.16 三平面相交

**3 个平面** 交于一点（除非平行）：

$$
\mathbf{n}_i \cdot \mathbf{p} = d_i,\quad i = 1, 2, 3
$$

**求解**：线性方程组（克莱姆法则或 LU 分解）。

```cpp
Vector<T, 3> threePlaneIntersect(const Plane& p1, const Plane& p2, const Plane& p3) {
    Matrix<T, 3, 3> M{
        p1.normal[0], p1.normal[1], p1.normal[2],
        p2.normal[0], p2.normal[1], p2.normal[2],
        p3.normal[0], p3.normal[1], p3.normal[2]
    };
    Vector<T, 3> d{p1.d, p2.d, p3.d};
    return M.inverse() * d;
}
```

---

## 与 EasyMath 项目关联

### 待实现（v1.1+ 路线图）

**几何类**（基础）：
```cpp
namespace EasyMath {
    // 基础几何
    struct Ray { Vector3 origin, direction; };
    struct Triangle { Vector3 v0, v1, v2; };
    struct Plane { Vector3 normal; float d; };

    // 包围体
    struct AABB { Vector3 min, max; };
    struct OBB { Vector3 center, halfExtents; Matrix3 axes; };
    struct Sphere { Vector3 center; float radius; };
    struct Frustum { Plane planes[6]; };

    // 构造
    AABB buildAABB(const std::vector<Vector3>& points);
    OBB buildOBB(const std::vector<Vector3>& points);
    Sphere buildSphere(const std::vector<Vector3>& points);

    // 视锥体
    Frustum extractFrustum(const Matrix<float, 4, 4>& vp);

    // 相交测试（已在 §22 实现）
    bool intersect(const Ray&, const AABB&, float& tMin, float& tMax);
    bool intersect(const Ray&, const Sphere&, float& t);
    bool intersect(const Ray&, const Triangle&, float& t, float& u, float& v);
    bool intersect(const Ray&, const OBB&, float& t);
    bool intersect(const AABB&, const Frustum&);
    bool intersect(const Sphere&, const Frustum&);
    bool intersect(const AABB&, const AABB&);
    bool intersect(const Sphere&, const AABB&);
    bool intersect(const Sphere&, const Sphere&);

    // 距离查询
    float distance(const AABB&, const AABB&);
    Vector3 closestPointOnAABB(const Vector3& p, const AABB& box);
}
```

### 公式参考

- **Möller-Trumbore**：见 `formulas.md §I.1`
- **射线-AABB slab**：见 `formulas.md §I.2`
- **射线-球二次方程**：见 `formulas.md §I.3`

## 关键术语索引

| 术语 | 释义 |
|------|------|
| 射线 | $R(t) = \mathbf{o} + t\mathbf{d}$ |
| Möller-Trumbore | 经典射线-三角形 |
| Plücker 坐标 | 几何代数表示 |
| Slab 方法 | AABB 射线测试 |
| 分离轴定理 | SAT 通用测试 |
| 视锥体 | 6 平面定义 |
| Gribb-Hartmann | 视锥平面提取 |
| 视锥剔除 | 用视锥体测试包围体 |
| p-vertex / n-vertex | AABB 在平面正/负侧的角 |
| 退化三角形 | 面积为零的三角形 |
| 拾取 | 鼠标 → 选 3D 对象 |
| ID 渲染 | 拾取的 GPU 方法 |
| $t_{MAX}$ | 射线最大有效距离 |

## 关键公式

见 `formulas.md §I`：
- Möller-Trumbore 射线-三角形
- 平板法射线-AABB
- 二次方程射线-球

## 抽取历史

- 2026-06-25：基于 `RTR4/_raw/ch022.txt`（109KB）抽取
- 来源章节号：Ch22
- 跳过内容：无（完整抽取）
- 置信度：高
