# Ch22 Intersection Test Methods 相交测试方法（**核心参考**）

**核心问题**：如何高效地判断几何对象是否相交？

---

## 22.1 GPU 加速的拾取

**拾取**：点击屏幕 → 选 3D 对象。

**实现**：
- 屏幕坐标 → 3D 射线
- 射线 vs 所有对象
- 找最近命中

**GPU 加速**：
- 用 transform feedback
- 计算着色器遍历
- 包围体树

---

## 22.2 定义和工具

- **射线** $R(t) = \mathbf{o} + t\mathbf{d}$
- **图元**：点、线、三角形、AABB、OBB、球
- **测试结果**：是/否 + 相交参数

---

## 22.3 创建包围体

### 22.3.1 AABB 和 k-DOP
- AABB：min/max
- k-DOP：固定方向的 $k/2$ 对法线

### 22.3.2 球
- **中心**：质心
- **半径**：到最远顶点的距离
- 紧密度：通常用迭代缩小

### 22.3.3 凸多面体
- **凸包**：从点集构建
- 紧致、但成本高

### 22.3.4 OBB
- **中心**：质心
- **轴**：PCA（主成分分析）
- **半轴长度**：每轴最远投影

---

## 22.4 几何概率

**Pick 法则**：随机向量、随机三角形、随机 ray-triangle 测试。

**目的**：比较算法效率。

---

## 22.5 经验法则

- **简单几何先测试**（AABB 比 OBB 快）
- **利用一致性**（spatial / temporal）
- **BVH 优于 BSP**（更易构造、缓存友好）

---

## 22.6 射线/球体相交

**数学解法**（一元二次）：
- 代入 $R(t)$ 到球方程
- 解 $at^2 + bt + c = 0$
- 判别式 $\Delta = b^2 - 4ac$

**优化解法**：
- 几何解法
- 避免 sqrt 直到必要时

**公式**（设 $L = \mathbf{o} - \mathbf{c}$）：
$$
a = \mathbf{d} \cdot \mathbf{d},\quad b = 2 \mathbf{d} \cdot L,\quad c = L \cdot L - r^2
$$

---

## 22.7 射线/Box 相交

### 22.7.1 平板法（Slab Method）

**核心**：每对平面定义一个"slab"，求射线与 slab 的相交区间，全局区间交集非空 → 命中。

**公式**（每轴 $i$）：
$$
t_{1,i} = \frac{b_{\min,i} - o_i}{d_i},\quad t_{2,i} = \frac{b_{\max,i} - o_i}{d_i}
$$
$$
t_{enter} = \max_i \min(t_{1,i}, t_{2,i}),\quad t_{exit} = \min_i \max(t_{1,i}, t_{2,i})
$$
命中：$t_{enter} \le t_{exit}$ 且 $t_{exit} \ge 0$。

### 22.7.2 射线斜率法

**核心**：避免除法，用倒数和整数比较。

---

## 22.8 射线/三角形相交

### 22.8.1 相交算法

**Möller-Trumbore**（标准）：
- 用重心坐标 + 克莱姆法则
- 无需预计算平面方程

**公式**（见 `formulas.md §I.1`）：
$$
\begin{pmatrix} t \\ u \\ v \end{pmatrix} = \frac{1}{\mathbf{d} \cdot (\mathbf{e}_1 \times \mathbf{e}_2)} \begin{pmatrix} \mathbf{s} \cdot (\mathbf{e}_1 \times \mathbf{e}_2) \\ \mathbf{d} \cdot (\mathbf{s} \times \mathbf{e}_2) \\ \mathbf{s} \cdot (\mathbf{d} \times \mathbf{e}_1) \end{pmatrix}
$$
其中 $\mathbf{e}_1 = \mathbf{v}_1 - \mathbf{v}_0$，$\mathbf{e}_2 = \mathbf{v}_2 - \mathbf{v}_0$，$\mathbf{s} = \mathbf{o} - \mathbf{v}_0$。

命中：$t \ge 0$，$u, v \ge 0$，$u + v \le 1$。

### 22.8.2 实现

**优化**：
- 预计算 $\mathbf{e}_1 \times \mathbf{e}_2$（面法线）
- **Watertight Ray-Triangle**：对退化三角形鲁棒
- **Plücker 坐标**：高效 SIMD

**C++ 片段**（EasyMath 待实现）：
```cpp
bool rayTriangleIntersect(
    const Vector<T, 3>& orig, const Vector<T, 3>& dir,
    const Vector<T, 3>& v0, const Vector<T, 3>& v1, const Vector<T, 3>& v2,
    T& tOut, T& uOut, T& vOut) {
    Vector<T, 3> e1 = v1 - v0;
    Vector<T, 3> e2 = v2 - v0;
    Vector<T, 3> pvec = cross(dir, e2);
    T det = dot(e1, pvec);
    if (std::abs(det) < epsilon) return false;
    T invDet = T(1) / det;
    Vector<T, 3> s = orig - v0;
    T u = dot(s, pvec) * invDet;
    if (u < 0 || u > 1) return false;
    Vector<T, 3> qvec = cross(s, e1);
    T v = dot(dir, qvec) * invDet;
    if (v < 0 || u + v > 1) return false;
    t = dot(e2, qvec) * invDet;
    return t >= 0;
}
```

---

## 22.9 射线/多边形相交

### 22.9.1 交叉点测试
- 用 Möller-Trumbore 测试每个三角形
- 找最近的

---

## 22.10 平面/Box 相交

### 22.10.1 AABB
- 平面 vs 每个轴的 slab
- 找 p-vertex 和 n-vertex

### 22.10.2 OBB
- 将平面变换到 OBB 局部空间
- 用 AABB 测试

---

## 22.11 三角形/三角形相交

**用途**：碰撞检测、CSG 渲染。

**算法**：
- Möller SAT 三角形
- 13 分离轴测试
- 分离轴：3 三角形的法线 + 9 边叉积

---

## 22.12 三角形/Box 相交

**用途**：布娃娃碰撞、毛发碰撞。

**算法**：
- 三角形 3 点 vs box
- 三角形边 vs box
- 三角形面 vs box
- 22 SAT 分离轴

---

## 22.13 BV/BV 相交

### 22.13.1 球-球
- 距离平方比较
- 简单

### 22.13.2 球-Box
- 找 box 上最近点
- 比较到球心距离

### 22.13.3 AABB-AABB
- 各轴独立比较
- 简单

### 22.13.4 k-DOP/k-DOP
- 每对法线独立比较
- $k$ 次 1D 测试

### 22.13.5 OBB-OBB
- SAT
- 15 分离轴

---

## 22.14 视锥体相交测试

### 22.14.1 提取视锥体平面
从 VP 矩阵提取 6 个平面。

### 22.14.2 视锥体/球
- 平面 vs 球心距离
- 简单

### 22.14.3 视锥体/Box
- 每个平面 vs box
- 简单 AABB 形式

---

## 22.15 线/线相交

### 22.15.1 二维
**方法 1**：解联立方程
**方法 2**：用叉积判断共线性

### 22.15.2 三维
- 最短距离线
- 求解 6 方程

---

## 22.16 三平面相交

**核心**：3 个平面求交点。

**唯一交点**（平行三平面时不存在）：
- 用线性方程组求解

---

## 与 EasyMath 关联

**高**——EasyMath 是相交测试的数学基础。

**未来可加**：
- AABB / OBB / Sphere 类
- 平面类（点 + 法线）
- 视锥体类（6 平面）
- 射线-三角形、Möller-Trumbore
- 视锥体-AABB 剔除
- BVH 构建算法

**EasyMath 待实现 API 概念**：
```cpp
// 概念性 API
namespace EasyMath {
    struct Plane { Vector3 normal; float d; };
    struct Ray { Vector3 origin, direction; };
    struct Triangle { Vector3 v0, v1, v2; };

    // 包围体类
    struct AABB { Vector3 min, max; };
    struct OBB { Vector3 center, halfExtents; Matrix3 axes; };
    struct Sphere { Vector3 center; float radius; };

    // 相交测试
    bool intersect(const Ray&, const AABB&, float& t);
    bool intersect(const Ray&, const Sphere&, float& t);
    bool intersect(const Ray&, const Triangle&, float& t);
    bool intersect(const Ray&, const OBB&, float& t);
    bool intersect(const AABB&, const Frustum&);
    bool intersect(const Sphere&, const Frustum&);

    // 视锥体构造
    Frustum extractFrustum(const Matrix<float, 4, 4>& vp);
}
```

## 关键术语

| 术语 | 释义 |
|------|------|
| Möller-Trumbore | 经典射线-三角形 |
| Plücker 坐标 | 几何代数表示 |
| Slab 方法 | AABB 射线测试 |
| 分离轴定理 | SAT 通用测试 |
| 视锥体 | 6 平面定义 |
| p-vertex / n-vertex | AABB 在平面正/负侧的角 |
| 退化三角形 | 面积为零的三角形 |

## 关键公式

- 射线-三角形（Möller-Trumbore）
- 射线-AABB（平板法）
- 射线-球（二次方程）
- 见 `formulas.md §I`
