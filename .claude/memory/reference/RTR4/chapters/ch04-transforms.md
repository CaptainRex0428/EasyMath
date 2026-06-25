# Ch4 Transform 变换（**核心参考章节**）

> **与 EasyMath 关联**：直接对应 `include/Matrix.h`、`include/Quaternion.h`、`include/EasyConversion.h`。本章节公式的 90% 已由 EasyMath 实现。

**核心问题**：如何用数学方法操控 3D 物体的位置、朝向、尺寸和形状？

---

## 4.0 仿射变换基础

**线性变换**：保持加法和标量乘法的变换 $L(\mathbf{a} + \mathbf{b}) = L(\mathbf{a}) + L(\mathbf{b})$。

**仿射变换**：线性变换 + 平移。形式为 $3\times 3$ 矩阵 + $3\times 1$ 平移向量，或 $4\times 4$ 齐次矩阵。

**仿射变换的关键性质**：
- 保持直线（直线变换后仍是直线）
- **保持平行性**（不保证角度和长度）
- 不一定保持手性（镜像会反转）

**实现**：4×4 齐次矩阵的最后一列存平移，最后一行存 $0, 0, 0, 1$。

---

## 4.1 基本变换

### 4.1.1 平移（Translation）

$$
T(t_x, t_y, t_z) = \begin{pmatrix} 1 & 0 & 0 & t_x \\ 0 & 1 & 0 & t_y \\ 0 & 0 & 1 & t_z \\ 0 & 0 & 0 & 1 \end{pmatrix}
$$

**C++ 片段**（EasyMath `MTXTranslation`）：
```cpp
auto M = MTXTranslation<float, 4>(tx, ty, tz);
// M(0,3) = tx, M(1,3) = ty, M(2,3) = tz
```

### 4.1.2 旋转（Rotation）

**2D 旋转矩阵**（推导：使用二角和差公式）：
$$
R(\theta) = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}
$$

**3D 绕坐标轴**（绕 $X$ 轴，$Y/Z$ 类似）：
$$
R_x(\theta) = \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & \cos\theta & -\sin\theta & 0 \\ 0 & \sin\theta & \cos\theta & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}
$$

**关键性质**：
- 行列式为 1
- 正交矩阵 $R^{-1} = R^T$
- 多个旋转矩阵的乘积仍是正交矩阵
- **逆矩阵 = 转置**（计算便宜）

**绕任意点旋转**（不是绕轴）：
$$M = T(\mathbf{p}) \cdot R(\theta) \cdot T(-\mathbf{p})$$

### 4.1.3 缩放（Scaling）

$$
S(s_x, s_y, s_z) = \begin{pmatrix} s_x & 0 & 0 & 0 \\ 0 & s_y & 0 & 0 \\ 0 & 0 & s_z & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}
$$

**均匀缩放**：$s_x = s_y = s_z$（保持形状）。
**非均匀缩放**：会改变法线方向（见 4.1.7）。
**反射**：负缩放因子 → 镜像翻转（会反转三角形顶点顺序 → 影响背面剔除）。

### 4.1.4 剪切（Shear）

将一个轴的位移作为另一个轴的函数：
$$
H_x = \begin{pmatrix} 1 & a & b & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix},\quad H_y = \begin{pmatrix} 1 & 0 & 0 & 0 \\ c & 1 & d & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}
$$

两个自由度的水平/垂直剪切：
$$
H = \begin{pmatrix} 1 & a & 0 & 0 \\ b & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}
$$

### 4.1.5 变换的连接（Concatenation）

**顺序很重要**——矩阵乘法不满足交换律。

**约定**：右乘先应用（post-multiplication），左乘后应用（pre-multiplication）。
- $M_{object} \cdot T$：先平移物体
- $T \cdot M_{object}$：后平移物体

**结合律可用**：$A \cdot B \cdot C = A \cdot (B \cdot C) = (A \cdot B) \cdot C$。

**性能优化**：组合多个静态变换为单个矩阵 → 减少运行时计算。

### 4.1.6 刚体变换（Rigid-body Transform）

只含平移和旋转的变换。**保持长度、角度、手性**。

任意刚体变换可写为：
$$
\mathbf{R} = \begin{pmatrix} \text{rotation} & \mathbf{t} \\ 0, 0, 0 & 1 \end{pmatrix}
$$

**逆矩阵**：
$$
\mathbf{R}^{-1} = \begin{pmatrix} R^T & -R^T \mathbf{t} \\ 0, 0, 0 & 1 \end{pmatrix}
$$

由于 $R^T$ 计算便宜，刚体逆矩阵计算量小。

### 4.1.7 法线变换（**重要**）

**为什么不能直接用模型矩阵变换法线？**

非均匀缩放或剪切会改变法线和切线之间的夹角——如果用同一个矩阵变换法线，就会指向错误方向。

**正确做法**：使用原始变换矩阵的**逆转置**（adjoint transpose）：

$$
\mathbf{n}' = (M^{-1})^T \mathbf{n}
$$

或等价地：伴随矩阵 $M^*$ 的转置。

**实用建议**：
- 纯旋转/平移（无缩放）时，$M^{-T} = M$，所以可以用模型矩阵
- 变换后**记得归一化**（长度可能改变）
- 如果从变换后的三角形边向量叉乘重算法线，就不需要这个修正

> **EasyMath 待办**：当前没有 `Vector.normalTransform(M)` 这样的方法。未来 PBR/光照计算时需要。

### 4.1.8 计算逆矩阵

三种方法（信息不同则方法不同）：
1. **已知类型**（平移/旋转/缩放/刚体）→ 闭式求逆（见 4.1.1-4.1.6）
2. **通用 4×4** → 伴随矩阵法（$\det(M) \cdot M^{-1} = \text{adj}(M)$，但贵）
3. **数值方法**（LU/Gauss-Jordan）→ O(n³)

**性能提示**：如果逆矩阵只用于变换向量，**只算左上 3×3 部分的转置**即可。

---

## 4.2 特殊矩阵变换

### 4.2.1 欧拉变换（Euler Transform）

**核心思想**：用三个角度（heading/pitch/roll 或 yaw/pitch/roll）描述一个方向。

**24 种组合**：6 个两轴旋转 × 2（内旋/外旋）= 24。

**典型顺序**：先 heading（绕 $Y$），再 pitch（绕 $X$），最后 roll（绕 $Z$）。

**矩阵形式**（heading $h$, pitch $p$, roll $r$）：
$$
M = R_z(r) \cdot R_x(p) \cdot R_y(h)
$$

### 4.2.2 从欧拉变换中提取参数

**步骤**：
1. 直接读出平移（最后一列前 3 行）
2. 从左上 $3\times 3$ 反解 3 个角度

**特殊处理**：万向锁（gimbal lock）—— 当某个中间角度为 $\pm 90°$ 时，两轴重合，丢失一个自由度。

**避免万向锁**：使用四元数或旋转矩阵表示。

### 4.2.3 矩阵分解

**目标**：从总变换矩阵 $M$ 提取子变换（平移、旋转、缩放、剪切）。

**方法**：
- 平移：直接读最后一列
- 旋转：从左上 $3\times 3$ 用 Gram-Schmidt 正交化
- 缩放：$M$ 的列向量的长度
- 剪切：从正交化前后差异得到

**应用**：
- 撤销已有变换
- 提取数据（如相机位置/朝向）
- 重用组件

### 4.2.4 绕任意轴旋转——Rodrigues 公式

**输入**：单位轴 $\hat{\mathbf{k}} = (k_x, k_y, k_z)$，角度 $\theta$。

**Rodrigues 紧凑形式**：
$$
R = I + \sin\theta \cdot K + (1-\cos\theta) \cdot K^2
$$

其中 $K$ 是 $\hat{\mathbf{k}}$ 的**反对称矩阵**（叉积矩阵）：
$$
K = \begin{pmatrix} 0 & -k_z & k_y \\ k_z & 0 & -k_x \\ -k_y & k_x & 0 \end{pmatrix}
$$

**展开形式**（直接构造 4×4 矩阵）：
$$
R(\hat{\mathbf{k}}, \theta) = \begin{pmatrix} c + k_x^2(1-c) & k_x k_y(1-c) - k_z s & k_x k_z(1-c) + k_y s & 0 \\ k_y k_x(1-c) + k_z s & c + k_y^2(1-c) & k_y k_z(1-c) - k_x s & 0 \\ k_z k_x(1-c) - k_y s & k_z k_y(1-c) + k_x s & c + k_z^2(1-c) & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}
$$

**C++ 片段**（EasyMath 待实现）：
```cpp
Matrix<T, 4, 4> MTXRotationAxis(const Vector<T, 3>& axis, T angle) {
    Vector<T, 3> k = axis.normalized();
    T c = std::cos(angle), s = std::sin(angle);
    T omc = T(1) - c;
    // 使用 4.2.4 公式填矩阵
}
```

---

## 4.3 四元数

### 4.3.1 数学背景

**定义**（Hamilton 1843）：
$$
\mathbf{q} = q_w + q_x \mathbf{i} + q_y \mathbf{j} + q_z \mathbf{k}
$$
其中 $\mathbf{i}^2 = \mathbf{j}^2 = \mathbf{k}^2 = \mathbf{ijk} = -1$。

**基本运算**：

**加法**（逐分量）：
$$
\mathbf{q}_0 + \mathbf{q}_1 = (w_0 + w_1,\ x_0 + x_1,\ y_0 + y_1,\ z_0 + z_1)
$$

**标量乘法**：
$$
s\mathbf{q} = (sw, sx, sy, sz)
$$

**乘法**（使用点积+叉积）：
$$
\mathbf{q}_0 \mathbf{q}_1 = (w_0 w_1 - \mathbf{v}_0 \cdot \mathbf{v}_1,\ w_0 \mathbf{v}_1 + w_1 \mathbf{v}_0 + \mathbf{v}_0 \times \mathbf{v}_1)
$$
其中 $\mathbf{q} = (w, \mathbf{v})$，$\mathbf{v} = (x, y, z)$。

**共轭**：
$$
\mathbf{q}^* = (w, -\mathbf{v}) = (w, -x, -y, -z)
$$

**范数**：
$$
||\mathbf{q}|| = \sqrt{w^2 + x^2 + y^2 + z^2} = \sqrt{\mathbf{q} \cdot \mathbf{q}^*}
$$

**逆**（对非零四元数）：
$$
\mathbf{q}^{-1} = \mathbf{q}^* / ||\mathbf{q}||^2
$$

**单位四元数**（$||\mathbf{q}|| = 1$）：$\mathbf{q}^{-1} = \mathbf{q}^*$。

**指数**（单位四元数）：
$$
\mathbf{q} = (\cos\phi, \sin\phi \mathbf{u}_q)
$$

表示绕单位轴 $\mathbf{u}_q$ 旋转 $2\phi$。

**对数**：
$$
\log(\mathbf{q}) = (0, \phi \mathbf{u}_q)
$$

**幂**：
$$
\mathbf{q}^t = (\cos t\phi, \sin t\phi \mathbf{u}_q)
$$

### 4.3.2 四元数变换

**旋转向量**（核心公式）：
$$
\mathbf{v}' = \mathbf{q} \mathbf{v} \mathbf{q}^{-1}
$$
其中 $\mathbf{v}$ 当作纯四元数 $(0, v_x, v_y, v_z)$。

**旋转构造**：$\mathbf{q} = (\cos(\theta/2), \sin(\theta/2) \hat{\mathbf{k}})$。

**C++ 片段**（EasyMath `Quaternion`，待实现 v1.1）：
```cpp
QuaternionF q{cosHalfTheta, sinHalfTheta * kx,
              sinHalfTheta * ky, sinHalfTheta * kz};
Vector<T, 3> rotated = q.rotate(vec);
```

#### 矩阵转换

四元数 → 旋转矩阵（**无三角函数**）：
$$
M = \begin{pmatrix} 1-2y^2-2z^2 & 2xy-2wz & 2xz+2wy & 0 \\ 2xy+2wz & 1-2x^2-2z^2 & 2yz-2wx & 0 \\ 2xz-2wy & 2yz+2wx & 1-2x^2-2y^2 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}
$$

**反向**（矩阵 → 四元数）：使用迹区分情况，**数值稳定性提示**：
- 先选最大对角元 $M_{ii}$ 作为基
- 避免除以过小数值

#### 球面线性插值（Slerp）

**核心问题**：在两个旋转之间平滑过渡。

**数学公式**：
$$
\text{slerp}(\mathbf{q}_0, \mathbf{q}_1; t) = \frac{\sin((1-t)\Omega)}{\sin\Omega} \mathbf{q}_0 + \frac{\sin(t\Omega)}{\sin\Omega} \mathbf{q}_1
$$

其中 $\cos\Omega = \mathbf{q}_0 \cdot \mathbf{q}_1$（点积）。

**软件实现技巧**：当 $\mathbf{q}_0 \cdot \mathbf{q}_1 < 0$ 时，先对 $\mathbf{q}_1$ 取反（等价旋转），避免绕过远路。

**多个方向插值**：使用样条（squad）而非链式 slerp。

**C++ 片段**（EasyMath 待实现）：
```cpp
Quaternion<T> slerp(const Quaternion<T>& q0, const Quaternion<T>& q1, T t) {
    T cosOmega = dot(q0, q1);
    if (cosOmega < 0) { cosOmega = -cosOmega; q1 = -q1; }  // 走近路
    T omega = std::acos(cosOmega);
    if (std::abs(omega) < epsilon) return q0;  // 防退化
    T sinOmega = std::sin(omega);
    return (std::sin((T(1) - t) * omega) / sinOmega) * q0
         + (std::sin(t * omega) / sinOmega) * q1;
}
```

#### 将一个向量旋转到另一个向量

给定 $\mathbf{u}, \mathbf{v}$（单位向量），构造一个将 $\mathbf{u}$ 转到 $\mathbf{v}$ 的四元数：

**注意**：$\mathbf{u} = \mathbf{v}$ 时返回单位四元数（无旋转）。
$\mathbf{u} = -\mathbf{v}$ 时需要选一个正交轴。

**中间变量**：
- $c = \mathbf{u} \cdot \mathbf{v}$（点积）
- $s = ||\mathbf{u} \times \mathbf{v}||$（叉积长度）

**结果**：
$$
\mathbf{q} = (s + c,\ v_x - u_x s/(1+c),\ ...\ )
$$

这种构造**避免了归一化**和三角函数——非常高效。

---

## 4.4 顶点混合（Skinning/Morphing）

**目的**：让网格顶点跟随多个骨骼（影响源）加权变换，常用于角色动画。

**数学**：
$$
\mathbf{v}' = \sum_{i=1}^{n} w_i(\mathbf{v}) \cdot M_i \mathbf{v}
$$
其中 $w_i$ 是第 $i$ 个骨骼对顶点 $\mathbf{v}$ 的权重（$\sum w_i = 1$），$M_i$ 是该骨骼的世界变换。

**双四元数蒙皮**（Dual Quaternion Skinning）：避免线性混合引起的"关节塌陷"。

---

## 4.5 变形（Morphing）

**核心问题**：从模型 A 平滑过渡到模型 B。

**前提**：必须有顶点对应（vertex correspondence）关系。

**简单方法**：插值对应顶点。

**变形目标**（morph target / blend shape）：预定义多个中间姿态，运行时加权混合（常用于面部表情）。

---

## 4.6 几何缓存回放

**思想**：将几何数据在 CPU 端预计算，存到 GPU 顶点缓冲区，运行时回放。
- 优势：每帧 CPU 计算量低
- 劣势：占用显存

---

## 4.7 投影（**EasyMath 已实现**）

### 4.7.1 正交投影（Orthographic）

**OpenGL 风格**（NDC $[-1, 1]$）：
$$
P_{ortho} = \begin{pmatrix} \frac{2}{r-l} & 0 & 0 & -\frac{r+l}{r-l} \\ 0 & \frac{2}{t-b} & 0 & -\frac{t+b}{t-b} \\ 0 & 0 & \frac{-2}{f-n} & -\frac{f+n}{f-n} \\ 0 & 0 & 0 & 1 \end{pmatrix}
$$

### 4.7.2 透视投影（Perspective）

$$
P_{persp} = \begin{pmatrix} \frac{1}{\text{aspect} \cdot \tan(\text{fov}/2)} & 0 & 0 & 0 \\ 0 & \frac{1}{\tan(\text{fov}/2)} & 0 & 0 \\ 0 & 0 & -\frac{f+n}{f-n} & -\frac{2fn}{f-n} \\ 0 & 0 & -1 & 0 \end{pmatrix}
$$

**C++ 片段**（EasyMath 已实现）：
```cpp
auto M = MTXPerspective<float>(fov, aspect, near, far);
auto O = MTXOrtho<float>(l, r, b, t, n, f);
```

---

## 与 EasyMath 项目关联

### 已实现 ✓
- 平移/旋转/缩放矩阵（`MTXTranslation` / `MTXRotationX/Y/Z` / `MTXScale`）
- 刚体变换
- 变换连接（`operator*` 链式）
- 矩阵分解（`submatrix` / `transpose` / `cofactor` / `adjugate` / `inverse`）
- 投影矩阵（`MTXPerspective` / `MTXOrtho` / `MTXLookAt`）

### 待实现（v1.1.0 路线图）
- **欧拉变换**（`fromEuler` / `toEuler`）
- **绕任意轴旋转**（`MTXRotationAxis`，用 Rodrigues）
- **四元数**：乘法、共轭、逆、`fromAxisAngle`、`toAxisAngle`、`slerp`、`nlerp`、`rotate` 向量
- **法线变换**（`Vector.normalTransform(M)`）
- **顶点混合**（蒙皮）

### EasyMath 已有但可参考 RTR4 改进的
- `MTXLookAt`：RTR4 的优化版本是"只归一化 2 次"（已实现）

## 关键术语索引

| 术语 | 释义 |
|------|------|
| 仿射变换 | 线性变换 + 平移 |
| 刚体变换 | 只含平移和旋转，保持长度/角度/手性 |
| 齐次坐标 | 4D 表示法，使平移可用矩阵乘法 |
| 万向锁 | 欧拉角中两轴重合，丢失自由度 |
| Slerp | 球面线性插值，四元数间平滑过渡 |
| 共轭四元数 | $(w, -\mathbf{v})$ |
| Rodrigues 公式 | 绕任意轴旋转矩阵的紧凑形式 |
| 法线变换 | $(M^{-1})^T \mathbf{n}$，不能直接用 $M$ |

## 重要参考

- 公式在 `RTR4-overview.md` → `formulas.md` 中有 LaTeX + C++ 摘要
- Shoemake 1985 是把四元数引入图形学的开创论文
- Goldman 1991 是绕任意轴旋转的另一种高效方法
