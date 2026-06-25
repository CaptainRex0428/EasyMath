# Ch4 Transform 变换（**核心参考**）

> **与 EasyMath 关联**：直接对应 `include/Matrix.h`、`include/Quaternion.h`、`include/EasyConversion.h`。本章节公式的 90% 已由 EasyMath 实现。
> **完成度**：P0（核心参考）— 完整抽取，2026-06-25

**核心问题**：如何用数学方法操控 3D 物体的位置、朝向、尺寸和形状？

**章节定位**：实时渲染中最常用的工具集。本章内容是几乎所有后续章节的基础。

---

## 4.0 仿射变换基础

### 4.0.1 线性变换（Linear Transform）

**定义**：保持向量加法和标量乘法的变换
$$
L(\mathbf{a} + \mathbf{b}) = L(\mathbf{a}) + L(\mathbf{b}) \\
L(c\mathbf{a}) = c \cdot L(\mathbf{a})
$$

**例**：缩放、旋转都是线性变换。

**关键限制**：3D 线性变换必须用 $3\times 3$ 矩阵（不是 $2\times 2$ 或 $4\times 4$）。

### 4.0.2 仿射变换（Affine Transform）

**核心思想**：线性变换 + 平移的组合。

**齐次符号**：用 $4 \times 4$ 矩阵表示。
- 方向向量：$\mathbf{v} = (v_x, v_y, v_z, 0)$
- 点：$\mathbf{p} = (p_x, p_y, p_z, 1)$

**为什么 $4\times 4$？** 因为平移是非线性的（$\mathbf{p} + \mathbf{t}$），用 $3\times 3$ 矩阵无法表达；齐次坐标用 $w$ 分量来区分点和方向。

### 4.0.3 仿射变换的关键性质

- 保持直线（直线变换后仍是直线）
- **保持平行性**（不保证角度和长度）
- 不一定保持手性（镜像会反转）

### 4.0.4 实时渲染中的仿射变换

平移、旋转、缩放、对称、剪切都是仿射变换。

### 4.0.5 行主序 vs 列主序

| 顺序 | 内存布局 | 应用 |
|------|---------|------|
| **列主序**（column-major） | 平移分量在第 4 列前 3 行（OpenGL 风格） | OpenGL、RTR4 |
| **行主序**（row-major） | 平移分量在第 4 行前 3 列（DirectX 风格） | DirectX、HLSL |

**EasyMath 风格**：列主序（与 RTR4、OpenGL 一致）。

---

## 4.1 基本变换

### 4.1.1 平移（Translation）

$$
T(t_x, t_y, t_z) = \begin{pmatrix} 1 & 0 & 0 & t_x \\ 0 & 1 & 0 & t_y \\ 0 & 0 & 1 & t_z \\ 0 & 0 & 0 & 1 \end{pmatrix}
$$

**关键性质**：
- 点会被平移
- 方向向量（$w=0$）**不会**被平移
- 逆矩阵 = $T(-t_x, -t_y, -t_z)$

**C++ 片段**（EasyMath `MTXTranslation`）：
```cpp
auto M = MTXTranslation<float, 4>(tx, ty, tz);
// M(0,3) = tx, M(1,3) = ty, M(2,3) = tz
```

### 4.1.2 旋转（Rotation）

**2D 旋转推导**：设向量 $\mathbf{v} = (r\cos\phi, r\sin\phi)$，旋转 $\theta$ 后得 $(r\cos(\phi+\theta), r\sin(\phi+\theta))$。用二角和差公式展开：

$$
R(\theta) = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}
$$

**3D 绕坐标轴**：

$$
R_x(\theta) = \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & \cos\theta & -\sin\theta & 0 \\ 0 & \sin\theta & \cos\theta & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}
$$

$$
R_y(\theta) = \begin{pmatrix} \cos\theta & 0 & \sin\theta & 0 \\ 0 & 1 & 0 & 0 \\ -\sin\theta & 0 & \cos\theta & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}
$$

$$
R_z(\theta) = \begin{pmatrix} \cos\theta & -\sin\theta & 0 & 0 \\ \sin\theta & \cos\theta & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}
$$

**关键性质**：
- 行列式为 +1
- **正交矩阵**（$R^{-1} = R^T$）
- 多个旋转矩阵的乘积仍是正交矩阵
- 旋转矩阵 $R(\theta)$ 与 $R(-\theta)$ 互为逆
- **旋转矩阵的迹**（与角度相关）：

$$
\text{trace}(R) = 1 + 2\cos\theta
$$

可用此从矩阵反解旋转角度（但要注意 $\cos\theta$ 的多解）。

**绕某点旋转**（不是绕轴）：
$$M = T(\mathbf{q}) \cdot R(\theta) \cdot T(-\mathbf{q})$$

**组合顺序**：右乘先应用（post-multiplication）。例如 `T(2,0,0) * Rz(90°) * v` 先旋转再平移。

### 4.1.3 缩放（Scaling）

$$
S(s_x, s_y, s_z) = \begin{pmatrix} s_x & 0 & 0 & 0 \\ 0 & s_y & 0 & 0 \\ 0 & 0 & s_z & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}
$$

**术语**：
- **均匀缩放**（uniform）：$s_x = s_y = s_z$（保持形状）
- **非均匀缩放**（nonuniform）：不同缩放因子（变形）
- **各向同性**（isotropic）= 均匀缩放
- **各向异性**（anisotropic）= 非均匀缩放

**逆矩阵**：$S^{-1} = S(1/s_x, 1/s_y, 1/s_z)$。

**缩放 $w$ 分量技巧**：将 $M(3,3)$ 设为 $1/k$ 也可以实现 $k$ 倍均匀缩放（通过齐次除法），但效率低，因为引入额外除法。

**反射（Reflection / Mirror）**：
- 1 个或 3 个分量为负 → 镜像
- 2 个分量为负 → 旋转 180°（不是镜像）

**检测反射**：计算左上 $3 \times 3$ 行列式。
- 正 → 正常缩放
- 负 → 包含反射

**为什么重要**：反射会反转三角形顶点顺序（CW ↔ CCW），影响背面剔除 → 需检测并特殊处理。

**按任意方向缩放**（不在坐标轴上）：
1. 找到目标方向的标准正交基 $\{\mathbf{u}, \mathbf{v}, \mathbf{w}\}$
2. 构造基底变换矩阵 $B = [\mathbf{u} | \mathbf{v} | \mathbf{w}]$
3. 缩放：$M = B \cdot S \cdot B^T$（$B$ 正交，$B^{-1} = B^T$）

**C++ 片段**（EasyMath 待实现）：
```cpp
Matrix<T, 4, 4> MTXRotationAxis(const Vector<T, 3>& axis, T angle) {
    // 沿任意方向缩放
    Vector<T, 3> u = axis.normalized();
    Vector<T, 3> v = findPerpendicular(u);  // 找任意正交向量
    Vector<T, 3> w = cross(u, v);
    Matrix<T, 3, 3> B{u, v, w};
    return B * S(scale) * B.transpose();
}
```

### 4.1.4 剪切（Shear）

**6 个基本剪切矩阵**（$H_{ij}$：$j$ 坐标剪切 $i$ 坐标）：

$$
H_{xy} = \begin{pmatrix} 1 & 0 & 0 & 0 \\ a & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix},\quad H_{xz} = \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ b & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix},\quad H_{yx} = \begin{pmatrix} 1 & c & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}
$$

（$H_{yz}, H_{zx}, H_{zy}$ 类似）

**关键性质**：任何剪切矩阵的行列式 = 1 → **体积保持**（volume-preserving）。

**应用**：游戏效果（迷幻扭曲、招牌倾斜）、木纹等。

### 4.1.5 变换的连接（Concatenation）

**关键事实**：矩阵乘法**不满足交换律**。

$$
M_1 M_2 \neq M_2 M_1
$$

**约定**：右乘 = 先应用（post-multiplication），左乘 = 后应用（pre-multiplication）。

**例**（3 个矩阵组合）：
$$
\mathbf{v}_{transformed} = T \cdot R \cdot S \cdot \mathbf{v}
$$
- $S$ 先作用（缩放）
- 然后 $R$（旋转）
- 最后 $T$（平移）

**为什么组合**：效率提升。把 3 个矩阵组合成 1 个 → 每顶点只做 1 次矩阵乘法（而不是 3 次）。

**结合律可用**（即使不满足交换律）：
$$
ABC = A(BC) = (AB)C
$$
可分组预计算（如先算刚体子矩阵 $TR$ 再组合）。

### 4.1.6 刚体变换（Rigid-body Transform）

**定义**：只含平移和旋转。**保持长度、角度、手性**。

$$
M_{\text{rigid}} = \begin{pmatrix} r_{00} & r_{01} & r_{02} & t_x \\ r_{10} & r_{11} & r_{12} & t_y \\ r_{20} & r_{21} & r_{22} & t_z \\ 0 & 0 & 0 & 1 \end{pmatrix}
$$

**逆矩阵**（利用 $R^{-1} = R^T$）：
$$
M_{\text{rigid}}^{-1} = \begin{pmatrix} r_{00} & r_{10} & r_{20} & -t_x \\ r_{01} & r_{11} & r_{21} & -t_y \\ r_{02} & r_{12} & r_{22} & -t_z \\ 0 & 0 & 0 & 1 \end{pmatrix}
$$

**应用**：相机视角、刚体物理（球、箱子）。

### 4.1.7 法线变换（**重要**）

**核心问题**：为什么不能直接用模型矩阵变换法线？

**原因**：非均匀缩放会改变法线和切线之间的夹角。例：$S = \text{diag}(2, 1, 1)$ 拉伸 $x$，但法线如果用 $S$ 变换就不垂直于拉伸后的切线了。

**正确做法**：用模型矩阵的**逆转置**（adjoint transpose）：

$$
\mathbf{n}' = (M^{-1})^T \mathbf{n}
$$

**为什么**：
- 切线变换：$\mathbf{t}' = M \mathbf{t}$
- 法线应满足：$\mathbf{t} \cdot \mathbf{n} = 0$（变换前正交）
- 变换后也要：$\mathbf{t}' \cdot \mathbf{n}' = 0$
- 即 $(M\mathbf{t})^T \mathbf{n}' = 0$
- 即 $\mathbf{t}^T M^T \mathbf{n}' = 0$
- 因为 $\mathbf{t}^T \mathbf{n} = 0$，所以需要 $M^T \mathbf{n}' = \mathbf{n}$
- 即 $\mathbf{n}' = (M^T)^{-1} \mathbf{n} = (M^{-1})^T \mathbf{n}$ ✓

**替代方案**：使用**伴随矩阵**（$\text{adj}(M) = \det(M) \cdot M^{-1}$）：
$$
\mathbf{n}' = \text{adj}(M)^T \mathbf{n}
$$
（避免显式求逆，更高效）

**实用建议**：
- 纯旋转/平移时 $M^{-T} = M$，可直接用模型矩阵
- 变换后记得**归一化**（长度可能改变）
- 从变换后的三角形边向量叉乘重算法线 → 不需要这个修正

**C++ 片段**（EasyMath 待实现）：
```cpp
Vector<T, 3> Vector<T, 3>::transformNormal(const Matrix<T, 4, 4>& M) const {
    // n' = (M^-1)^T * n，假设 n 是方向向量
    Matrix<T, 4, 4> inv = M.inverse();
    Vector<T, 4> r = inv.transpose() * Vector<T, 4>(data[0], data[1], data[2], T(0));
    return Vector<T, 3>{r[0], r[1], r[2]}.normalized();
}
```

### 4.1.8 计算逆矩阵

**三种方法**（按可用信息分类）：

1. **已知类型**（平移/旋转/缩放/刚体）→ 闭式求逆（见 4.1.1-4.1.6）
2. **通用 4×4** → 伴随矩阵法 $\det(M) \cdot M^{-1} = \text{adj}(M)$
3. **数值方法**（LU/Gauss-Jordan）→ O(n³)

**性能提示**：
- 如果逆矩阵**只用于变换向量** → 只算左上 3×3 部分的转置即可
- 刚体矩阵的逆很便宜（$R^T$ 计算 + 简单平移取反）
- 通用伴随矩阵法贵（9 个 3×3 余子式 + 行列式）

---

## 4.2 特殊矩阵变换

### 4.2.1 欧拉变换（Euler Transform）

**核心思想**：用三个角度（heading/pitch/roll 或 yaw/pitch/roll）描述一个方向。

**24 种组合**：
- 6 个两轴旋转（XY、YZ、ZX、XZ、YX、ZY）× 2（内旋/外旋）= 12
- 6 个三轴旋转 × 2 = 12
- 总 24 种

**典型约定**：heading 绕 $Y$，pitch 绕 $X$，roll 绕 $Z$（head→pitch→roll 顺序）

**矩阵形式**：
$$
M_{\text{euler}} = R_z(\text{roll}) \cdot R_x(\text{pitch}) \cdot R_y(\text{heading})
$$
（最后应用的旋转在**最左边**——本例是 roll）

**直觉**：
- 改变 heading → 摇头
- 改变 pitch → 点头
- 改变 roll → 歪头

### 4.2.2 从欧拉变换中提取参数

**步骤**：
1. 直接读出平移（最后一列前 3 行）
2. 从左上 $3\times 3$ 反解 3 个角度

**以 head/pitch/roll 为例**：

$$
M = R_z(r) \cdot R_x(p) \cdot R_y(h) = \begin{pmatrix} m_{00} & m_{01} & m_{02} \\ m_{10} & m_{11} & m_{12} \\ m_{20} & m_{21} & m_{22} \end{pmatrix}
$$

**展开后**：
- $m_{21} = \sin p$（如果 $p \in [-90°, 90°]$，可直接 $\arcsin$）
- $m_{20} = -\cos p \sin h$
- $m_{22} = \cos p \cos h$

**反解**：
- $p = \arcsin(m_{21})$
- $h = \arctan2(-m_{20}, m_{22})$
- $r = \arctan2(-m_{12}, m_{11})$
- （假设 $m_{11} = \cos p \cos r$，$m_{12} = -\cos p \sin r$）

**特殊情况**：$\cos p = 0$（pitch = ±90°）→ **万向锁**。

**万向锁（Gimbal Lock）**：当中间角度为 $\pm 90°$ 时，相邻两轴重合，丢失一个自由度。**导致无数解**。

**避免方法**：
- 用四元数
- 用轴角表示
- 接受万向锁（保持解的连续性时可用）

### 4.2.3 矩阵分解

**目标**：从总变换矩阵 $M$ 提取子变换（平移、旋转、缩放、剪切）。

**步骤**：
1. **平移**：直接读 $M$ 的第 4 列前 3 行
2. **旋转**：用 Gram-Schmidt 正交化 $M$ 的左上 $3\times 3$
3. **缩放**：$M$ 的列向量的长度
4. **剪切**：从正交化前后差异得到

**应用**：
- 撤销已有变换
- 提取相机位置/朝向
- 重用组件

**算法参考**：Thomas 1997、Goldman 2005、Shoemake 1994。

### 4.2.4 绕任意轴旋转——Rodrigues 公式

**输入**：单位轴 $\hat{\mathbf{k}} = (k_x, k_y, k_z)$，角度 $\theta$。

**推导思路**：
1. 找到任意将 $\hat{\mathbf{k}}$ 转到 $z$ 轴的旋转 $A$
2. 绕 $z$ 轴旋转 $\theta$：$R_z(\theta)$
3. 反向旋转 $A^{-1} = A^T$（正交矩阵）
4. 组合：$R_{\hat{\mathbf{k}}}(\theta) = A^T \cdot R_z(\theta) \cdot A$

**Rodrigues 紧凑形式**（**关键公式**）：
$$
R = I + \sin\theta \cdot K + (1 - \cos\theta) \cdot K^2
$$
其中 $K$ 是 $\hat{\mathbf{k}}$ 的**反对称矩阵**（叉积矩阵）：
$$
K = \begin{pmatrix} 0 & -k_z & k_y \\ k_z & 0 & -k_x \\ -k_y & k_x & 0 \end{pmatrix}
$$

**展开形式**（4×4 齐次）：
$$
R(\hat{\mathbf{k}}, \theta) = \begin{pmatrix}
c + k_x^2(1-c) & k_x k_y(1-c) - k_z s & k_x k_z(1-c) + k_y s & 0 \\
k_y k_x(1-c) + k_z s & c + k_y^2(1-c) & k_y k_z(1-c) - k_x s & 0 \\
k_z k_x(1-c) - k_y s & k_z k_y(1-c) + k_x s & c + k_z^2(1-c) & 0 \\
0 & 0 & 0 & 1
\end{pmatrix}
$$
其中 $c = \cos\theta, s = \sin\theta$。

**关键性质**：
- $K$ 反对称：$K^T = -K$
- $K^2$ 对称
- $K^3 = -K$（幂循环）
- 因此 $K$ 的泰勒展开有限项

**C++ 片段**（EasyMath 待实现）：
```cpp
Matrix<T, 4, 4> MTXRotationAxis(const Vector<T, 3>& axis, T angle) {
    Vector<T, 3> k = axis.normalized();
    T c = std::cos(angle);
    T s = std::sin(angle);
    T omc = T(1) - c;

    Matrix<T, 4, 4> R;
    R(0,0) = c + k[0]*k[0]*omc;
    R(0,1) = k[0]*k[1]*omc - k[2]*s;
    R(0,2) = k[0]*k[2]*omc + k[1]*s;
    R(1,0) = k[1]*k[0]*omc + k[2]*s;
    R(1,1) = c + k[1]*k[1]*omc;
    R(1,2) = k[1]*k[2]*omc - k[0]*s;
    R(2,0) = k[2]*k[0]*omc - k[1]*s;
    R(2,1) = k[2]*k[1]*omc + k[0]*s;
    R(2,2) = c + k[2]*k[2]*omc;
    R(3,3) = T(1);
    return R;
}
```

**Goldman 替代方案**（更高效）：用两向量叉积法构造正交基，再合成旋转。

---

## 4.3 四元数（**核心**）

### 4.3.1 数学背景

**定义**（Hamilton 1843）：
$$
\mathbf{q} = q_w + q_x \mathbf{i} + q_y \mathbf{j} + q_z \mathbf{k}
$$
其中 $\mathbf{i}^2 = \mathbf{j}^2 = \mathbf{k}^2 = \mathbf{ijk} = -1$。

**复数扩展**：复数 $\mathbb{C}$ = 实部 + 虚部；四元数 $\mathbb{H}$ = 1 实部 + 3 虚部（$i, j, k$）。

**Robinson 1958**：四元数已用于刚体模拟（早于 Shoemake 引入图形学）。

#### 基本运算

**加法**（逐分量）：
$$
\mathbf{q}_0 + \mathbf{q}_1 = (w_0 + w_1,\ x_0 + x_1,\ y_0 + y_1,\ z_0 + z_1)
$$

**标量乘法**：
$$
s\mathbf{q} = (sw, sx, sy, sz)
$$

**点积**（类似 4D 向量）：
$$
\mathbf{q}_0 \cdot \mathbf{q}_1 = w_0 w_1 + x_0 x_1 + y_0 y_1 + z_0 z_1
$$

**乘法**（用点积+叉积）：
$$
\mathbf{q}_0 \mathbf{q}_1 = (w_0 w_1 - \mathbf{v}_0 \cdot \mathbf{v}_1,\ w_0 \mathbf{v}_1 + w_1 \mathbf{v}_0 + \mathbf{v}_0 \times \mathbf{v}_1)
$$
其中 $\mathbf{q} = (w, \mathbf{v})$，$\mathbf{v} = (x, y, z)$。

**共轭**（conjugate）：
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

**单位四元数**（$||\mathbf{q}|| = 1$）：
- 满足 $\mathbf{q}^{-1} = \mathbf{q}^*$（实部不变，虚部取反即可）
- 满足 $\mathbf{q}^* = \mathbf{q}^{-1}$ 与 $\mathbf{q}^{-1} = \mathbf{q}^*$ 等价

**指数与对数**：

**指数**：给定纯四元数 $\mathbf{q} = (0, \phi \mathbf{u}_q)$（$\mathbf{u}_q$ 单位）：
$$
e^{\mathbf{q}} = (\cos\phi, \sin\phi \mathbf{u}_q)
$$
（恰好是单位四元数！）

**对数**（单位四元数的逆）：
$$
\log(\mathbf{q}) = (0, \phi \mathbf{u}_q)
$$
其中 $\phi = \arccos(w)$，$\mathbf{u}_q = \mathbf{v}/||\mathbf{v}||$。

**幂**（用于 slerp）：
$$
\mathbf{q}^t = e^{t \log \mathbf{q}} = (\cos t\phi, \sin t\phi \mathbf{u}_q)
$$

#### 单位四元数与旋转

**关键事实**：单位四元数 $e^{\mathbf{q}} = (\cos\phi, \sin\phi \mathbf{u}_q)$ 表示绕单位轴 $\mathbf{u}_q$ 旋转 $2\phi$。

**构造**（从轴角）：
$$
\mathbf{q} = (\cos(\theta/2), \sin(\theta/2) \cdot \hat{\mathbf{k}})
$$

**等价表示**：$\mathbf{q}$ 和 $-\mathbf{q}$ 表示**相同旋转**（$\theta$ 与 $\theta + 2\pi$ 相同）。

**3b1b 互动视频**（<https://eater.net/quaternions>）：四元数可视化最佳资源。

### 4.3.2 四元数变换

#### 旋转向量（核心公式）

$$
\mathbf{v}' = \mathbf{q} \mathbf{v} \mathbf{q}^{-1}
$$
其中 $\mathbf{v}$ 当作纯四元数 $(0, v_x, v_y, v_z)$。

**计算展开**（$\mathbf{q} = (w, \mathbf{u})$，$\mathbf{v}$ 是 3D 向量）：
$$
\mathbf{v}' = (w^2 - \mathbf{u} \cdot \mathbf{u})\mathbf{v} + 2(\mathbf{u} \cdot \mathbf{v})\mathbf{u} + 2w(\mathbf{u} \times \mathbf{v})
$$

**C++ 片段**（EasyMath `Quaternion`，v1.1 待实现）：
```cpp
template<typename T>
Vector<T, 3> Quaternion<T>::rotate(const Vector<T, 3>& v) const {
    Vector<T, 3> u(data[1], data[2], data[3]);  // 虚部
    T w = data[0];
    return (w*w - dot(u, u)) * v
         + T(2) * dot(u, v) * u
         + T(2) * w * cross(u, v);
}
```

#### 矩阵转换

**四元数 → 旋转矩阵**（**无三角函数**，仅 2 次乘法和减法）：
$$
M = \begin{pmatrix}
1-2y^2-2z^2 & 2xy-2wz & 2xz+2wy & 0 \\
2xy+2wz & 1-2x^2-2z^2 & 2yz-2wx & 0 \\
2xz-2wy & 2yz+2wx & 1-2x^2-2y^2 & 0 \\
0 & 0 & 0 & 1
\end{pmatrix}
$$

**反向**（矩阵 → 四元数）：用迹区分情况，**数值稳定性**很关键。

**方法**：先选最大对角元 $M_{ii}$（最稳定的），由此推算其他分量。

设 $\mathbf{q} = (w, x, y, z)$，$4w^2 = 1 + m_{00} + m_{11} + m_{22}$。

- 若 $4w^2$ 最大：$w = \sqrt{4w^2}/2$，$x = (m_{21} - m_{12})/(4w)$，$y = (m_{02} - m_{20})/(4w)$，$z = (m_{10} - m_{01})/(4w)$
- 若 $4x^2$ 最大：$x = \sqrt{4x^2}/2$，依此类推

**核心原则**：避免除以过小数值。

#### 球面线性插值（Slerp）

**核心问题**：在两个旋转之间平滑过渡。

**数学公式**：
$$
\text{slerp}(\mathbf{q}_0, \mathbf{q}_1; t) = \frac{\sin((1-t)\Omega)}{\sin\Omega} \mathbf{q}_0 + \frac{\sin(t\Omega)}{\sin\Omega} \mathbf{q}_1
$$
其中 $\cos\Omega = \mathbf{q}_0 \cdot \mathbf{q}_1$（点积）。

**实现技巧**：
1. 当 $\mathbf{q}_0 \cdot \mathbf{q}_1 < 0$ 时，先对 $\mathbf{q}_1$ 取反（绕"近路"）
2. 当 $\Omega \approx 0$（两四元数几乎相同）→ 退化为 lerp（避免除零）
3. 当 $\Omega$ 接近 $\pi$（两四元数相对）→ 仍能工作

**nlerp**（归一化线性插值）：$n\mathbf{q} = (1-t)\mathbf{q}_0 + t\mathbf{q}_1$，然后归一化。**更快但角度速度不恒定**。

**多四元数插值**：用 **squad**（球面二次插值）而非链式 slerp。

**C++ 片段**（EasyMath 待实现）：
```cpp
template<typename T>
Quaternion<T> Quaternion<T>::slerp(const Quaternion<T>& q0,
                                   const Quaternion<T>& q1, T t) {
    T cosOmega = q0.dot(q1);
    Quaternion<T> q1Adj = q1;
    if (cosOmega < T(0)) {
        cosOmega = -cosOmega;
        q1Adj = -q1;  // 走近路
    }
    T omega = std::acos(std::min(T(1), std::max(T(-1), cosOmega)));
    if (std::abs(omega) < T(1e-6)) {
        // 退化：退化为 lerp
        return Quaternion<T>((T(1) - t) * q0 + t * q1Adj).normalized();
    }
    T sinOmega = std::sin(omega);
    return (std::sin((T(1) - t) * omega) / sinOmega) * q0
         + (std::sin(t * omega) / sinOmega) * q1Adj;
}
```

#### 将一个向量旋转到另一个向量

**核心问题**：构造将 $\mathbf{u}$ 转到 $\mathbf{v}$ 的四元数（最短路径）。

**特殊情况**：
- $\mathbf{u} = \mathbf{v}$ → 单位四元数（无旋转）
- $\mathbf{u} = -\mathbf{v}$ → 需选一个正交轴

**构造公式**：
1. 归一化 $\mathbf{u}$ 和 $\mathbf{v}$
2. 中间变量：
   - $c = \mathbf{u} \cdot \mathbf{v}$（点积）
   - 轴方向 $\mathbf{axis} = \mathbf{u} \times \mathbf{v}$（叉积）
3. 单位四元数：

$$
\mathbf{q} = \frac{(1 + c,\ \mathbf{axis}_x,\ \mathbf{axis}_y,\ \mathbf{axis}_z)}{\sqrt{(1+c)^2 + \|\mathbf{axis}\|^2}}
$$

**注意**：$\|\mathbf{axis}\|^2 = 1 - c^2$（叉积长度平方等于 $1 - c^2$），所以分母 = $\sqrt{2(1+c)} = 2\cos(\Omega/2)$。

**对反向情况**（$\mathbf{u} = -\mathbf{v}$）：选任意与 $\mathbf{u}$ 正交的轴。

**优势**：**无 sqrt 和三角函数**（分母可化简为 $\sqrt{2(1+c)}$）→ 比欧拉角法（用 $\arccos$ 和叉积归一化）更快。

---

## 4.4 顶点混合（Skinning）

**目的**：让网格顶点跟随多个骨骼（影响源）加权变换。

**数学**：
$$
\mathbf{v}' = \sum_{i=1}^{n} w_i(\mathbf{v}) \cdot M_i \mathbf{v}
$$
其中 $w_i$ 是第 $i$ 个骨骼对顶点 $\mathbf{v}$ 的权重（$\sum w_i = 1$），$M_i$ 是该骨骼的世界变换。

**问题**：权重超过范围（$\sum w_i > 1$）→ 体积膨胀。

**双四元数蒙皮**（Dual Quaternion Skinning）：避免线性混合引起的"关节塌陷"。

---

## 4.5 变形（Morphing）

**核心问题**：从模型 A 平滑过渡到模型 B。

**前提**：必须有顶点对应（vertex correspondence）关系。

**简单方法**：插值对应顶点（顶点必须一一对应）。

**变形目标 / Blend Shape**：
- 预定义多个中间姿态
- 运行时加权混合
- 常用于面部表情

**变体**：Maya 的 Blend Shape、Blender 的 Shape Keys。

---

## 4.6 几何缓存回放

**思想**：将几何数据在 CPU 端预计算，存到 GPU 顶点缓冲区，运行时回放。

**优势**：每帧 CPU 计算量低

**劣势**：占用显存

**应用**：骨骼动画、几何变形

---

## 4.7 投影（**EasyMath 已实现**）

### 4.7.1 正交投影（Orthographic）

**OpenGL 风格**（NDC $[-1, 1]$）：
$$
P_{ortho} = \begin{pmatrix} \frac{2}{r-l} & 0 & 0 & -\frac{r+l}{r-l} \\ 0 & \frac{2}{t-b} & 0 & -\frac{t+b}{t-b} \\ 0 & 0 & \frac{-2}{f-n} & -\frac{f+n}{f-n} \\ 0 & 0 & 0 & 1 \end{pmatrix}
$$

**变体**：DirectX 风格用 NDC $[0, 1]$（用于 z 缓冲）—— 等价但矩阵第 3 行不同。

### 4.7.2 透视投影（Perspective）

$$
P_{persp} = \begin{pmatrix} \frac{1}{\text{aspect} \cdot \tan(\text{fov}/2)} & 0 & 0 & 0 \\ 0 & \frac{1}{\tan(\text{fov}/2)} & 0 & 0 \\ 0 & 0 & -\frac{f+n}{f-n} & -\frac{2fn}{f-n} \\ 0 & 0 & -1 & 0 \end{pmatrix}
$$

**特点**：透视除法后 $w = -z$（OpenGL）或 $w = z$（DirectX）。

**C++ 片段**（EasyMath 已实现）：
```cpp
auto M = MTXPerspective<float>(fov, aspect, near, far);
auto O = MTXOrtho<float>(l, r, b, t, n, f);
```

### 4.7.3 投影重要参数

- **FOV**（Field of View）：垂直视野角度
- **aspect**：宽高比（width/height）
- **near**：近裁剪面（> 0）
- **far**：远裁剪面（> near）

**Z-fighting 风险**：near/far 范围越大 → 深度精度越差。

**Z 范围**：
- OpenGL：$z \in [-1, 1]$（NDC）
- DirectX：$z \in [0, 1]$（NDC）
- 透视除法后都是 [0, 1]（最终深度缓冲值）

---

## 4.8 表 4.1：变换矩阵速查

| 变换 | 矩阵 | 性质 | 逆 |
|------|------|------|------|
| 平移 | $T$ | 保长度（仅位移） | $T(-\mathbf{t})$ |
| 旋转 | $R$ | 正交、$\det = 1$ | $R^T$ |
| 缩放 | $S$ | 体积变化 $\prod s_i$ | $S(1/s_i)$ |
| 反射 | 负缩放因子 | $\det = -1$（行序反） | 自己 |
| 剪切 | $H$ | 体积不变 ($\det = 1$) | $H^{-1}$ 显式 |
| 刚体 | $T \cdot R$ | 正交前 3×3 | $R^T \cdot T(-)$ |

**正交矩阵的核心优势**：$R^{-1} = R^T$ → 求逆极快。

---

## 与 EasyMath 项目关联

### 已实现 ✓
- 平移/旋转/缩放矩阵（`MTXTranslation` / `MTXRotationX/Y/Z` / `MTXScale`）
- 刚体变换
- 变换连接（`operator*` 链式）
- 矩阵分解（`submatrix` / `transpose` / `cofactor` / `adjugate` / `inverse`）
- 投影矩阵（`MTXPerspective` / `MTXOrtho` / `MTXLookAt`）

### 待实现（v1.1.0 路线图）
- **欧拉变换**（`fromEuler` / `toEuler`）—— 注意万向锁
- **绕任意轴旋转**（`MTXRotationAxis`，用 Rodrigues 公式）
- **按任意方向缩放**（用基底变换）
- **四元数**：乘法、共轭、逆、`fromAxisAngle`、`toAxisAngle`、`slerp`、`nlerp`、`rotate` 向量
- **法线变换**（`Vector.normalTransform(M)`，用 $(M^{-1})^T$）
- **顶点混合 / 蒙皮**
- **变形目标（Blend Shape）**

### EasyMath 已有但可参考 RTR4 改进的
- `MTXLookAt`：RTR4 的优化版本是"只归一化 2 次"（已实现）
- 法线变换：当前没有 → v1.1 应加

## 关键术语索引

| 术语 | 释义 |
|------|------|
| 仿射变换 | 线性变换 + 平移 |
| 刚体变换 | 只含平移和旋转，保持长度/角度/手性 |
| 齐次坐标 | 4D 表示法，使平移可用矩阵乘法 |
| 万向锁 | 欧拉角中两轴重合，丢失自由度 |
| Slerp | 球面线性插值，四元数间平滑过渡 |
| nlerp | 归一化线性插值，快但不恒角速 |
| 共轭四元数 | $(w, -\mathbf{v})$ |
| Rodrigues 公式 | 绕任意轴旋转矩阵的紧凑形式 |
| 法线变换 | $(M^{-1})^T \mathbf{n}$，不能直接用 $M$ |
| 反对称矩阵 | $K^T = -K$，用于叉积 $K\mathbf{v} = \mathbf{a} \times \mathbf{v}$ |
| 旋转矩阵的迹 | $1 + 2\cos\theta$ |
| 行列式检测反射 | $\det < 0$ → 含反射 |
| 体积保持 | $\det = 1$（剪切、旋转） |
| 列主序 vs 行主序 | OpenGL vs DirectX |

## 重要参考

- 公式在 `RTR4-overview.md` → `formulas.md` 中有 LaTeX + C++ 摘要
- Shoemake 1985 是把四元数引入图形学的开创论文
- Goldman 1991 是绕任意轴旋转的另一种高效方法
- Möller & Hughes 1999：高效四元数到矩阵转换

## 抽取历史

- 2026-06-25：基于 `RTR4/_raw/ch004.txt`（102KB 干净文本）抽取
- 来源章节号：Ch4
- 跳过内容：无（完整抽取）
- 置信度：高（与 EasyMath API 直接对应）
