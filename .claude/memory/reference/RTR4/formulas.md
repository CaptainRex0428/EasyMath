# RTR4 关键公式汇总

> **约定**：每条公式标注 `来源章节 §小节`；`C++ 片段` 用 EasyMath 风格（`namespace EM`）展示实现思路（不复制源码，给出关键行）。

---

## A. 变换矩阵（来源：Ch4）

### A.1 平移矩阵（Ch4 §4.1.1）

$$
T(\mathbf{t}) = \begin{pmatrix} 1 & 0 & 0 & t_x \\ 0 & 1 & 0 & t_y \\ 0 & 0 & 1 & t_z \\ 0 & 0 & 0 & 1 \end{pmatrix}
$$

**C++ 片段**（EasyMath `MTXTranslation`）：
```cpp
auto T = MTXTranslation<float, 4>(tx, ty, tz);  // OpenGL 行主序
// 矩阵(0,3)=tx, (1,3)=ty, (2,3)=tz, (3,3)=1
```

### A.2 绕 X/Y/Z 轴旋转（Ch4 §4.1.2）

$$
R_x(\theta) = \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & c & -s & 0 \\ 0 & s & c & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix},\quad R_y(\theta) = \begin{pmatrix} c & 0 & s & 0 \\ 0 & 1 & 0 & 0 \\ -s & 0 & c & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix},\quad R_z(\theta) = \begin{pmatrix} c & -s & 0 & 0 \\ s & c & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}
$$

其中 $c=\cos\theta,\ s=\sin\theta$。

**C++ 片段**（EasyMath）：
```cpp
auto rx = MTXRotationX<float, 4>(theta);
auto ry = MTXRotationY<float, 4>(theta);
auto rz = MTXRotationZ<float, 4>(theta);
```

### A.3 缩放矩阵（Ch4 §4.1.3）

$$
S(\mathbf{s}) = \begin{pmatrix} s_x & 0 & 0 & 0 \\ 0 & s_y & 0 & 0 \\ 0 & 0 & s_z & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}
$$

### A.4 法线变换矩阵（Ch4 §4.1.7）

非均匀缩放或剪切变换下，法线不能用模型矩阵直接变换。正确做法是：

$$
\mathbf{n}' = (M^{-1})^T \mathbf{n}
$$

或者使用伴随矩阵共轭转置 $\mathbf{n}' = (M^{-T})\mathbf{n}$。

### A.5 绕任意轴旋转——Rodrigues 公式（Ch4 §4.2.4）

设单位轴 $\hat{\mathbf{k}} = (k_x, k_y, k_z)$，角度 $\theta$：

$$
R(\hat{\mathbf{k}}, \theta) = \begin{pmatrix} c+k_x^2(1-c) & k_x k_y (1-c) - k_z s & k_x k_z (1-c) + k_y s & 0 \\ k_y k_x (1-c) + k_z s & c + k_y^2(1-c) & k_y k_z(1-c) - k_x s & 0 \\ k_z k_x (1-c) - k_y s & k_z k_y(1-c) + k_x s & c + k_z^2(1-c) & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}
$$

Rodrigues 紧凑形式（无平移）：
$$
R = I + \sin\theta \cdot K + (1 - \cos\theta) \cdot K^2
$$
其中 $K$ 是 $\hat{\mathbf{k}}$ 的反对称矩阵（叉积矩阵）。

### A.6 投影矩阵（Ch4 §4.7）

**正交投影**（OpenGL NDC $[-1,1]$）：
$$
P_{ortho} = \begin{pmatrix} \frac{2}{r-l} & 0 & 0 & -\frac{r+l}{r-l} \\ 0 & \frac{2}{t-b} & 0 & -\frac{t+b}{t-b} \\ 0 & 0 & \frac{-2}{f-n} & -\frac{f+n}{f-n} \\ 0 & 0 & 0 & 1 \end{pmatrix}
$$

**透视投影**（OpenGL 风格，$w'=-z$）：
$$
P_{persp} = \begin{pmatrix} \frac{1}{\text{aspect} \cdot \tan(\text{fov}/2)} & 0 & 0 & 0 \\ 0 & \frac{1}{\tan(\text{fov}/2)} & 0 & 0 \\ 0 & 0 & -\frac{f+n}{f-n} & -\frac{2fn}{f-n} \\ 0 & 0 & -1 & 0 \end{pmatrix}
$$

**C++ 片段**（EasyMath）：
```cpp
auto proj_persp = MTXPerspective<float>(fov, aspect, near, far);
auto proj_ortho = MTXOrtho<float>(l, r, b, t, n, f);
```

### A.7 Look-At 视图矩阵（Ch2 §2.3 / Ch4 隐含）

```cpp
Vector3 forward = (target - eye).normalized();
Vector3 right   = cross(forward, up).normalized();
Vector3 camUp   = cross(right, forward);  // 已正交

return Matrix{
    right.x,  right.y,  right.z,  -dot(right, eye),
    camUp.x,  camUp.y,  camUp.z,  -dot(camUp, eye),
   -forward.x,-forward.y,-forward.z, dot(forward, eye),
    0,        0,        0,         1
};
```

---

## B. 四元数（来源：Ch4 §4.3）

### B.1 四元数定义

$$
\mathbf{q} = (q_w, q_x, q_y, q_z) = q_w + q_x \mathbf{i} + q_y \mathbf{j} + q_z \mathbf{k}
$$

3D 旋转四元数：$\mathbf{q} = (\cos(\theta/2), \sin(\theta/2) \cdot \hat{\mathbf{k}})$

### B.2 四元数乘法

$$
\mathbf{q}_1 \mathbf{q}_2 = \begin{pmatrix} w_1 w_2 - \mathbf{v}_1 \cdot \mathbf{v}_2 \\ w_1 \mathbf{v}_2 + w_2 \mathbf{v}_1 + \mathbf{v}_1 \times \mathbf{v}_2 \end{pmatrix}
$$

### B.3 旋转向量

$$
\mathbf{v}' = \mathbf{q} \mathbf{v} \mathbf{q}^{-1}
$$

（$\mathbf{v}$ 当作纯四元数 $(0, v_x, v_y, v_z)$ 处理）

**C++ 片段**（EasyMath `Quaternion`）：
```cpp
QuaternionF q{cosHalf, sinHalf*kx, sinHalf*ky, sinHalf*kz};
auto rotated = q.rotate(vec);  // TODO: 待 v1.1 实现
```

### B.4 Slerp 球面线性插值

$$
\text{slerp}(\mathbf{q}_0, \mathbf{q}_1; t) = \frac{\sin((1-t)\Omega)}{\sin\Omega} \mathbf{q}_0 + \frac{\sin(t\Omega)}{\sin\Omega} \mathbf{q}_1
$$

其中 $\cos\Omega = \mathbf{q}_0 \cdot \mathbf{q}_1$（点积）。

---

## C. 矩阵分解（来源：Ch4 §4.2.3）

任意可逆仿射矩阵可分解为平移、旋转、缩放、剪切（QR/SVD 分解等）。实用方法：

- **极分解（Polar Decomposition）**：$M = R \cdot S$，其中 $R$ 是旋转，$S$ 是对称正定缩放/剪切。
- **特征值分解**：用于提取主轴和缩放。

---

## D. 颜色与光（来源：Ch8）

### D.1 Rec.709 亮度系数

$$
Y = 0.2126 R + 0.7152 G + 0.0722 B
$$

> EasyMath `Color::luminance()` 实现的 Rec709 系数即此。

### D.2 sRGB ↔ Linear（Ch8 §8.1.4）

**sRGB → Linear**（解码）：
$$
C_{lin} = \begin{cases} C_{srgb} / 12.92, & C_{srgb} \le 0.04045 \\ \left(\frac{C_{srgb} + 0.055}{1.055}\right)^{2.4}, & C_{srgb} > 0.04045 \end{cases}
$$

**Linear → sRGB**（编码）：
$$
C_{srgb} = \begin{cases} 12.92 \cdot C_{lin}, & C_{lin} \le 0.0031308 \\ 1.055 \cdot C_{lin}^{1/2.4} - 0.055, & C_{lin} > 0.0031308 \end{cases}
$$

### D.3 WCAG 对比度

$$
\text{contrast}(C_1, C_2) = \frac{L_1 + 0.05}{L_2 + 0.05}
$$

其中 $L_1 \ge L_2$ 是两个颜色的 Rec.709 相对亮度（在线性空间计算）。

> EasyMath `Color::contrastRatio()` 即此。

---

## E. 着色模型（来源：Ch5, Ch9）

### E.1 Phong 模型（Ch5 §5.1）

$$
I = k_a I_a + \sum_i \left( k_d (\mathbf{n} \cdot \mathbf{l}_i) + k_s (\mathbf{r}_i \cdot \mathbf{v})^{\alpha} \right) I_i
$$

### E.2 Blinn-Phong 半向量模型

$$
I_{spec} = k_s (\mathbf{n} \cdot \mathbf{h})^{\alpha},\quad \mathbf{h} = \frac{\mathbf{l} + \mathbf{v}}{|\mathbf{l} + \mathbf{v}|}
$$

### E.3 Cook-Torrance 微表面 BRDF（Ch9 §9.8）

$$
f_r(\mathbf{l}, \mathbf{v}) = \frac{F(\mathbf{h}, \mathbf{l}) \cdot G(\mathbf{l}, \mathbf{v}, \mathbf{n}) \cdot D(\mathbf{h}, \alpha)}{4 (\mathbf{n} \cdot \mathbf{l})(\mathbf{n} \cdot \mathbf{v})}
$$

- $F$：菲涅尔项
- $G$：几何遮蔽
- $D$：法线分布函数（NDF）
- $\mathbf{h}$：半向量 $= (\mathbf{l}+\mathbf{v})/|\mathbf{l}+\mathbf{v}|$

### E.4 菲涅尔——Schlick 近似（Ch9 §9.5）

$$
F(\mathbf{h}, \mathbf{l}) \approx F_0 + (1 - F_0) (1 - \mathbf{h} \cdot \mathbf{l})^5
$$

其中 $F_0$ 是法向入射反射率：
- 常见电介质：$F_0 \approx 0.04$
- 金属：$F_0$ 即金属本身的反射颜色（金/银/铜各有 RGB 值）

### E.5 GGX 法线分布（Ch9 §9.8.1，各向同性）

$$
D(\mathbf{h}) = \frac{\alpha^2}{\pi ((\mathbf{n} \cdot \mathbf{h})^2 (\alpha^2 - 1) + 1)^2}
$$

其中 $\alpha = \text{roughness}^2$（UE/Disney 重映射）。

---

## F. 阴影（来源：Ch7）

### F.1 阴影贴图基础

1. 以光源为相机渲染深度到 shadow map
2. 主相机渲染时，对每个片元做 lookup，比较深度

### F.2 PCF（Percentage-Closer Filtering，Ch7 §7.5）

软阴影：对阴影贴图查询做 $N \times N$ 采样（如 $3\times3$），按"通过/失败"比例决定阴影强度。

### F.3 PCSS（Percentage-Closer Soft Shadows，Ch7 §7.6）

改进 PCF：根据遮挡物距离动态调整采样半径 → 接近遮挡物时硬、远离时软。

---

## G. 纹理（来源：Ch6）

### G.1 Mipmap LOD

$$
D = \log_2 \left( \max\left( \sqrt{\left(\frac{du}{dx}\right)^2 + \left(\frac{dv}{dx}\right)^2}, \sqrt{\left(\frac{du}{dy}\right)^2 + \left(\frac{dv}{dy}\right)^2} \right) \right)
$$

- 屏幕空间像素 $(x,y)$ 对应纹理坐标 $(u,v)$ 的偏导
- $D$ 是 mipmap 层级（$0$ = 原图，$D$ 越大图越小）

### G.2 各向异性过滤

当 $\frac{du}{dx}$ 和 $\frac{du}{dy}$ 差异大时（如观察斜面），用椭圆足迹采样 → RipMap、Earth mover's distance 等方案。

---

## H. 加速结构（来源：Ch19）

### H.1 包围球

球心 $c$，半径 $r$。点 $p$ 在球内当且仅当 $\|p - c\| \le r$。

### H.2 AABB（轴对齐包围盒）

点 $p$ 在 AABB $[b_{min}, b_{max}]$ 内当且仅当 $\forall i: b_{min,i} \le p_i \le b_{max,i}$。

### H.3 k-DOP（k-Discrete Oriented Polytope）

AABB 的推广：$k/2$ 对方向固定的法线。$k=6$ 即 AABB；$k=14, 18, 26$ 等用于更紧的包围。

### H.4 OBB（有向包围盒）

由中心 + 3 个正交轴方向 + 3 个半轴长度定义。点 $p$ 在 OBB 内当且仅当将 $p$ 变换到 OBB 局部坐标系后落在 $[-1,1]^3$ 内。

### H.5 BVH（Bounding Volume Hierarchy）

层次包围体：树形结构，每个内部节点是子节点的包围体（球/AABB/OBB）。

构造方法：自顶向下（SAH 启发式 / 中点切分）、自底向上。

---

## I. 相交测试（来源：Ch22）

### I.1 Möller-Trumbore 射线-三角形相交（Ch22 §22.8）

$$
\mathbf{t} = \begin{pmatrix} t \\ u \\ v \end{pmatrix} = \frac{1}{\mathbf{d} \cdot (\mathbf{e}_1 \times \mathbf{e}_2)} \begin{pmatrix} \mathbf{s} \cdot (\mathbf{e}_1 \times \mathbf{e}_2) \\ \mathbf{d} \cdot (\mathbf{s} \times \mathbf{e}_2) \\ \mathbf{s} \cdot (\mathbf{d} \times \mathbf{e}_1) \end{pmatrix}
$$

其中：
- $\mathbf{o}$：射线起点，$\mathbf{d}$：方向
- $\mathbf{v}_0, \mathbf{v}_1, \mathbf{v}_2$：三角形顶点
- $\mathbf{e}_1 = \mathbf{v}_1 - \mathbf{v}_0$，$\mathbf{e}_2 = \mathbf{v}_2 - \mathbf{v}_0$
- $\mathbf{s} = \mathbf{o} - \mathbf{v}_0$

命中条件：$t \ge 0$，$u, v \ge 0$，$u + v \le 1$。

### I.2 射线-AABB 相交（平板法，Ch22 §22.7.1）

对每个轴 $i$：
$$
t_{1,i} = \frac{b_{min,i} - o_i}{d_i},\quad t_{2,i} = \frac{b_{max,i} - o_i}{d_i}
$$
取 $t_{enter} = \max_i \min(t_{1,i}, t_{2,i})$，$t_{exit} = \min_i \max(t_{1,i}, t_{2,i})$。
命中当且仅当 $t_{enter} \le t_{exit}$ 且 $t_{exit} \ge 0$。

### I.3 射线-球 相交

$$
d^2 (\mathbf{o} + t\mathbf{d} - \mathbf{c}) = r^2
$$

展开为一元二次方程 $at^2 + bt + c = 0$，判别式 $\Delta = b^2 - 4ac$。

---

## J. 性能公式（来源：Ch18）

### J.1 SAH（Surface Area Heuristic）

$$
C_{split} = C_T + \frac{S_A}{S_{total}} N_A \cdot C_i + \frac{S_B}{S_{total}} N_B \cdot C_i
$$

用于 BVH 划分时评估哪个切分点最优（最小化期望遍历代价）。

---

## 提取完成度

- [x] A. 变换矩阵（来自 Ch4）
- [x] B. 四元数（来自 Ch4）
- [x] C. 矩阵分解（来自 Ch4）
- [x] D. 颜色与光（来自 Ch8）
- [x] E. 着色模型（来自 Ch5/Ch9）
- [x] F. 阴影（来自 Ch7）
- [x] G. 纹理（来自 Ch6）
- [x] H. 加速结构（来自 Ch19）
- [x] I. 相交测试（来自 Ch22）
- [x] J. 性能公式（来自 Ch18）
- [ ] K. 体积渲染（Ch14，待补）
- [ ] L. 碰撞检测（Ch25，待补）
- [ ] M. 光线追踪（Ch26，待补）
