# Matrix Demo 设计文档

**最后更新**: 2026-06-22
**对应 C++ 文件**: `include/Matrix.h`

## Demo 概述

Matrix Demo 展示 3D 图形学中最核心的矩阵变换管道，帮助学习者理解：
1. **模型变换** - 平移、旋转、缩放
2. **相机变换** - LookAt 视图矩阵
3. **投影变换** - 透视投影 vs 正交投影
4. **矩阵组合** - MVP 矩阵的正确乘法顺序

---

## 管道架构

### 完整变换管道

```
[Local Space]
     │
     │ M (Model Matrix)
     │ MTXTranslation * MTXRotation * MTXScale
     ↓
[World Space]
     │
     │ V (View Matrix)
     │ MTXLookAt(eye, target, up)
     ↓
[Camera Space]
     │
     │ P (Projection Matrix)
     │ MTXPerspective(fov, aspect, near, far)
     │ 或 MTXOrtho(left, right, bottom, top, near, far)
     ↓
[Clip Space]
     │
     │ Perspective Division (自动)
     ↓
[NDC Space]
     │
     │ Viewport Transform (自动)
     ↓
[Screen Space]
```

### 矩阵组合公式

```
MVP = P × V × M

FinalVertex = MVP × LocalVertex
```

**重要**: 矩阵乘法顺序是**从右到左**！

```javascript
const MVP = P.multiply(V).multiply(M);
// 或者
const MVP = P * (V * M);
```

---

## 实现细节

### 1. Model Matrix (模型矩阵)

**用途**: 将顶点从局部空间变换到世界空间

**组成**: 缩放 × 旋转 × 平移

```javascript
const S = MTXScale(sx, sy, sz);
const Rx = MTXRotationX(rx * DEG_TO_RAD);
const Ry = MTXRotationY(ry * DEG_TO_RAD);
const Rz = MTXRotationZ(rz * DEG_TO_RAD);
const T = MTXTranslation(tx, ty, tz);

// 组合顺序：T * R * S
const M = T.multiply(Rz).multiply(Ry).multiply(Rx).multiply(S);
```

**为什么是 TRS 而不是 SRT？**
- 先缩放（S）
- 然后旋转（R）
- 最后平移（T）
- 这样物体先在原点变换，再移动到目标位置

**对应 C++ 函数**:
```cpp
// Matrix.h
template<typename T>
Matrix<T, 4, 4> MTXTranslation(T x, T y, T z);

template<typename T>
Matrix<T, 4, 4> MTXRotationX(T radians);
Matrix<T, 4, 4> MTXRotationY(T radians);
Matrix<T, 4, 4> MTXRotationZ(T radians);

template<typename T>
Matrix<T, 4, 4> MTXScale(T x, T y, T z);
```

---

### 2. View Matrix (视图矩阵)

**用途**: 将顶点从世界空间变换到相机空间

**实现**: LookAt 函数

```javascript
const eye = [ex, ey, ez];
const target = [ltx, lty, ltz];
const up = [0, 1, 0]; // 固定向上方向

const V = MTXLookAt(eye, target, up);
```

**LookAt 原理**:
1. 计算相机方向向量 `forward = normalize(target - eye)`
2. 计算右向量 `right = normalize(cross(forward, up))`
3. 计算相机上向量 `camUp = cross(right, forward)`
4. 构建旋转矩阵将世界对齐到相机坐标系

**对应 C++ 函数**:
```cpp
// Matrix.h
template<typename T>
Matrix<T, 4, 4> MTXLookAt(
    const Vector<T, 3>& eye,
    const Vector<T, 3>& target,
    const Vector<T, 3>& up
);
```

---

### 3. Projection Matrix (投影矩阵)

**用途**: 将顶点从相机空间变换到裁剪空间

#### 透视投影 (Perspective)

**特点**: 近大远小，符合人眼视觉

```javascript
const aspect = canvas.width / canvas.height;
const fov = 60 * DEG_TO_RAD;
const near = 0.1;
const far = 100;

const P = MTXPerspective(fov, aspect, near, far);
```

**透视投影原理**:
1. 将视锥体 (Frustum) 映射到标准立方体
2. 深度值非线性分布（更多精度给近处物体）
3. 透视除法后 x, y ∈ [-1, 1], z ∈ [-1, 1]

**对应 C++ 函数**:
```cpp
// Matrix.h
template<typename T>
Matrix<T, 4, 4> MTXPerspective(T fov, T aspect, T near, T far);
```

#### 正交投影 (Orthographic)

**特点**: 平行投影，不产生透视效果

```javascript
const size = 4;
const aspect = canvas.width / canvas.height;
const halfW = size * aspect;
const halfH = size;

const P = MTXOrtho(-halfW, halfW, -halfH, halfH, 0.1, 100);
```

**正交投影用途**:
- 2D 游戏
- CAD 制图
- UI 渲染
- 等轴测视角游戏

**对应 C++ 函数**:
```cpp
// Matrix.h
template<typename T>
Matrix<T, 4, 4> MTXOrtho(
    T left, T right,
    T bottom, T top,
    T near, T far
);
```

---

## WebGL 实现

### 矩阵存储

**JavaScript (行主序)**:
```javascript
class Matrix {
    constructor(data) {
        // data 是长度为 16 的数组
        // 存储顺序：行主序
        // [0, 1, 2, 3]      = 第一行
        // [4, 5, 6, 7]      = 第二行
        // [8, 9, 10, 11]    = 第三行
        // [12, 13, 14, 15]  = 第四行
        this.data = data;
    }
}
```

**WebGL Shader (列主序)**:
```glsl
uniform mat4 uMVP;

// shader 中矩阵是列主序存储
// 需要转置！
```

**传递给 Shader**:
```javascript
// WebGL uniformMatrix4fv 参数：
// location, transpose, value
gl.uniformMatrix4fv(
    uMVPLocation,
    false,        // 不转置（我们手动转置）
    transpose(MVP) // 转置后的数据
);
```

### 转置函数

```javascript
function transpose(m) {
    const out = new Array(16);
    for (let i = 0; i < 4; i++) {
        for (let j = 0; j < 4; j++) {
            out[j * 4 + i] = m[i * 4 + j];
        }
    }
    return out;
}
```

---

## UI 设计

### 控制面板

```
┌─────────────────────────────┐
│ 🔢 Pure WebGL · Matrix      │
├─────────────────────────────┤
│ Model — 平移 Translation     │
│ tx [━━━○━━━] 0.0           │
│ ty [━━━○━━━] 0.0           │
│ tz [━━━○━━━] 0.0           │
├─────────────────────────────┤
│ Model — 旋转 Rotation (°)   │
│ rx [━━━○━━━] 0°            │
│ ry [━━━○━━━] 0°            │
│ rz [━━━○━━━] 0°            │
├─────────────────────────────┤
│ Model — 缩放 Scale          │
│ sx [━━━○━━━] 1.0           │
│ sy [━━━○━━━] 1.0           │
│ sz [━━━○━━━] 1.0           │
└─────────────────────────────┘
```

### 矩阵显示

```
┌─ Model Matrix ─────────────┐
│ MTXTranslation * MTX...    │
│ ┌────────────────────────┐ │
│ │ 0.00  0.00  0.00  0.00│ │
│ │ 0.00  0.00  0.00  0.00│ │
│ │ 0.00  0.00  0.00  0.00│ │
│ │ 0.00  0.00  0.00  1.00│ │
│ └────────────────────────┘ │
└────────────────────────────┘
```

### 管道可视化

```
[Local] → [World] → [Camera] → [Clip] → [NDC]
   M         Model      View        Proj
```

---

## 渲染对象

### 立方体

**顶点数据**:
```javascript
// 8 个顶点
const vertices = new Float32Array([
  // 前面
  -0.5, -0.5,  0.5,   // 0
   0.5, -0.5,  0.5,   // 1
   0.5,  0.5,  0.5,   // 2
  -0.5,  0.5,  0.5,   // 3
  // 后面
  -0.5, -0.5, -0.5,   // 4
   0.5, -0.5, -0.5,   // 5
   0.5,  0.5, -0.5,   // 6
  -0.5,  0.5, -0.5,   // 7
]);
```

**索引**:
```javascript
const indices = new Uint16Array([
  // 前面
  0, 1, 2,  0, 2, 3,
  // 后面
  4, 5, 6,  4, 6, 7,
  // ... 其他面
]);
```

**颜色**:
```javascript
const colors = new Float32Array([
  // 每个顶点一个 RGB 颜色
  1.0, 0.0, 0.0,  // 顶点 0: 红
  0.0, 1.0, 0.0,  // 顶点 1: 绿
  // ...
]);
```

### 地面网格

**用途**: 参考坐标系，帮助理解空间关系

**实现**:
```javascript
// 使用 VP 而不是 MVP！
// 这样网格固定在世界空间
const VP = P.multiply(V);
gl.uniformMatrix4fv(uMVP, false, transpose(VP));

// 绘制网格线
for (let i = -10; i <= 10; i++) {
    drawLine([-10, 0, i], [10, 0, i]);  // 平行于 X 轴
    drawLine([i, 0, -10], [i, 0, 10]);  // 平行于 Z 轴
}
```

---

## Shader 设计

### Vertex Shader

```glsl
attribute vec3 aPosition;
attribute vec3 aColor;

uniform mat4 uMVP;

varying vec3 vColor;

void main() {
    // 变换顶点位置
    gl_Position = uMVP * vec4(aPosition, 1.0);

    // 传递颜色到片段着色器
    vColor = aColor;
}
```

### Fragment Shader

```glsl
precision mediump float;

varying vec3 vColor;

void main() {
    gl_FragColor = vec4(vColor, 1.0);
}
```

---

## 交互设计

### 滑块控制

每个滑块都有：
- **标签** - 参数名（tx, rx, fov 等）
- **滑块** - 可拖动的范围控制
- **数值显示** - 当前值

### 投影切换

```javascript
function setProj(type) {
  projType = type;
  updateButtons();
  recalculateMatrices();
}

function updateButtons() {
  document.getElementById('btn-persp').className =
    projType === 'persp' ? 'toggle-btn active' : 'toggle-btn';
  document.getElementById('btn-ortho').className =
    projType === 'ortho' ? 'toggle-btn active' : 'toggle-btn';
}
```

### 实时更新

```javascript
function onSliderChange() {
  // 更新数值显示
  updateValueDisplays();

  // 重新计算矩阵
  recalculateMatrices();

  // 更新矩阵显示
  updateMatrixDisplay();
}

// 所有滑块绑定此事件
document.getElementById('tx').addEventListener('input', onSliderChange);
```

---

## 常见问题

### Q: 为什么立方体看起来变形了？

A: 检查：
1. 矩阵乘法顺序是否正确（P × V × M）
2. 透视投影的 aspect 是否正确
3. 是否正确传递了 MVP 给 shader

### Q: 为什么网格跟着立方体转？

A: 网格应该使用 VP 而不是 MVP：
```javascript
// 错误
drawGrid(MVP);

// 正确
drawGrid(VP);
```

### Q: 透视投影看起来不对？

A: 检查：
1. FOV 是否合理（推荐 45-75 度）
2. near/far 比例是否过大（影响深度精度）
3. aspect 是否正确（width / height）

### Q: 正交投影和透视看起来一样？

A: 检查：
1. 是否真的切换了投影类型
2. 正交投影的视口大小是否合适
3. 相机距离是否足够远（透视效果不明显）

---

## 性能优化

### 1. 避免每帧重新分配内存

```javascript
// ❌ 错误：每帧分配新数组
function render() {
  const MVP = new Matrix(16);  // 每帧分配！
  // ...
}

// ✅ 正确：复用数组
const MVP = new Matrix(16);
function render() {
  MVP.update();  // 复用
  // ...
}
```

### 2. 缓存矩阵结果

```javascript
// 如果某些参数不变，缓存结果
let cachedP = null;

function getProjection() {
  if (!cachedP || projChanged) {
    cachedP = MTXPerspective(fov, aspect, near, far);
  }
  return cachedP;
}
```

### 3. 减少矩阵乘法

```javascript
// ❌ 错误：重复计算
const MVP = P.multiply(V).multiply(M);
gl.uniformMatrix4fv(uMVP, false, transpose(MVP));
gl.uniformMatrix4fv(uMVP2, false, transpose(MVP));

// ✅ 正确：复用结果
const MVP = P.multiply(V).multiply(M);
const transposed = transpose(MVP);
gl.uniformMatrix4fv(uMVP, false, transposed);
gl.uniformMatrix4fv(uMVP2, false, transposed);
```

---

## 扩展方向

### 未来可以添加

1. **更多变换**
   - 剪切变换 (Shear)
   - 反射变换 (Reflection)
   - 自定义轴旋转

2. **更多投影**
   - 斜投影 (Oblique)
   - 立体投影 (Stereoscopic)

3. **动画演示**
   - 关键帧动画
   - 路径动画
   - 弹簧动画

4. **高级可视化**
   - 矩阵分解显示
   - 变换管道动画
   - 顶点轨迹显示

---

## 相关资源

### C++ 源码

- `include/Matrix.h` - 矩阵实现
- `include/Vector.h` - 向量运算
- `sandbox/sandbox.cpp` - C++ 测试

### Demo 文件

- `demo/src/pages/matrix/index.astro` - 统一 Matrix Demo（WebGL + Three.js，2026-06-22 重构）
- `demo/src/pages/matrix/webgl.astro` - 重定向到 `/matrix/?tab=webgl`
- `demo/src/pages/matrix/threejs.astro` - 重定向到 `/matrix/?tab=threejs`

### 文档

- `docs/Matrix_API.md` - Matrix API 文档
- `docs/Matrix_CN.md` - Matrix 中文文档

---

**维护者**: CaptainRex
**最后更新**: 2026-06-22
