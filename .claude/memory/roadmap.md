# 开发路线图

## v1.1.0 开发计划

### 🔴 高优先级（核心功能）

#### 1. 完善 Quaternion API
**目标**: 使四元数成为完整可用的旋转表示

- [ ] **基础运算**
  - [ ] 四元数乘法（`operator*`）
  - [ ] 共轭（`conjugate()`）
  - [ ] 逆（`inverse()`）

- [ ] **创建与转换**
  - [ ] 从轴角创建（`fromAxisAngle(axis, angle)`）
  - [ ] 转为轴角（`toAxisAngle() -> (axis, angle)`）
  - [ ] 从旋转矩阵创建（`fromMatrix(Matrix3x3)`）
  - [ ] 转为旋转矩阵（`toMatrix() -> Matrix3x3/4x4`）

- [ ] **插值与旋转**
  - [ ] 球面线性插值（`slerp(q1, q2, t)`）
  - [ ] 归一化线性插值（`nlerp(q1, q2, t)` - 性能优化版本）
  - [ ] 旋转向量（`rotate(Vector3) -> Vector3`）

- [ ] **辅助功能**
  - [ ] 点积（`dot(q1, q2)`）
  - [ ] 角度差（`angleBetween(q1, q2)`）

#### 2. 欧拉角支持
**目标**: 提供欧拉角作为另一种旋转表示

- [ ] **创建 EulerAngles 类**
  ```cpp
  template<typename T>
  class EulerAngles {
      T pitch; // X 轴旋转
      T yaw;   // Y 轴旋转
      T roll;  // Z 轴旋转
      RotationOrder order; // XYZ, ZYX, etc.
  };
  ```

- [ ] **支持多种旋转顺序**
  - [ ] XYZ (Pitch-Yaw-Roll)
  - [ ] ZYX (Yaw-Pitch-Roll)
  - [ ] XZY, YXZ, YZX, ZXY

- [ ] **万向锁处理**
  - [ ] 检测万向锁情况
  - [ ] 提供警告或替代方案

#### 3. 旋转表示转换矩阵
**目标**: 实现所有旋转表示之间的无缝转换

| From ↓ / To → | Quaternion | EulerAngles | AxisAngle | Matrix3x3 |
|---------------|------------|-------------|-----------|-----------|
| **Quaternion**    | -          | [ ]         | [ ]       | [ ]       |
| **EulerAngles**   | [ ]        | -           | [ ]       | [ ]       |
| **AxisAngle**     | [ ]        | [ ]         | -         | [ ]       |
| **Matrix3x3**     | [ ]        | [ ]         | [ ]       | -         |

**实现函数**:
- [ ] `Quaternion::fromEuler(euler) -> Quaternion`
- [ ] `Quaternion::toEuler() -> EulerAngles`
- [ ] `EulerAngles::fromQuaternion(q) -> EulerAngles`
- [ ] `EulerAngles::toQuaternion() -> Quaternion`
- [ ] `AxisAngle::fromQuaternion(q) -> AxisAngle`
- [ ] `AxisAngle::toQuaternion() -> Quaternion`

---

### 🟡 中优先级（功能增强）

#### 4. Color 运算扩展
- [ ] **颜色混合模式**
  - [ ] Multiply
  - [ ] Screen
  - [ ] Overlay
  - [ ] Add
  - [ ] Subtract
  - [ ] Darken
  - [ ] Lighten

- [ ] **颜色插值**
  - [ ] 线性 RGB 插值
  - [ ] HSV 空间插值
  - [ ] HSL 空间插值
  - [ ] 感知均匀插值（LAB 空间）

- [ ] **颜色调整工具**
  - [ ] 色调偏移（`shiftHue(degrees)`）
  - [ ] 饱和度调整（`adjustSaturation(factor)`）
  - [ ] 亮度调整（`adjustBrightness(factor)`）
  - [ ] 对比度调整（`adjustContrast(factor)`）

#### 5. Matrix 投影支持
- [ ] 透视投影矩阵（`MTXPerspective(fov, aspect, near, far)`）
- [ ] 正交投影矩阵（`MTXOrtho(left, right, bottom, top, near, far)`）
- [ ] Look-At 矩阵（`MTXLookAt(eye, target, up)`）

---

### 🟢 低优先级（性能与质量）

#### 6. 性能优化
- [ ] **SIMD 支持**
  - [ ] SSE 优化（Vector3/4, Matrix4x4）
  - [ ] AVX 优化（批量运算）
  - [ ] NEON 优化（ARM 平台）

- [ ] **性能基准测试**
  - [ ] 建立基准测试框架
  - [ ] 对比其他数学库（GLM, Eigen）

- [ ] **内存布局优化**
  - [ ] 对齐优化（`alignas(16/32)`）
  - [ ] Cache-friendly 数据结构

#### 7. 测试与文档
- [ ] **单元测试**
  - [ ] Vector 单元测试
  - [ ] Matrix 单元测试
  - [ ] Quaternion 单元测试
  - [ ] Color 单元测试

- [ ] **API 文档补充**
  - [ ] Quaternion.md（中英文）
  - [ ] Color.md（中英文）
  - [ ] EulerAngles.md（中英文）

- [ ] **示例代码**
  - [ ] 旋转示例
  - [ ] 颜色处理示例
  - [ ] 矩阵变换示例

#### 8. 跨平台支持
- [ ] Linux 构建支持
- [ ] macOS 构建支持
- [ ] CMake 构建系统（作为 Premake 的替代）

---

## 未来版本规划

### v1.2.0 - 高级功能
- [ ] 样条曲线（Bezier, Catmull-Rom）
- [ ] 噪声函数（Perlin, Simplex）
- [ ] 包围盒（AABB, OBB）
- [ ] 碰撞检测基础

### v2.0.0 - 大版本重构
- [ ] 模块化设计（独立的 Vector/Matrix/Color 模块）
- [ ] 完整的 SIMD 抽象层
- [ ] GPU 计算支持（CUDA/Metal/Vulkan Compute）

---

## 开发优先级总结
1. **立即开始**: Quaternion 完善（v1.1.0 核心）
2. **紧随其后**: EulerAngles + 旋转转换
3. **功能增强**: Color 运算 + Matrix 投影
4. **质量提升**: 测试 + 文档 + 性能优化
