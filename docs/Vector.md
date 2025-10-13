向量 API 说明
============

概览
----
`EM::Vector<T, N>` 是维度泛化、定长的向量模板，`T` 必须为算术类型。支持元素访问、swizzle、算术运算、范数与归一化、投影与反射、插值、齐次坐标转换以及多种图形学实用函数；在 Debug 下通过 `assert` 做边界与前置条件检查。

模板参数
--------
- `T`：元素类型（float/double/int 等）
- `N`：维度（N > 0）

构造
----
- 默认零初始化：`Vector()`
- 标量填充：`explicit Vector(const T& value)`
- 列表构造：`Vector(std::initializer_list<T>)`（长度必须等于 `N`）


元素访问
--------
- `operator[](size_t)` 读写 / 常量访问
- 边界检查：`at(size_t)` 读写 / 常量访问
- 迭代器：`begin()/end()`（含常量版本）
- 原始数据：`T* Data()` / `const T* Data() const`
- 尺寸与维度：`static constexpr size_t size()`、`getDimension()`、`Dimension()`

Swizzle
-------
- 一维别名：`VectorFilterDimension1D`（x,y,z,w / r,g,b,a，含大小写）
- 二维 swizzle：`operator[](VectorFilterDimension2D)`，返回 `Vector<T,2>`（xy, xz, yz, xw, yw, zw 等）
- 三维 swizzle：`operator[](VectorFilterDimension3D)`，返回 `Vector<T,3>`（xyz, xyw, yzw 等）

矩阵转换
--------
- 列向量：`Matrix<T, N, 1> toColMatrix()`
- 行向量：`Matrix<T, 1, N> toRowMatrix()`

算术运算（逐元素，除非特别说明）
------------------------------
- 原地：`+=`、`-=`、`*=`、`/=`（向量或标量）
- 友元：`+`、一元 `-`、`-`、`*`、`/`（向量或标量）
- 比较：`==`、`!=`（浮点类型使用 `NEARZERO_THRESHOLD` 作为 epsilon）

范数与归一化
------------
- `length(bool dimensionalityReduction = false) const`
- `constexpr lengthSquared(bool dimensionalityReduction = false) const`
- `normalized(bool dimensionalityReduction = false) const`
- `normalize(bool dimensionalityReduction = false)`
- `isZero(T epsilon = T{NEARZERO_THRESHOLD}) const`
- `isNormalized(T epsilon = T{NEARZERO_THRESHOLD}, bool dimensionalityReduction = true) const`
  - 当 `dimensionalityReduction = true` 时，范数仅计算前三个分量；`normalized()` 会保持其余维度不变

向量间运算
--------------
- 点积：`dot(a, b)`
- 叉积（3D）：`cross(a, b)`，返回 `Vector<T,3>`
- 2D 叉积（标量）：`cross2D(a, b)`
- 距离：`distance(a, b, bool dimensionalityReduction = true)`
- 距离平方：`distanceSquared(a, b, bool dimensionalityReduction = true)`
- 夹角：`angle(a, b, bool dimensionalityReduction = true)`（弧度）

插值
----
- 线性插值：成员 `lerp(other, t)` 与自由函数 `lerp(a, b, t)`
- 球面插值（单位向量）：`slerp(a, b, t)`

投影与反射
----------
- `project(const Vector& onto)` / `project(const Vector& onto, bool dimensionalityReduction)`
- `reflect(const Vector& normal)`（要求单位法向，N ≥ 2）

齐次坐标
--------
- 转齐次：`template<size_t newDim = N+1> toHomogeneous(T w = T{1}) const`
- 还原齐次：`template<size_t newDim = N-1> fromHomogeneous() const`（`w==0` 视为方向向量，不做透视除法）

图形学辅助
----------
- 平移矩阵（N==3）：`toTranslationMatrix(bool usedWithOrient=false)` -> `Matrix<T,4,4>`
- 反对称矩阵：`skewSymmetric()`（N=2/3 有专门构造，N>3 提供通用构造）
- 二维专用反对称矩阵：`skewSymmetric_2D()` -> `Matrix<T,3,3>`（仅 N==2）

输出
----
- `print(std::ostream&) const`，形如 `(x, y, z, ...) (V mode: xyzw)`

常用类型别名
------------
- `Vector2/3/4`（float）
- `Vector2f/3f/4f`、`Vector2d/3d/4d`、`Vector2i/3i/4i`

示例
----
```cpp
using namespace EM;

Vector3 a{1,2,3};
Vector3 b{2,0,1};

float d = dot(a,b);
auto c = cross(a,b);
auto an = a.normalized();

auto mat = a.toColMatrix();
auto angleRad = angle(a,b);
auto r = a.reflect(Vector3{0,1,0});

auto h = a.toHomogeneous();
auto eu = h.fromHomogeneous();
```

注意事项
--------
- `dimensionalityReduction` 仅对前三个分量进行范数相关计算，`normalized()` 会保留其他维度
- 反射需要单位法向量；不确定时先对法向做归一化
- 浮点比较使用固定 epsilon（`NEARZERO_THRESHOLD`）

扩展计划
--------
- 更多 swizzle 组合、clamp/saturate、随机与单位向量生成
- SIMD 专用化与小尺寸向量优化


