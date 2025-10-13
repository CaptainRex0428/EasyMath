矩阵 API 说明
============

概览
----
`EM::Matrix<T, rows, cols>` 为行主序、定长的矩阵模板，`T` 必须为算术类型（`std::is_arithmetic_v<T>`）。在 Debug 下通过 `assert` 做越界与前置条件检查。

模板参数
--------
- `T`：元素类型（float/double/int 等）
- `rows`：行数（>0）
- `cols`：列数（>0）

构造
----
- 默认构造：零初始化存储
- 列表构造：`Matrix(std::initializer_list<T>)`（行主序，长度需等于 `rows*cols`）

元素访问
--------
- 线性下标（行主序）：`operator[]`
- 二维访问：`operator()(row, col)`
- 维度：`Rows()` / `Cols()`

格式化输出
----------
- `print(std::ostream&)`；`operator<<` 转发到 `print`，默认定宽并保留 3 位小数

子矩阵
------
- `submatrix(size_t r, size_t c)`：去除第 `r` 行和第 `c` 列后的 `(rows-1)×(cols-1)` 矩阵

行列式 / 代数余子式 / 伴随矩阵 / 逆
------------------------------
- `determinant()`（仅方阵）：1×1 直接返回元素；2×2 使用闭式；更高阶按首行拉普拉斯展开
- `cofactor(r, c)`：代数余子式
- `cofactorMatrix()`：余子式矩阵
- `adjugate()`：伴随矩阵（余子式矩阵的转置）
- `inverse()`：需 `det != 0`；实现为 `adjugate() / det`

代数运算
--------
- 原地：`+=`、`-=`、`*=`（标量）
- 友元：标量乘、矩阵加减、一元负号、矩阵×矩阵

矩阵与向量乘法
--------------
- `Matrix<T, rows, cols> * Vector<T, cols> -> Vector<T, rows>`
- `Vector<T, rows> * Matrix<T, rows, cols> -> Vector<T, cols>`

转置
----
- `transpose()`：返回 `cols×rows`

全局辅助函数
------------
- 单位矩阵：`MTXIdentity<T, N>()`
- 3D 旋转（弧度）：`MTXRotationX/Y/Z<T, N>(T radians)`
- 3D 平移：`MTXTranslation<T, N>(x, y, z, bool usedWithOrient=false)`（若 `usedWithOrient=true`，将 (3,3) 置 0 以便与包含齐次项的朝向矩阵组合）
- 3D 缩放：`MTXScale<T, N>(sx, sy, sz)`

存储与约定
----------
- 行主序：线性索引 `i = r*cols + c`
- 初始化列表为行主序
- 旋转角度为弧度
- Debug 断言：越界、维度匹配、方阵要求、求逆非奇异

复杂度说明
----------
- `determinant()` 的拉普拉斯展开在最坏情况下为阶乘复杂度，适合小规模（≤4）；更大规模可考虑 LU 分解（未来计划）
- 矩阵乘法复杂度为 O(rows*cols*inner)

示例
----
```cpp
using namespace EM;

auto rx = MTXRotationX<float, 4>(3.14159f * 0.5f);
auto t  = MTXTranslation<float, 4>(1.0f, 2.0f, 3.0f);
auto m  = rx * t;

auto inv = m.inverse();
auto tr  = m.transpose();

Vector<float,4> v{1,2,3,1};
auto v_out = m * v;
```

注意事项与扩展计划
------------------
- `inverse()` 在行列式为 0 时断言失败
- 浮点数数值稳定性取决于条件数；若来自噪声数据，可对旋转列/行做规范化
- `print()` 默认输出 3 位小数
- 计划：基于 LU 的行列式/逆、更多构造器、分块/视图工具


