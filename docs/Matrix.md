# `Matrix`类

[![Static Badge](https://img.shields.io/badge/Back_to_README-gray)](../README.md "点击返回README") | [![Static Badge](https://img.shields.io/badge/Language-zh-red)](Matrix.en.md "click to english version")

---

![Static Badge](https://img.shields.io/badge/Type-class-green)
![Static Badge](https://img.shields.io/badge/State-completed-blue)


Matrix是一个长度可变、类型可变的泛型矩阵

| Attribute     | Value                                                              |
| ------------- | ------------------------------------------------------------------ |
| **Namespace** | **`EasyMath`**                                                     |
| **File**      | [**`include\Matrix.h`**](../include/Matrix.h)                      |
| **Parent**    | -                                                                  |
| **Feature**   | ![Static Badge](https://img.shields.io/badge/Feature-template-red) |


## 声明
```C++
template <typename T, size_t rows, size_t cols, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
class Matrix
```
### 模板
| Type     | Name   | Description                                            |
| -------- | ------ | ------------------------------------------------------ |
| typename | `T`    | 元素类型，受到 `std::is_arithmetic_v` （算数类型）限制 |
| size_t   | `rows` | 矩阵行数                                               |
| size_t   | `cols` | 矩阵列数                                               |
### 构建函数
| Access   | Modifier                                                                       | Return | Name                                            | Parameter                                                                                 | Description                                            |
| -------- | ------------------------------------------------------------------------------ | ------ | ----------------------------------------------- | ----------------------------------------------------------------------------------------- | ------------------------------------------------------ |
| `public` | -                                                                              |        | `Matrix()`                                      | -                                                                                         | 默认构建，全部元素初始化为0                            |
| `public` | -                                                                              |        | `Matrix(std::initializer_list<T>)`              | `std::initializer_list<T> InitializeList`-类型为`T`的初始化数列                           | 如果数列初始化标量元素数量不等于 `rows * cols`则会报错 |
| `public` | -                                                                              |        | `Matrix(std::initializer_list<Vector<T,cols>>)` | `std::initializer_list<Vector<T,cols>> InitializeList`-类型为`Vector<T,cols>`的初始化数列 | 如果数列初始化向量元素数量不等于 `rows`则会报错        |
| `public` | `template<typename T2, typename = std::enable_if_t<std::is_arithmetic_v<T2>>>` |        | `Matrix(const Matrix<T2,rows,cols>&)`           | `const Matrix<T2,rows,cols>& other`-同维度另一类型的矩阵                                  |                                                        |
```C++
EasyMath::Matrix<float,3,3> MatrixA;
EasyMath::Matrix<float,4,4> MatrixB{ 1, 2, 3, 4,
                               5, 6, 7, 8,
                               9,10,11,12,
                              13,14,15,16};

EasyMath::Vector<float,4> VectorA{1,2,3,4};
EasyMath::Vector<float,4> VectorB{5,6,7,8};
EasyMath::Vector<float,4> VectorC{9,10,11,12};
EasyMath::Vector<float,4> VectorD{13,14,15,16};

EasyMath::Matrix<float,4,4> MatrixC{VectorA,VectorB,VectorC,VectorD};

EasyMath::Matrix<float,3,4> MatrixD_T{VectorA,VectorB,VectorC};
EasyMath::Matrix<float,4,3> MatrixD = MatrixD_T.transpose();
```
## 数据访问
| Access   | Modifier                        | Return     | Name                        | Parameter                                                   | Description                                  |
| -------- | ------------------------------- | ---------- | --------------------------- | ----------------------------------------------------------- | -------------------------------------------- |
| `public` | `virtual`                       | `T&`       | `operator[](size_t)`        | `size_t idx`-元素索引位置                                   | -                                            |
| `public` | `virtual` `const`               | `const T&` | `operator[](size_t)`        | `size_t idx`-元素索引位置                                   | -                                            |
| `public` | `virtual`                       | `T&`       | `at(size_t)`                | `size_t idx`-元素索引位置                                   | -                                            |
| `public` | `virtual` `const`               | `const T&` | `at(size_t)`                | `size_t idx`-元素索引位置                                   | -                                            |
| `public` | `virtual`                       | `T&`       | `operator()(size_t,size_t)` | `size_t row`-行元素索引位置,<br>`size_t col`-列元素索引位置 | -                                            |
| `public` | `virtual` `const`               | `const T&` | `operator()(size_t,size_t)` | `size_t row`-行元素索引位置,<br>`size_t col`-列元素索引位置 | -                                            |
| `public` | `static` `constexpr`            | `size_t`   | `Rows()`                    | -                                                           | 矩阵行数                                     |
| `public` | `static` `constexpr`            | `size_t`   | `Cols()`                    | -                                                           | 矩阵列数                                     |
| `public` | `virtual` `noexcept`            | `T*`       | `Data()`                    | -                                                           | 直接返回data数列的指针                       |
| `public` | `virtual` `const` `noexcept`    | `const T*` | `Data()`                    | -                                                           | 直接返回data数列的指针                       |
| `public` | `virtual` `noexcept`            | `T*`       | `begin()`                   | -                                                           | 返回data数列第一个元素的指针                 |
| `public` | `virtual` `const` `noexcept`    | `const T*` | `begin()`                   | -                                                           | 返回data数列第一个元素的指针                 |
| `public` | `virtual` `noexcept`            | `T*`       | `end()`                     | -                                                           | 返回data数列的末尾指针（最后一个元素索引+1） |
| `public` | `virtual` `const` `noexcept`    | `const T*` | `end()`                     | -                                                           | 返回data数列的末尾指针（最后一个元素索引+1） |
| `public` | `static` `constexpr` `noexcept` | `size_t`   | `Size()`                    | -                                                           | 返回向量的维度（数组的维度）                 |
## 运算支持
### 复合赋值运算
| Access   | Modifier  | Return                   | Name                                       | Parameter                                                    | Description |
| -------- | --------- | ------------------------ | ------------------------------------------ | ------------------------------------------------------------ | ----------- |
| `public` | `virtual` | `Matrix<T, rows, cols>&` | `operator+=(const Matrix<T, rows, cols>&)` | `const Matrix<T, rows, cols>& matrix`-同维度矩阵             | -           |
| `public` | `virtual` | `Matrix<T, rows, cols>&` | `operator+=(const T&)`                     | `const T& scalar`-常量                                       | -           |
| `public` | `virtual` | `Matrix<T, rows, cols>&` | `operator-=(const Matrix<T, rows, cols>&)` | `const Matrix<T, rows, cols>& matrix`-同维度矩阵             | -           |
| `public` | `virtual` | `Matrix<T, rows, cols>&` | `operator-=(const T&)`                     | `const T& scalar`-常量                                       | -           |
| `public` | `virtual` | `Matrix<T, rows, cols>&` | `operator*=(const Matrix<T, cols, cols>&)` | `const Matrix<T, cols, cols>& matrix`-矩阵（累乘必须为方阵） | -           |
| `public` | `virtual` | `Matrix<T, rows, cols>&` | `operator*=(const T&)`                     | `const T& scalar`-常量                                       | -           |
| `public` | `virtual` | `Matrix<T, rows, cols>&` | `operator/=(const T&)`                     | `const T& scalar`-常量                                       | -           |
### 二元运算
| Access             | Modifier                             | Return                   | Name                                                                     | Parameter                                                                              | Description                                                                       |
| ------------------ | ------------------------------------ | ------------------------ | ------------------------------------------------------------------------ | -------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------- |
| `private` `friend` | -                                    | `Matrix<T, rows, cols>`  | `operator*(const Matrix<T, rows, cols>&, const T&)`                      | `const Matrix<T, rows, cols>& matrix`-矩阵,<br>`const T& scalar`-标量                  | -                                                                                 |
| `private` `friend` | -                                    | `Matrix<T, rows, cols>`  | `operator*(const T&, const Matrix<T, rows, cols>&)`                      | `const T& scalar`-标量,<br>`const Matrix<T, rows, cols>& matrix`-矩阵                  | -                                                                                 |
| `private` `friend` | `template<typename T, size_t cols2>` | `Matrix<T, rows, cols2>` | `operator*(const Matrix<T, rows, cols>&, const Matrix<T, cols, cols2>&)` | `const Matrix<T, rows, cols>& lhs`-矩阵A,<br>`const Matrix<T, cols, cols2>& rhs`-矩阵B | -                                                                                 |
| `private` `friend` | -                                    | `Matrix<T, rows, 1>`     | `operator*(const Matrix<T, rows, cols>&, const Vector<T, cols>&)`        | `const Matrix<T, rows, cols>& matrix`-矩阵,<br>`const Vector<T, cols>& vec`-向量       | 矩阵 × 向量（列向量）: Matrix<T, rows, cols> × Vector<T, cols> -> Vector<T, rows> |
| `private` `friend` | -                                    | `Matrix<T, 1, cols>`     | `operator*( const Vector<T, cols>&, const Matrix<T, rows, cols>&)`       | `const Vector<T, cols>& vec`-向量,<br>`const Matrix<T, rows, cols>& matrix`-矩阵       | 向量（行向量）× 矩阵: Vector<T, rows> × Matrix<T, rows, cols> -> Vector<T, cols>  |
| `private` `friend` | -                                    | `Matrix<T, rows, cols>`  | `operator+(const Matrix<T, rows, cols>&, const Matrix<T, rows, cols>&)`  | `const Matrix<T, rows, cols>& lhs`-矩阵A,<br>`const Matrix<T, rows, cols>& rhs`-矩阵B  | -                                                                                 |
| `private` `friend` | -                                    | `Matrix<T, rows, cols>`  | `operator+(const Matrix<T, rows, cols>&, const T&)`                      | `const Matrix<T, rows, cols>& lhs`-矩阵,<br>`const T& scalar`-标量                     | -                                                                                 |
| `private` `friend` | -                                    | `Matrix<T, rows, cols>`  | `operator-(const Matrix<T, rows, cols>&, const Matrix<T, rows, cols>&)`  | `const Matrix<T, rows, cols>& lhs`-矩阵A,<br>`const Matrix<T, rows, cols>& rhs`-矩阵B  | -                                                                                 |
| `private` `friend` | -                                    | `Matrix<T, rows, cols>`  | `operator－(const Matrix<T, rows, cols>&, const T&)`                     | `const Matrix<T, rows, cols>& lhs`-矩阵,<br>`const T& scalar`-标量                     | -                                                                                 |
```C++
EasyMath::Vector<float,4> VectorA{1,2,3,4};
EasyMath::Vector<float,4> VectorB{5,6,7,8};
EasyMath::Vector<float,4> VectorC{9,10,11,12};
EasyMath::Vector<float,4> VectorD{13,14,15,16};
EasyMath::Matrix<float,4,4> MatrixC{VectorA,VectorB,VectorC,VectorD};

std::cout << VectorA * MatrixC << std::endl;
std::cout << MatrixC * VectorA  << std::endl;
```
输出
```
Matrix 1x4
-                                -
| 90.000 100.000 110.000 120.000 |
-                                -

Matrix 4x1
-        -
| 30.000 |
| 70.000 |
|110.000 |
|150.000 |
-        -
```

### 取反
| Access             | Modifier | Return                  | Name                                      | Parameter                                  | Description |
| ------------------ | -------- | ----------------------- | ----------------------------------------- | ------------------------------------------ | ----------- |
| `private` `friend` | -        | `Matrix<T, rows, cols>` | `operator-(const Matrix<T, rows, cols>&)` | `const Matrix<T, rows, cols>& matrix`-矩阵 | -           |
```C++
EasyMath::Matrix<float,2,2> matrixA {1,2,3,4};
std::cout << -matrixA << std::endl;
```
### 比较运算
| Access             | Modifier                                                                                                   | Return | Name                                                                    | Parameter                                                                                 | Description |
| ------------------ | ---------------------------------------------------------------------------------------------------------- | ------ | ----------------------------------------------------------------------- | ----------------------------------------------------------------------------------------- | ----------- |
| `private` `friend` | `template<typename T2, size_t rows2, size_t cols2, typename = std::enable_if_t<std::is_arithmetic_v<T2>>>` | `bool` | `operator==(const Matrix<T,rows,cols>&, const Matrix<T2,rows2,cols2>&)` | `const Matrix<T,rows,cols>& matrixA`-矩阵A, `const Matrix<T2,rows2,cols2>& matrixB`-矩阵B | -           |
| `private` `friend` | `template<typename T2, size_t rows2, size_t cols2, typename = std::enable_if_t<std::is_arithmetic_v<T2>>>` | `bool` | `operator!=(const Matrix<T,rows,cols>&, const Matrix<T2,rows2,cols2>&)` | `const Matrix<T,rows,cols>& matrixA`-矩阵A, `const Matrix<T2,rows2,cols2>& matrixB`-矩阵B | -           |
## 输出
| Access   | Modifier                                 | Return          | Name                                                     | Parameter                                                                    | Description                        |
| -------- | ---------------------------------------- | --------------- | -------------------------------------------------------- | ---------------------------------------------------------------------------- | ---------------------------------- |
| `public` | `virtual` `const`                        | `std::ostream&` | `print(std::ostream&)`                                   | `std::ostream& out`-输出流                                                   | 在控制台输出（这个函数为成员函数） |
| -        | `template<typename T, size_t dimension>` | `std::ostream&` | `operator<<(std::ostream&, const Matrix<T,rows, cols>&)` | `std::ostream& out`-输出流,<br>`const Matrix<T,rows, cols>& matrix`-输出矩阵 | 在控制台输出（这个函数为全局函数） |
```C++
EasyMath::Matrix<float,4,4> MatrixA{ 1, 2, 3, 4,
                               5, 6, 7, 8,
                               9,10,11,12,
                              13,14,15,16};

std::cout << MatrixA << std::endl;
```
输出
```
Matrix 4x4
-                            -
| 1.000  2.000  3.000  4.000 |
| 5.000  6.000  7.000  8.000 |
| 9.000 10.000 11.000 12.000 |
|13.000 14.000 15.000 16.000 |
-                            -
```
## 成员方法
| Access   | Modifier                                                                                                      | Return                          | Name                        | Parameter                                                           | Description                                                                                                                                             |
| -------- | ------------------------------------------------------------------------------------------------------------- | ------------------------------- | --------------------------- | ------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `public` | `const`                                                                                                       | `Matrix<T, rows - 1, cols - 1>` | `submatrix(size_t, size_t)` | `size_t row`-指定删除的行的索引,<br>`size_t col`-指定删除的列的索引 | 创建去掉指定行列的子矩阵                                                                                                                                |
| `public` | `const`                                                                                                       | `Matrix<T, cols, rows>`         | `transpose()`               | -                                                                   | 计算矩阵转置                                                                                                                                            |
| `public` | `template<size_t R = rows, size_t C = cols> typename std::enable_if_t<R == C, T>` `const`                     | `T`                             | `determinant()`             | -                                                                   | 计算行列式<br>**Matrix必须为方阵**                                                                                                                      |
| `public` | `template<size_t R = rows, size_t C = cols> typename std::enable_if_t<R == C && R >= 2, T>` `const`           | `T`                             | `cofactor(size_t, size_t)`  | `size_t row`-指定删除的行的索引,<br>`size_t col`-指定删除的列的索引 | 计算指定位置的代数余子式<br>代数余子式$$ = (-1)^{(i+j)} \cdot det(M_{ij})$$<br>其中$$M_{ij}$$是去掉第i行第j列后的子矩阵的行列式<br>**Matrix必须为方阵** |
| `public` | `template<size_t R = rows, size_t C = cols> typename std::enable_if_t<R == C, Matrix<T, rows, cols>>` `const` | `Matrix<T, rows, cols>`         | `cofactorMatrix()`          | -                                                                   | 计算代数余子式矩阵(伴随矩阵的转置) <br>**Matrix必须为方阵**                                                                                             |
| `public` | `template<size_t R = rows, size_t C = cols> typename std::enable_if_t<R == C, Matrix<T, rows, cols>>` `const` | `Matrix<T, rows, cols>`         | `adjugate()`                | -                                                                   | 计算伴随矩阵（代数余子式矩阵的转置）<br>**Matrix必须为方阵**                                                                                            |
| `public` | `template<size_t R = rows, size_t C = cols> typename std::enable_if_t<R == C, Matrix<T, rows, cols>>` `const` | `Matrix<T, rows, cols>`         | `inverse()`                 | -                                                                   | 计算逆矩阵（使用伴随矩阵方法）<br>**Matrix必须为方阵**                                                                                                  |
```C++
EasyMath::Matrix<double,4,3> matrixA {1,2,3,4,5,6,7,8,9,10,11,12};
EasyMath::Matrix<double,4,4> matrixB {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};

matrixA.submatrix(0,1);
matrixA.transpose();

matrixB.submatrix(1,2);
matrixB.transpose();
matrixB.determinant();
matrixB.cofactorMatrix();
matrixB.adjugate();
matrixB.inverse();
```
## 全局方法
| Access | Modifier                                                                               | Return            | Name                                                 | Parameter                                                                                                                                                   | Description                 |
| ------ | -------------------------------------------------------------------------------------- | ----------------- | ---------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------- |
| -      | `template<typename T, size_t N, typename = std::enable_if_t<std::is_arithmetic_v<T>>>` | `Matrix<T, N, N>` | `MTXIdentity()`                                      | -                                                                                                                                                           | 返回单位矩阵                |
| -      | `template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>`           | `Matrix<T, 4, 4>` | `MTXRotationX(const T&)`                             | `const T& radians`-旋转弧度                                                                                                                                 | 返回绕x轴旋转固定弧度的矩阵 |
| -      | `template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>`           | `Matrix<T, 4, 4>` | `MTXRotationY(const T&)`                             | `const T& radians`-旋转弧度                                                                                                                                 | 返回绕y轴旋转固定弧度的矩阵 |
| -      | `template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>`           | `Matrix<T, 4, 4>` | `MTXRotationZ(const T&)`                             | `const T& radians`-旋转弧度                                                                                                                                 | 返回绕z轴旋转固定弧度的矩阵 |
| -      | `template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>`           | `Matrix<T, 4, 4>` | `MTXTranslation(const T&, const T&, const T&, bool)` | `const T& x`-x方向上位移的距离,<br>`const T& y`-y方向上位移的距离,<br>`const T& z`-z方向上位移的距离,<br>`bool usedWithOrient = false`-是否用于平移方向向量 | 返回平移矩阵                |
| -      | `template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>`           | `Matrix<T, 4, 4>` | `MTXScale(const T&, const T&, const T&)`             | `const T& x`-x方向上的缩放单位,<br>`const T& y`-y方向上的缩放单位,<br>`const T& z`-z方向上的缩放单位                                                        | 返回缩放矩阵                |
| -      | `template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>`           | `Matrix<T, 4, 4>` | `MTXReflection(ReflectionPlane)`                     | `ReflectionPlane plane`-反射平面类型（`ReflectionPlane::YZ/XZ/XY/Origin`）                                                                                  | 返回基础反射矩阵（关于坐标平面或原点） |
| -      | `template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>`           | `Matrix<T, 4, 4>` | `MTXReflection(const Vector<T,3>&, const Vector<T,3>&)` | `const Vector<T,3>& normal`-平面法向量（自动归一化）,<br>`const Vector<T,3>& pointOnPlane`-平面上任意一点                                                    | 返回任意平面反射矩阵（Householder 反射） |

---

## `MTXReflection` 反射矩阵

### 函数签名

#### 重载 A：基础反射（枚举形式）
```cpp
template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
Matrix<T, 4, 4> MTXReflection(ReflectionPlane plane)
```

#### 重载 B：任意平面反射
```cpp
template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
Matrix<T, 4, 4> MTXReflection(
    const Vector<T, 3>& normal,
    const Vector<T, 3>& pointOnPlane)
```

### 参数说明

#### 重载 A 参数
| 参数 | 类型 | 描述 |
|------|------|------|
| `plane` | `ReflectionPlane` | 反射平面类型，可选值：<br>• `ReflectionPlane::YZ` - 关于 YZ 平面（X=0）反射，翻转 X 坐标<br>• `ReflectionPlane::XZ` - 关于 XZ 平面（Y=0）反射，翻转 Y 坐标<br>• `ReflectionPlane::XY` - 关于 XY 平面（Z=0）反射，翻转 Z 坐标<br>• `ReflectionPlane::Origin` - 关于原点反射，翻转所有坐标 |

#### 重载 B 参数
| 参数 | 类型 | 描述 |
|------|------|------|
| `normal` | `const Vector<T, 3>&` | 平面的法向量（**不需要归一化**，函数会自动归一化；零向量会触发断言） |
| `pointOnPlane` | `const Vector<T, 3>&` | 平面上任意一点（用于确定平面位置） |

### 返回值

- 返回 `Matrix<T, 4, 4>` 类型的 4×4 反射矩阵
- **行列式 det = -1**（反射会翻转手性）
- **3×3 旋转部分正交**：左上 3×3 块是 Householder 矩阵（正交）
- **自逆性**（满足 `M = M⁻¹`，因为反射两次等于原状）

### 数学原理

#### Householder 反射公式
任意平面反射矩阵的公式为：

```
M = I - 2 n n^T + 2 (n·P) n ⊗ e4
```

其中：
- `n` 是单位法向量
- `P` 是平面上一点
- `e4 = (0, 0, 0, 1)^T`
- `⊗` 表示外积

#### 矩阵展开

```
M = | 1-2nx²    -2nxny    -2nxnz     2nx·d |
    | -2nynx     1-2ny²    -2ny nz     2ny·d |
    | -2nznx    -2nzny     1-2nz²     2nz·d |
    |    0          0          0         1   |
```

其中 `d = n·P` 是平面到原点的有向距离。

#### 基础反射（枚举形式）

| 枚举值 | 几何意义 | 对角线矩阵 |
|--------|----------|------------|
| `YZ` | 关于 YZ 平面反射（X=0） | `diag(-1, 1, 1, 1)` |
| `XZ` | 关于 XZ 平面反射（Y=0） | `diag(1, -1, 1, 1)` |
| `XY` | 关于 XY 平面反射（Z=0） | `diag(1, 1, -1, 1)` |
| `Origin` | 关于原点反射 | `diag(-1, -1, -1, 1)` |

### 不变量与性质

1. **自逆性**：`M = M⁻¹`（反射矩阵的逆等于自身，因为反射两次等于原状）
2. **行列式**：`det(M) = -1`（反射会翻转手性）
3. **3×3 旋转部分正交**：左上 3×3 块 `R = I - 2 n n^T` 是 Householder 矩阵，满足 `R * R^T = I`
   - 注：完整 4×4 反射矩阵包含平移分量，整体不是 4×4 正交矩阵（这是仿射反射的正常行为）
4. **不动点**：平面上的点反射后位置不变

### 使用场景

#### 1. 平面镜渲染
```cpp
// 创建镜面反射矩阵
Vector3 mirrorNormal{0, 1, 0};        // 镜面法向量（水平面）
Vector3 mirrorPoint{0, 0, 0};         // 镜面经过原点
auto reflectMatrix = MTXReflection(mirrorNormal, mirrorPoint);

// 应用反射
auto reflectedModel = reflectMatrix * modelMatrix;
// 渲染镜像场景
```

#### 2. 镜像复制（Half-model 建模）
```cpp
// 只建模一半（如汽车右半侧），然后镜像生成完整模型
auto mirrorModel = MTXReflection(ReflectionPlane::XZ) * halfModel;
```

#### 3. 物理模拟（碰撞反弹）
```cpp
// 计算球撞墙后的反弹方向
Vector3 wallNormal{1, 0, 0};     // 墙的法向量
Vector3 wallPoint{10, 0, 0};     // 墙的位置
auto reflectWall = MTXReflection(wallNormal, wallPoint);

auto bouncedDirection = reflectWall * velocityVector;
```

#### 4. Stencil 阴影（平面阴影）
```cpp
// 在地面上渲染阴影
Vector3 groundNormal{0, 1, 0};
Vector3 groundPoint{0, -2, 0};
auto shadowMatrix = MTXReflection(groundNormal, groundPoint);
```

### 使用示例

#### 示例 1：基础反射（枚举形式）
```cpp
// 关于 YZ 平面反射（翻转 X）
auto reflectYZ = MTXReflection(ReflectionPlane::YZ);
// 输出：diag(-1, 1, 1, 1)

Vector4 point{2, 3, 4, 1};
auto reflected = reflectYZ * point;
// 结果：(-2, 3, 4, 1)
```

#### 示例 2：任意平面反射
```cpp
// 关于 Y = 5 平面反射
Vector3 normal{0, 1, 0};          // 法向量指向 Y 轴正方向
Vector3 pointOnPlane{0, 5, 0};    // 平面经过点 (0, 5, 0)

auto reflectMatrix = MTXReflection(normal, pointOnPlane);

Vector4 testPoint{2, 8, 4, 1};     // 测试点 (2, 8, 4)
auto reflectedPoint = reflectMatrix * testPoint;
// 结果：(2, 2, 4, 1) —— Y 坐标镜像：5 - (8-5) = 2
```

#### 示例 3：验证性质
```cpp
Vector3 normal{0, 1, 0};
Vector3 pointOnPlane{0, 5, 0};
auto reflectMat = MTXReflection(normal, pointOnPlane);

// 验证行列式 = -1
float det = reflectMat.determinant();  // det = -1.0f

// 验证 3×3 上三角块正交性（左上 3×3 块是 Householder 矩阵）
Matrix<float, 3, 3> upper3x3;
for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
        upper3x3(i, j) = reflectMat(i, j);
    }
}
auto orthoCheck = upper3x3 * upper3x3.transpose();
// orthoCheck ≈ MTXIdentity<float, 3>()（仅线性部分正交）

// 验证自逆性：M == M⁻¹
auto inverseMat = reflectMat.inverse();
bool isSelfInverse = (reflectMat == inverseMat);  // true
```

#### 示例 4：一致性验证（枚举 vs 通用）
```cpp
// 枚举形式：关于 YZ 平面反射
auto reflectYZ_enum = MTXReflection(ReflectionPlane::YZ);

// 通用形式：法向量 (-1, 0, 0) 过原点
Vector3 normal{-1, 0, 0};
Vector3 point{0, 0, 0};
auto reflectYZ_generic = MTXReflection(normal, point);

// 两者结果相同
bool isEqual = (reflectYZ_enum == reflectYZ_generic);  // true
```

### 注意事项

1. **法向量自动归一化**：`MTXReflection` 会自动归一化法向量，无需手动归一化
2. **零向量检查**：如果传入零法向量，会触发断言：`"Reflection normal vector cannot be zero-length"`
3. **负法向量等价**：传 `-n` 与传 `+n` 会得到相同的反射矩阵（公式中符号会抵消）
4. **矩阵存储顺序**：EasyMath 使用**行主序**（DirectX 风格），反射矩阵公式与存储顺序无关
5. **矩阵×向量**：`Matrix × Vector` 返回 `Matrix<T, 4, 1>`，需要用 `[0]`/`[1]`/`[2]` 访问分量
6. **整数类型支持**：虽然支持整数类型（满足 `is_arithmetic_v`），但实际应用应使用浮点类型

### 与 `Vector::reflect` 的关系

| 特性 | `MTXReflection` | `Vector::reflect` |
|------|-----------------|-------------------|
| 形式 | 矩阵形式 | 向量形式 |
| 适用场景 | 变换整个物体（模型矩阵） | 计算反射方向（光线、碰撞） |
| 参数 | 平面法向量 + 平面上一点 | 入射向量 + 平面法向量 |
| 返回值 | 4×4 变换矩阵 | 反射后的向量 |

两者互补：`MTXReflection` 用于物体级变换，`Vector::reflect` 用于向量级计算。

--- 

[![Static Badge](https://img.shields.io/badge/Back_to_Top-gray)](#matrix类 "点击返回顶部") | [![Static Badge](https://img.shields.io/badge/Back_to_README-gray)](../README.md "点击返回README")