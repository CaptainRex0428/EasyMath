# `Matrix`类

[![Static Badge](https://img.shields.io/badge/Back_to_README-gray)](../README.md "点击返回README") | [![Static Badge](https://img.shields.io/badge/Language-zh-red)](Matrix.en.md "click to english version")

---

![Static Badge](https://img.shields.io/badge/Type-class-green)
![Static Badge](https://img.shields.io/badge/State-completed-blue)


Matrix是一个长度可变、类型可变的泛型矩阵

| Attribute     | Value                                                              |
| ------------- | ------------------------------------------------------------------ |
| **Namespace** | **`EM`**                                                           |
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
EM::Matrix<float,3,3> MatrixA;
EM::Matrix<float,4,4> MatrixB{ 1, 2, 3, 4,
                               5, 6, 7, 8,
                               9,10,11,12,
                              13,14,15,16};

EM::Vector<float,4> VectorA{1,2,3,4};
EM::Vector<float,4> VectorB{5,6,7,8};
EM::Vector<float,4> VectorC{9,10,11,12};
EM::Vector<float,4> VectorD{13,14,15,16};

EM::Matrix<float,4,4> MatrixC{VectorA,VectorB,VectorC,VectorD};

EM::Matrix<float,3,4> MatrixD_T{VectorA,VectorB,VectorC};
EM::Matrix<float,4,3> MatrixD = MatrixD_T.transpose();
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
EM::Vector<float,4> VectorA{1,2,3,4};
EM::Vector<float,4> VectorB{5,6,7,8};
EM::Vector<float,4> VectorC{9,10,11,12};
EM::Vector<float,4> VectorD{13,14,15,16};
EM::Matrix<float,4,4> MatrixC{VectorA,VectorB,VectorC,VectorD};

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
EM::Matrix<float,2,2> matrixA {1,2,3,4};
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
EM::Matrix<float,4,4> MatrixA{ 1, 2, 3, 4,
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
EM::Matrix<double,4,3> matrixA {1,2,3,4,5,6,7,8,9,10,11,12};
EM::Matrix<double,4,4> matrixB {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};

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

--- 

[![Static Badge](https://img.shields.io/badge/Back_to_Top-gray)](#matrix类 "点击返回顶部") | [![Static Badge](https://img.shields.io/badge/Back_to_README-gray)](../README.md "点击返回README")