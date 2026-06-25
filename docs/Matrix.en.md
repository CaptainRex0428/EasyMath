# Class `Matrix`

[![Static Badge](https://img.shields.io/badge/Back_to_README-gray)](../README.md "Click to back to README") | [![Static Badge](https://img.shields.io/badge/Language-en-red)](Matrix.md "点击切换中文版")

---

![Static Badge](https://img.shields.io/badge/Type-class-green)
![Static Badge](https://img.shields.io/badge/State-completed-blue)


`Matrix` is a generic matrix with variable length and element type

| Attribute     | Value                                                              |
| ------------- | ------------------------------------------------------------------ |
| **Namespace** | **`EasyMath`**                                                     |
| **File**      | [**`include\Matrix.h`**](../include/Matrix.h)                      |
| **Parent**    | -                                                                  |
| **Feature**   | ![Static Badge](https://img.shields.io/badge/Feature-template-red) |


## Declaration
```C++
template <typename T, size_t rows, size_t cols, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
class Matrix
```
### Templates
| Type     | Name   | Description                                                     |
| -------- | ------ | --------------------------------------------------------------- |
| typename | `T`    | Element type constrained by `std::is_arithmetic_v` (arithmetic) |
| size_t   | `rows` | Number of rows                                                  |
| size_t   | `cols` | Number of columns                                               |
### Constructors
| Access   | Modifier                                                                       | Return | Name                                            | Parameter                                                                                 | Description                                                                                  |
| -------- | ------------------------------------------------------------------------------ | ------ | ----------------------------------------------- | ----------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------- |
| `public` | -                                                                              |        | `Matrix()`                                      | -                                                                                         | Default construction, all elements initialized to 0                                          |
| `public` | -                                                                              |        | `Matrix(std::initializer_list<T>)`              | `std::initializer_list<T> InitializeList`-Initializer list of type `T`                    | If the number of scalar elements in the list is not equal to `rows * cols`, it will error    |
| `public` | -                                                                              |        | `Matrix(std::initializer_list<Vector<T,cols>>)` | `std::initializer_list<Vector<T,cols>> InitializeList`-Initializer list of `Vector<T,cols>` | If the number of vector elements in the list is not equal to `rows`, it will error           |
| `public` | `template<typename T2, typename = std::enable_if_t<std::is_arithmetic_v<T2>>>` |        | `Matrix(const Matrix<T2,rows,cols>&)`           | `const Matrix<T2,rows,cols>& other`-Matrix of the same dimension but different element type |                                                                                              |
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
## Data access
| Access   | Modifier                        | Return     | Name                        | Parameter                                                   | Description                                     |
| -------- | ------------------------------- | ---------- | --------------------------- | ----------------------------------------------------------- | ----------------------------------------------- |
| `public` | `virtual`                       | `T&`       | `operator[](size_t)`        | `size_t idx`-element index                                  | -                                               |
| `public` | `virtual` `const`               | `const T&` | `operator[](size_t)`        | `size_t idx`-element index                                  | -                                               |
| `public` | `virtual`                       | `T&`       | `at(size_t)`                | `size_t idx`-element index                                  | -                                               |
| `public` | `virtual` `const`               | `const T&` | `at(size_t)`                | `size_t idx`-element index                                  | -                                               |
| `public` | `virtual`                       | `T&`       | `operator()(size_t,size_t)` | `size_t row`-row index,<br>`size_t col`-column index        | -                                               |
| `public` | `virtual` `const`               | `const T&` | `operator()(size_t,size_t)` | `size_t row`-row index,<br>`size_t col`-column index        | -                                               |
| `public` | `static` `constexpr`            | `size_t`   | `Rows()`                    | -                                                           | Number of rows                                  |
| `public` | `static` `constexpr`            | `size_t`   | `Cols()`                    | -                                                           | Number of columns                               |
| `public` | `virtual` `noexcept`            | `T*`       | `Data()`                    | -                                                           | Return the pointer to the data array            |
| `public` | `virtual` `const` `noexcept`    | `const T*` | `Data()`                    | -                                                           | Return the pointer to the data array            |
| `public` | `virtual` `noexcept`            | `T*`       | `begin()`                   | -                                                           | Pointer to the first element of the data array  |
| `public` | `virtual` `const` `noexcept`    | `const T*` | `begin()`                   | -                                                           | Pointer to the first element of the data array  |
| `public` | `virtual` `noexcept`            | `T*`       | `end()`                     | -                                                           | End pointer (last element index + 1)            |
| `public` | `virtual` `const` `noexcept`    | `const T*` | `end()`                     | -                                                           | End pointer (last element index + 1)            |
| `public` | `static` `constexpr` `noexcept` | `size_t`   | `Size()`                    | -                                                           | Return the vector dimension (array dimension)   |
## Operations support
### Compound assignment
| Access   | Modifier  | Return                   | Name                                       | Parameter                                                    | Description |
| -------- | --------- | ------------------------ | ------------------------------------------ | ------------------------------------------------------------ | ----------- |
| `public` | `virtual` | `Matrix<T, rows, cols>&` | `operator+=(const Matrix<T, rows, cols>&)` | `const Matrix<T, rows, cols>& matrix`-Matrix of the same dim | -           |
| `public` | `virtual` | `Matrix<T, rows, cols>&` | `operator+=(const T&)`                     | `const T& scalar`-scalar                                     | -           |
| `public` | `virtual` | `Matrix<T, rows, cols>&` | `operator-=(const Matrix<T, rows, cols>&)` | `const Matrix<T, rows, cols>& matrix`-Matrix of the same dim | -           |
| `public` | `virtual` | `Matrix<T, rows, cols>&` | `operator-=(const T&)`                     | `const T& scalar`-scalar                                     | -           |
| `public` | `virtual` | `Matrix<T, rows, cols>&` | `operator*=(const Matrix<T, cols, cols>&)` | `const Matrix<T, cols, cols>& matrix`-Matrix (accumulation must be square) | -           |
| `public` | `virtual` | `Matrix<T, rows, cols>&` | `operator*=(const T&)`                     | `const T& scalar`-scalar                                     | -           |
| `public` | `virtual` | `Matrix<T, rows, cols>&` | `operator/=(const T&)`                     | `const T& scalar`-scalar                                     | -           |
### Binary operators
| Access             | Modifier                             | Return                   | Name                                                                     | Parameter                                                                              | Description                                                                        |
| ------------------ | ------------------------------------ | ------------------------ | ------------------------------------------------------------------------ | -------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------- |
| `private` `friend` | -                                    | `Matrix<T, rows, cols>`  | `operator*(const Matrix<T, rows, cols>&, const T&)`                      | `const Matrix<T, rows, cols>& matrix`-matrix,<br>`const T& scalar`-scalar              | -                                                                                  |
| `private` `friend` | -                                    | `Matrix<T, rows, cols>`  | `operator*(const T&, const Matrix<T, rows, cols>&)`                      | `const T& scalar`-scalar,<br>`const Matrix<T, rows, cols>& matrix`-matrix              | -                                                                                  |
| `private` `friend` | `template<typename T, size_t cols2>` | `Matrix<T, rows, cols2>` | `operator*(const Matrix<T, rows, cols>&, const Matrix<T, cols, cols2>&)` | `const Matrix<T, rows, cols>& lhs`-Matrix A,<br>`const Matrix<T, cols, cols2>& rhs`-Matrix B | -                                                                             |
| `private` `friend` | -                                    | `Matrix<T, 1, cols>`     | `operator*(const Matrix<T, rows, cols>&, const Vector<T, cols>&)`        | `const Matrix<T, rows, cols>& matrix`-matrix,<br>`const Vector<T, cols>& vec`-vector   | Matrix × Vector (column vector): Matrix<T, rows, cols> × Vector<T, cols> -> Vector<T, rows> |
| `private` `friend` | -                                    | `Matrix<T, rows, 1>`     | `operator*( const Vector<T, cols>&, const Matrix<T, rows, cols>&)`       | `const Vector<T, cols>& vec`-vector,<br>`const Matrix<T, rows, cols>& matrix`-matrix   | Vector (row vector) × Matrix: Vector<T, rows> × Matrix<T, rows, cols> -> Vector<T, cols>    |
| `private` `friend` | -                                    | `Matrix<T, rows, cols>`  | `operator+(const Matrix<T, rows, cols>&, const Matrix<T, rows, cols>&)`  | `const Matrix<T, rows, cols>& lhs`-Matrix A,<br>`const Matrix<T, rows, cols>& rhs`-Matrix B  | -                                                                              |
| `private` `friend` | -                                    | `Matrix<T, rows, cols>`  | `operator+(const Matrix<T, rows, cols>&, const T&)`                      | `const Matrix<T, rows, cols>& lhs`-matrix,<br>`const T& scalar`-scalar                  | -                                                                                  |
| `private` `friend` | -                                    | `Matrix<T, rows, cols>`  | `operator-(const Matrix<T, rows, cols>&, const Matrix<T, rows, cols>&)`  | `const Matrix<T, rows, cols>& lhs`-Matrix A,<br>`const Matrix<T, rows, cols>& rhs`-Matrix B  | -                                                                              |
| `private` `friend` | -                                    | `Matrix<T, rows, cols>`  | `operator－(const Matrix<T, rows, cols>&, const T&)`                     | `const Matrix<T, rows, cols>& lhs`-matrix,<br>`const T& scalar`-scalar                  | -                                                                                  |
```C++
EasyMath::Vector<float,4> VectorA{1,2,3,4};
EasyMath::Vector<float,4> VectorB{5,6,7,8};
EasyMath::Vector<float,4> VectorC{9,10,11,12};
EasyMath::Vector<float,4> VectorD{13,14,15,16};
EasyMath::Matrix<float,4,4> MatrixC{VectorA,VectorB,VectorC,VectorD};

std::cout << VectorA * MatrixC << std::endl;
std::cout << MatrixC * VectorA  << std::endl;
```
Output
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

### Unary
| Access             | Modifier | Return                  | Name                                      | Parameter                                  | Description |
| ------------------ | -------- | ----------------------- | ----------------------------------------- | ------------------------------------------ | ----------- |
| `private` `friend` | -        | `Matrix<T, rows, cols>` | `operator-(const Matrix<T, rows, cols>&)` | `const Matrix<T, rows, cols>& matrix`-matrix | -           |
```C++
EasyMath::Matrix<float,2,2> matrixA {1,2,3,4};
std::cout << -matrixA << std::endl;
```
### Comparison
| Access             | Modifier                                                                                                   | Return | Name                                                                    | Parameter                                                                                 | Description |
| ------------------ | ---------------------------------------------------------------------------------------------------------- | ------ | ----------------------------------------------------------------------- | ----------------------------------------------------------------------------------------- | ----------- |
| `private` `friend` | `template<typename T2, size_t rows2, size_t cols2, typename = std::enable_if_t<std::is_arithmetic_v<T2>>>` | `bool` | `operator==(const Matrix<T,rows,cols>&, const Matrix<T2,rows2,cols2>&)` | `const Matrix<T,rows,cols>& matrixA`-Matrix A, `const Matrix<T2,rows2,cols2>& matrixB`-Matrix B | -     |
| `private` `friend` | `template<typename T2, size_t rows2, size_t cols2, typename = std::enable_if_t<std::is_arithmetic_v<T2>>>` | `bool` | `operator!=(const Matrix<T,rows,cols>&, const Matrix<T2,rows2,cols2>&)` | `const Matrix<T,rows,cols>& matrixA`-Matrix A, `const Matrix<T2,rows2,cols2>& matrixB`-Matrix B | -     |
## Output
| Access   | Modifier                                 | Return          | Name                                                     | Parameter                                                                    | Description                                 |
| -------- | ---------------------------------------- | --------------- | -------------------------------------------------------- | ---------------------------------------------------------------------------- | ------------------------------------------- |
| `public` | `virtual` `const`                        | `std::ostream&` | `print(std::ostream&)`                                   | `std::ostream& out`-output                                                   | Output to console (this is a member function) |
| -        | `template<typename T, size_t dimension>` | `std::ostream&` | `operator<<(std::ostream&, const Matrix<T,rows, cols>&)` | `std::ostream& out`-output,<br>`const Matrix<T,rows, cols>& matrix`-matrix   | Output to console (this is a global function) |
```C++
EasyMath::Matrix<float,4,4> MatrixA{ 1, 2, 3, 4,
                               5, 6, 7, 8,
                               9,10,11,12,
                              13,14,15,16};

std::cout << MatrixA << std::endl;
```
Output
```
Matrix 4x4
-                            -
| 1.000  2.000  3.000  4.000 |
| 5.000  6.000  7.000  8.000 |
| 9.000 10.000 11.000 12.000 |
|13.000 14.000 15.000 16.000 |
-                            -
```
## Member methods
| Access   | Modifier                                                                                                      | Return                          | Name                        | Parameter                                                           | Description                                                                                                                                             |
| -------- | ------------------------------------------------------------------------------------------------------------- | ------------------------------- | --------------------------- | ------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `public` | `const`                                                                                                       | `Matrix<T, rows - 1, cols - 1>` | `submatrix(size_t, size_t)` | `size_t row`-row index to remove,<br>`size_t col`-column index to remove | Create the submatrix that removes the specified row and column                                                   |
| `public` | `const`                                                                                                       | `Matrix<T, cols, rows>`         | `transpose()`               | -                                                                   | Compute the transpose                                                                                                                                    |
| `public` | `template<size_t R = rows, size_t C = cols> typename std::enable_if_t<R == C, T>` `const`                     | `T`                             | `determinant()`             | -                                                                   | Compute the determinant<br>**Matrix must be square**                                                                                                     |
| `public` | `template<size_t R = rows, size_t C = cols> typename std::enable_if_t<R == C && R >= 2, T>` `const`           | `T`                             | `cofactor(size_t, size_t)`  | `size_t row`-row index,<br>`size_t col`-column index                | Compute the cofactor at the specified position<br>Cofactor $$ = (-1)^{(i+j)} \cdot det(M_{ij})$$<br>Where $$M_{ij}$$ is the determinant of the submatrix after removing row i and column j<br>**Matrix must be square** |
| `public` | `template<size_t R = rows, size_t C = cols> typename std::enable_if_t<R == C, Matrix<T, rows, cols>>` `const` | `Matrix<T, rows, cols>`         | `cofactorMatrix()`          | -                                                                   | Compute the cofactor matrix (transpose of the adjugate) <br>**Matrix must be square**                                                                    |
| `public` | `template<size_t R = rows, size_t C = cols> typename std::enable_if_t<R == C, Matrix<T, rows, cols>>` `const` | `Matrix<T, rows, cols>`         | `adjugate()`                | -                                                                   | Compute the adjugate matrix (transpose of the cofactor matrix) <br>**Matrix must be square**                                                             |
| `public` | `template<size_t R = rows, size_t C = cols> typename std::enable_if_t<R == C, Matrix<T, rows, cols>>` `const` | `Matrix<T, rows, cols>`         | `inverse()`                 | -                                                                   | Compute the inverse matrix (using the adjugate method) <br>**Matrix must be square**                                                                     |
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
## Global methods
| Access | Modifier                                                                               | Return            | Name                                                 | Parameter                                                                                                                                                   | Description                 |
| ------ | -------------------------------------------------------------------------------------- | ----------------- | ---------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------- |
| -      | `template<typename T, size_t N, typename = std::enable_if_t<std::is_arithmetic_v<T>>>` | `Matrix<T, N, N>` | `MTXIdentity()`                                      | -                                                                                                                                                           | Return the identity matrix  |
| -      | `template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>`           | `Matrix<T, 4, 4>` | `MTXRotationX(const T&)`                             | `const T& radians`-rotation radians                                                                                                                         | Return the rotation matrix around the x-axis |
| -      | `template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>`           | `Matrix<T, 4, 4>` | `MTXRotationY(const T&)`                             | `const T& radians`-rotation radians                                                                                                                         | Return the rotation matrix around the y-axis |
| -      | `template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>`           | `Matrix<T, 4, 4>` | `MTXRotationZ(const T&)`                             | `const T& radians`-rotation radians                                                                                                                         | Return the rotation matrix around the z-axis |
| -      | `template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>`           | `Matrix<T, 4, 4>` | `MTXTranslation(const T&, const T&, const T&, bool)` | `const T& x`-translation distance along x,<br>`const T& y`-translation distance along y,<br>`const T& z`-translation distance along z,<br>`bool usedWithOrient = false`-Whether used with orientation vectors | Return the translation matrix |
| -      | `template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>`           | `Matrix<T, 4, 4>` | `MTXScale(const T&, const T&, const T&)`             | `const T& x`-scale along x,<br>`const T& y`-scale along y,<br>`const T& z`-scale along z                                                                      | Return the scaling matrix   |

--- 

[![Static Badge](https://img.shields.io/badge/Back_to_Top-gray)](#class-matrix "点击返回顶部") | [![Static Badge](https://img.shields.io/badge/Back_to_README-gray)](../README.md "点击返回README")

