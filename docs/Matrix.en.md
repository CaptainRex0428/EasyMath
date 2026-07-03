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
| -      | `template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>`           | `Matrix<T, 4, 4>` | `MTXReflection(ReflectionPlane)`                     | `ReflectionPlane plane`-reflection plane type (`ReflectionPlane::YZ/XZ/XY/Origin`)                                                                               | Return the basic reflection matrix (about coordinate planes or origin) |
| -      | `template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>`           | `Matrix<T, 4, 4>` | `MTXReflection(const Vector<T,3>&, const Vector<T,3>&)` | `const Vector<T,3>& normal`-plane normal vector (auto-normalized),<br>`const Vector<T,3>& pointOnPlane`-any point on the plane                                      | Return the arbitrary plane reflection matrix (Householder reflection) |

---

## `MTXReflection` Reflection Matrix

### Function Signatures

#### Overload A: Basic Reflection (Enum Form)
```cpp
template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
Matrix<T, 4, 4> MTXReflection(ReflectionPlane plane)
```

#### Overload B: Arbitrary Plane Reflection
```cpp
template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
Matrix<T, 4, 4> MTXReflection(
    const Vector<T, 3>& normal,
    const Vector<T, 3>& pointOnPlane)
```

### Parameters

#### Overload A Parameters
| Parameter | Type | Description |
|-----------|------|-------------|
| `plane` | `ReflectionPlane` | Reflection plane type, options:<br>• `ReflectionPlane::YZ` - Reflect across YZ plane (X=0), flip X coordinate<br>• `ReflectionPlane::XZ` - Reflect across XZ plane (Y=0), flip Y coordinate<br>• `ReflectionPlane::XY` - Reflect across XY plane (Z=0), flip Z coordinate<br>• `ReflectionPlane::Origin` - Reflect across origin, flip all coordinates |

#### Overload B Parameters
| Parameter | Type | Description |
|-----------|------|-------------|
| `normal` | `const Vector<T, 3>&` | Plane normal vector (**does not need to be normalized**, function will auto-normalize; zero vector triggers assertion) |
| `pointOnPlane` | `const Vector<T, 3>&` | Any point on the plane (used to determine plane position) |

### Return Value

- Returns `Matrix<T, 4, 4>` type 4×4 reflection matrix
- **Determinant det = -1** (reflection flips chirality)
- **3×3 rotation part is orthogonal**: The upper-left 3×3 block is a Householder matrix (orthogonal)
- **Self-inverse** (satisfies `M = M⁻¹`, since reflecting twice returns to original state)

### Mathematical Principles

#### Householder Reflection Formula
The formula for an arbitrary plane reflection matrix is:

```
M = I - 2 n n^T + 2 (n·P) n ⊗ e4
```

Where:
- `n` is the unit normal vector
- `P` is a point on the plane
- `e4 = (0, 0, 0, 1)^T`
- `⊗` denotes outer product

#### Matrix Expansion

```
M = | 1-2nx²    -2nxny    -2nxnz     2nx·d |
    | -2nynx     1-2ny²    -2ny nz     2ny·d |
    | -2nznx    -2nzny     1-2nz²     2nz·d |
    |    0          0          0         1   |
```

Where `d = n·P` is the signed distance from the plane to the origin.

#### Basic Reflection (Enum Form)

| Enum Value | Geometric Meaning | Diagonal Matrix |
|------------|-------------------|-----------------|
| `YZ` | Reflect across YZ plane (X=0) | `diag(-1, 1, 1, 1)` |
| `XZ` | Reflect across XZ plane (Y=0) | `diag(1, -1, 1, 1)` |
| `XY` | Reflect across XY plane (Z=0) | `diag(1, 1, -1, 1)` |
| `Origin` | Reflect across origin | `diag(-1, -1, -1, 1)` |

### Invariants and Properties

1. **Self-inverse**: `M = M⁻¹` (reflection matrix's inverse equals itself, since reflecting twice returns to original state)
2. **Determinant**: `det(M) = -1` (reflection flips chirality)
3. **3×3 rotation part is orthogonal**: The upper-left 3×3 block `R = I - 2 n n^T` is a Householder matrix, satisfying `R * R^T = I`
   - Note: The full 4×4 reflection matrix contains a translation component and is not a 4×4 orthogonal matrix as a whole (this is normal behavior for an affine reflection)
4. **Fixed points**: Points on the plane remain unchanged after reflection

### Use Cases

#### 1. Planar Mirror Rendering
```cpp
// Create mirror reflection matrix
Vector3 mirrorNormal{0, 1, 0};        // Mirror normal vector (horizontal plane)
Vector3 mirrorPoint{0, 0, 0};         // Mirror passes through origin
auto reflectMatrix = MTXReflection(mirrorNormal, mirrorPoint);

// Apply reflection
auto reflectedModel = reflectMatrix * modelMatrix;
// Render mirrored scene
```

#### 2. Mirror Duplication (Half-model Modeling)
```cpp
// Only model one half (e.g., right side of car), then mirror to generate complete model
auto mirrorModel = MTXReflection(ReflectionPlane::XZ) * halfModel;
```

#### 3. Physics Simulation (Collision Bounce)
```cpp
// Calculate bounce direction after ball hits wall
Vector3 wallNormal{1, 0, 0};     // Wall normal vector
Vector3 wallPoint{10, 0, 0};     // Wall position
auto reflectWall = MTXReflection(wallNormal, wallPoint);

auto bouncedDirection = reflectWall * velocityVector;
```

#### 4. Stencil Shadows (Planar Shadows)
```cpp
// Render shadow on ground
Vector3 groundNormal{0, 1, 0};
Vector3 groundPoint{0, -2, 0};
auto shadowMatrix = MTXReflection(groundNormal, groundPoint);
```

### Usage Examples

#### Example 1: Basic Reflection (Enum Form)
```cpp
// Reflect across YZ plane (flip X)
auto reflectYZ = MTXReflection(ReflectionPlane::YZ);
// Output: diag(-1, 1, 1, 1)

Vector4 point{2, 3, 4, 1};
auto reflected = reflectYZ * point;
// Result: (-2, 3, 4, 1)
```

#### Example 2: Arbitrary Plane Reflection
```cpp
// Reflect across Y = 5 plane
Vector3 normal{0, 1, 0};          // Normal vector pointing along positive Y axis
Vector3 pointOnPlane{0, 5, 0};    // Plane passes through point (0, 5, 0)

auto reflectMatrix = MTXReflection(normal, pointOnPlane);

Vector4 testPoint{2, 8, 4, 1};     // Test point (2, 8, 4)
auto reflectedPoint = reflectMatrix * testPoint;
// Result: (2, 2, 4, 1) — Y coordinate mirrored: 5 - (8-5) = 2
```

#### Example 3: Verify Properties
```cpp
Vector3 normal{0, 1, 0};
Vector3 pointOnPlane{0, 5, 0};
auto reflectMat = MTXReflection(normal, pointOnPlane);

// Verify determinant = -1
float det = reflectMat.determinant();  // det = -1.0f

// Verify 3×3 upper block orthogonality (left-upper 3×3 is a Householder matrix)
Matrix<float, 3, 3> upper3x3;
for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
        upper3x3(i, j) = reflectMat(i, j);
    }
}
auto orthoCheck = upper3x3 * upper3x3.transpose();
// orthoCheck ≈ MTXIdentity<float, 3>() (only the linear part is orthogonal)

// Verify self-inverse: M == M⁻¹
auto inverseMat = reflectMat.inverse();
bool isSelfInverse = (reflectMat == inverseMat);  // true
```

#### Example 4: Consistency Verification (Enum vs Generic)
```cpp
// Enum form: reflect across YZ plane
auto reflectYZ_enum = MTXReflection(ReflectionPlane::YZ);

// Generic form: normal (-1, 0, 0) passing through origin
Vector3 normal{-1, 0, 0};
Vector3 point{0, 0, 0};
auto reflectYZ_generic = MTXReflection(normal, point);

// Both produce identical results
bool isEqual = (reflectYZ_enum == reflectYZ_generic);  // true
```

### Notes

1. **Auto-normalization**: `MTXReflection` automatically normalizes the normal vector; manual normalization is not required
2. **Zero vector check**: Passing a zero normal vector triggers assertion: `"Reflection normal vector cannot be zero-length"`
3. **Negative normal equivalence**: Passing `-n` produces the same reflection matrix as `+n` (signs cancel out in the formula)
4. **Matrix storage order**: EasyMath uses **row-major order** (DirectX style), reflection matrix formula is storage-order independent
5. **Matrix × Vector**: `Matrix × Vector` returns `Matrix<T, 4, 1>`, access components with `[0]`/`[1]`/`[2]`
6. **Integer type support**: While integer types are supported (satisfy `is_arithmetic_v`), floating-point types should be used in practice

### Relationship with `Vector::reflect`

| Feature | `MTXReflection` | `Vector::reflect` |
|----------|-----------------|-------------------|
| Form | Matrix form | Vector form |
| Use case | Transform entire object (model matrix) | Calculate reflection direction (lighting, collision) |
| Parameters | Plane normal + point on plane | Incident vector + plane normal |
| Return value | 4×4 transformation matrix | Reflected vector |

They are complementary: `MTXReflection` for object-level transforms, `Vector::reflect` for vector-level calculations.

--- 

[![Static Badge](https://img.shields.io/badge/Back_to_Top-gray)](#class-matrix "点击返回顶部") | [![Static Badge](https://img.shields.io/badge/Back_to_README-gray)](../README.md "点击返回README")

