Matrix API Reference
====================

Overview
--------
`EM::Matrix<T, rows, cols>` is a row-major, fixed-size matrix template. `T` must be arithmetic (`std::is_arithmetic_v<T>`). Debug builds use `assert` for bounds and precondition checks.

Template parameters
-------------------
- `T`: element type (float/double/int...)
- `rows`: number of rows (> 0)
- `cols`: number of columns (> 0)

Construction
------------
- `Matrix()` — zero-initialized storage
- `Matrix(std::initializer_list<T> values)` — row‑major order; size must be `rows * cols`

Element access
--------------
- Linear index (row‑major):
  - `T& operator[](size_t idx)` / `const T& operator[](size_t idx) const`
- 2D access:
  - `T& operator()(size_t r, size_t c)` / `const T& operator()(size_t r, size_t c) const`
- Dimensions:
  - `static constexpr size_t Rows()`
  - `static constexpr size_t Cols()`

Formatting
----------
- `std::ostream& print(std::ostream& out) const`
- `operator<<` forwards to `print` with fixed width and 3 decimal precision

Submatrix
---------
- `Matrix<T, rows-1, cols-1> submatrix(size_t r, size_t c) const`
  - Returns the minor with row `r` and column `c` removed

Determinant / Cofactor / Adjugate / Inverse
-------------------------------------------
- `T determinant() const` (square matrices only)
  - 1×1: the element itself
  - 2×2: closed form `ad - bc`
  - N×N: Laplace expansion along first row (recursive)
- `T cofactor(size_t r, size_t c) const`
- `Matrix<T, rows, cols> cofactorMatrix() const`
- `Matrix<T, rows, cols> adjugate() const` — transpose of the cofactor matrix
- `Matrix<T, rows, cols> inverse() const`
  - Requires `determinant() != 0`; uses `adjugate() / det`

Algebraic operations
--------------------
- In‑place:
  - `Matrix& operator+=(const Matrix& rhs)`
  - `Matrix& operator-=(const Matrix& rhs)`
  - `Matrix& operator*=(Scalar s)` (scalar multiply)
- Friends:
  - Scalar multiply: `matrix * s`, `s * matrix`
  - Add/Sub: `lhs + rhs`, `lhs - rhs`
  - Unary minus: `-matrix`
  - Matrix × Matrix: `lhs * rhs` (dimension‑checked)

Matrix × Vector
---------------
- `Matrix<T, rows, cols> * Vector<T, cols> -> Vector<T, rows>`
- `Vector<T, rows> * Matrix<T, rows, cols> -> Vector<T, cols>`

Transpose
---------
- `Matrix<T, cols, rows> transpose() const`

Helpers (global)
----------------
- Identity: `Matrix<T, N, N> MTXIdentity<T, N>()`
- 3D Rotations (radians):
  - `MTXRotationX<T, N>(T radians)`
  - `MTXRotationY<T, N>(T radians)`
  - `MTXRotationZ<T, N>(T radians)`
- 3D Translation:
  - `MTXTranslation<T, N>(T x, T y, T z, bool usedWithOrient = false)`
    - If `usedWithOrient` is true, sets (3,3) to 0 for composition with orientation matrices
- 3D Scale:
  - `MTXScale<T, N>(T sx, T sy, T sz)`

Storage & Conventions
---------------------
- Row‑major layout: `data[r*cols + c]`
- Initializer lists are row‑major
- Rotation angles are radians
- Debug assertions guard bounds, size compatibility, and non‑singular inverse

Complexity notes
----------------
- `determinant()` by Laplace expansion is factorial in worst case; suitable for small N (≤ 4)
- Matrix multiplication: O(rows * cols * inner)

Examples
--------
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

Planned extensions
------------------
- LU‑based determinant/inverse for performance
- Additional constructors, block operations, lightweight views


