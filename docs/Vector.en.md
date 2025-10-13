Vector API Reference
====================

Overview
--------
`EM::Vector<T, N>` is a dimension‑generic, fixed‑size vector template. `T` must be arithmetic. It supports element access, swizzles, arithmetic, norms, projections, interpolation, homogeneous conversions, and utilities for graphics math. Debug builds use `assert` for bounds and preconditions.

Template parameters
-------------------
- `T`: element type (float/double/int...)
- `N`: dimension (N > 0)

Construction
------------
- Default zero‑init: `Vector()`
- Fill with scalar: `explicit Vector(const T& value)`
- Initializer list: `Vector(std::initializer_list<T>)` (size must equal `N`)
- Copy/move: defaulted
- From `Quaternion<T>` (only when `N == 4`)

Element access
--------------
- `T& operator[](size_t i)` / `const T& operator[](size_t i) const`
- Bounds‑checked access: `at(size_t i)` / `const at(size_t i) const`
- Iterators: `begin()/end()` const and non‑const
- Raw data: `T* Data()` / `const T* Data() const`
- Size/introspection: `static constexpr size_t size()` / `getDimension()` / `Dimension()`

Swizzles
--------
- 1D component aliases via `VectorFilterDimension1D` (x,y,z,w / r,g,b,a with case variants)
- 2D swizzle to `Vector<T,2>` via `operator[](VectorFilterDimension2D)` (xy, xz, yz, xw, yw, zw, and color aliases)
- 3D swizzle to `Vector<T,3>` via `operator[](VectorFilterDimension3D)` (xyz, xyw, yzw, and color aliases)

Matrix conversions
------------------
- Column vector: `Matrix<T, N, 1> toColMatrix()`
- Row vector: `Matrix<T, 1, N> toRowMatrix()`

Arithmetic (element‑wise unless noted)
-------------------------------------
- In‑place: `+=`, `-=`, `*=`, `/=` with vector or scalar variants
- Friends: `+`, unary `-`, `-`, `*`, `/` with vector or scalar
- Comparisons: `==`, `!=` (floating‑point uses epsilon via `NEARZERO_THRESHOLD`)

Norms and normalization
-----------------------
- `T length(bool dimensionalityReduction = false) const`
- `constexpr T lengthSquared(bool dimensionalityReduction = false) const`
- `Vector normalized(bool dimensionalityReduction = false) const`
- `void normalize(bool dimensionalityReduction = false)`
- `bool isZero(T epsilon = T{NEARZERO_THRESHOLD}) const`
- `bool isNormalized(T epsilon = T{NEARZERO_THRESHOLD}, bool dimensionalityReduction = true) const`
  - When `dimensionalityReduction = true`, only the first 3 components contribute to the norm; higher dims are preserved during `normalized()`.

Products, distances, angles
---------------------------
- Dot: `T dot(const Vector& a, const Vector& b)`
- Cross (3D): `Vector<T,3> cross(const Vector<T,3>& a, const Vector<T,3>& b)`
- Cross2D (scalar): `T cross2D(const Vector<T,2>& a, const Vector<T,2>& b)`
- Distance: `T distance(const Vector& a, const Vector& b, bool dimensionalityReduction = true)`
- Distance squared: `T distanceSquared(const Vector& a, const Vector& b, bool dimensionalityReduction = true)`
- Angle: `T angle(const Vector& a, const Vector& b, bool dimensionalityReduction = true)` (returns radians)

Interpolation
-------------
- Linear: member `lerp(const Vector& other, T t) const`, free `lerp(a, b, t)`
- Spherical (unit vectors): `slerp(a, b, t)`

Projection & reflection
-----------------------
- `Vector project(const Vector& onto) const`
- `Vector project(const Vector& onto, bool dimensionalityReduction) const`
- `Vector reflect(const Vector& normal) const` (expects unit normal; asserts N ≥ 2)

Homogeneous coordinates
-----------------------
- To homogeneous (append w): `template<size_t newDim = N+1> Vector<T, newDim> toHomogeneous(T w = T{1}) const`
- From homogeneous (perspective divide on last component): `template<size_t newDim = N-1> Vector<T, newDim> fromHomogeneous() const`
  - Handles `w == 0` as direction (no divide)

Graphics helpers
----------------
- Translation matrix (for 3D vectors only): `toTranslationMatrix(bool usedWithOrient = false)` -> `Matrix<T,4,4>`
- Skew‑symmetric matrix: `skewSymmetric()` returns an antisymmetric `N×N` matrix (specialized for N=2,3)
- `skewSymmetric_2D()` -> `Matrix<T,3,3>` available only when `N == 2`

Output
------
- `std::ostream& print(std::ostream&) const` — prints `(x, y, z, ...) (V mode: xyzw)`

Common aliases
--------------
- `Vector2/3/4` (float)
- `Vector2f/3f/4f`, `Vector2d/3d/4d`, `Vector2i/3i/4i`

Examples
--------
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

Notes & caveats
---------------
- `dimensionalityReduction` reduces norm computations to the first 3 components; higher dimensions are preserved in `normalized()`
- Reflection expects a unit normal; call `normalized()` on the normal if unsure
- Comparisons of floating‑point vectors use a fixed epsilon (`NEARZERO_THRESHOLD`)

Planned extensions
------------------
- Additional swizzle packs, clamp/saturate helpers, random/unit generators
- SIMD specializations and small‑vector optimizations


