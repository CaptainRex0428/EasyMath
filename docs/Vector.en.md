# Class `Vector`

[![Static Badge](https://img.shields.io/badge/Back_to_README-gray)](../README.md "Click to back to README") | [![Static Badge](https://img.shields.io/badge/language-en-red)](Vector.md "点击切换中文版")

--- 

![Static Badge](https://img.shields.io/badge/type-class-green)
![Static Badge](https://img.shields.io/badge/state-completed-blue)


`Vector` is a generic vector class with variable length and element type

| Attribute     | Value                                                                                                                                  |
| ------------- | -------------------------------------------------------------------------------------------------------------------------------------- |
| **Namespace** | **`EM`**                                                                                                                               |
| **File**      | [**`include\vector.h`**](../include/Vector.h)                                                                                          |
| **Parent**    | -                                                                                                                                      |
| **Feature**   | ![Static Badge](https://img.shields.io/badge/Feature-virtual-green) ![Static Badge](https://img.shields.io/badge/Feature-template-red) |

## Declaration

```C++
template<typename T, size_t dimension, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
class Vector
```

### Templates

| Type     | Name        | Description                                                     |
| -------- | ----------- | --------------------------------------------------------------- |
| typename | `T`         | Element type constrained by `typename = std::enable_if_t<>>` (arithmetic) |
| size_t   | `dimension` | Vector dimension/length                                         |

### Constructors

| Access   | Modifier   | Return | Name                                           | Parameter                                                 | Description                                                      |
| -------- | ---------- | ------ | ---------------------------------------------- | --------------------------------------------------------- | ---------------------------------------------------------------- |
| `public` | -          |        | `Vector()`                                     | -                                                         | Default construction, all elements initialized to 0              |
| `public` | `explicit` |        | `Vector(const T& value)`                       | Single value `const T& value` of type `T`                 | Initialize all vector elements with the single input value       |
| `public` | -          |        | `Vector(std::initializer_list<T> InitializeList)` | `std::initializer_list<T> InitializeList`                | Initialize from list of `T`; if size != `dimension`, it errors   |

```C++
EM::Vector<float,2> vectorA;
EM::Vector<float,3> vectorB(1);
EM::Vector<float,4> vectorC{ 1, 2, 3, 4};
```

## Data access

### General data access

| Access   | Modifier                                 | Return     | Name                 | Parameter                | Description                                  |
| -------- | ---------------------------------------- | ---------- | -------------------- | ------------------------ | -------------------------------------------- |
| `public` | `virtual`                                | `T&`       | `operator[](size_t)` | `size_t idx` element index | -                                            |
| `public` | `virtual` `const`                        | `const T&` | `operator[](size_t)` | `size_t idx` element index | -                                            |
| `public` | `virtual`                                | `T&`       | `at(size_t)`         | `size_t idx` element index | -                                            |
| `public` | `virtual` `const`                        | `const T&` | `at(size_t)`         | `size_t idx` element index | -                                            |
| `public` | `virtual` `noexcept`                     | `T*`       | `Data()`             | -                        | Return the pointer to the data array         |
| `public` | `virtual` `const` `noexcept`             | `const T*` | `Data()`             | -                        | Return the pointer to the data array         |
| `public` | `virtual` `noexcept`                     | `T*`       | `begin()`            | -                        | Pointer to the first element of the data     |
| `public` | `virtual` `const` `noexcept`             | `const T*` | `begin()`            | -                        | Pointer to the first element of the data     |
| `public` | `virtual` `noexcept`                     | `T*`       | `end()`              | -                        | End pointer (last element index + 1)         |
| `public` | `virtual` `const` `noexcept`             | `const T*` | `end()`              | -                        | End pointer (last element index + 1)         |
| `public` | `virtual` `constexpr` `const` `noexcept` | `size_t`   | `getDimension()`     | -                        | Return the vector dimension (array dimension) |
| `public` | `static` `constexpr` `noexcept`          | `size_t`   | `Size()`             | -                        | Return the vector dimension (array dimension) |
| `public` | `static` `constexpr` `noexcept`          | `size_t`   | `Dimension()`        | -                        | Return the vector dimension (array dimension) |


### Swizzle
Element access supports using xyzw or rgba like in GLSL/HLSL to access and modify elements
```C++
EM::Vector<float,4> vectorA{1,2,3,4};
std::cout << vectorA << std::endl;

vectorA.x = 5;
vectorA.zy = EM::Vector<float,2> {6,7};
std::cout << vectorA << std::endl;

std::cout << vectorA.xyz << std::endl;
```
```
// Output
(1, 2, 3, 4) (V mode: xyzw)
(5, 7, 6, 4) (V mode: xyzw)
(5, 7, 6) (V mode: xyz)
```

## Operations

### Compound assignment
| Access   | Modifier  | Return                  | Name                                      | Parameter                                                   | Description |
| -------- | --------- | ----------------------- | ----------------------------------------- | ----------------------------------------------------------- | ----------- |
| `public` | `virtual` | `Vector<T, dimension>&` | `operator=(const Vector<T, dimension>&)`  | `const Vector<T, dimension>& other` another vector of same type and dimension | -           |
| `public` | `virtual` | `Vector<T, dimension>&` | `operator=(Vector<T, dimension>&)`        | `Vector<T, dimension>& other` another vector of same type and dimension       | -           |
| `public` | `virtual` | `Vector<T, dimension>&` | `operator+=(const Vector<T, dimension>&)` | `const Vector<T, dimension>& other` another vector of same type and dimension | -           |
| `public` | `virtual` | `Vector<T, dimension>&` | `operator+=(const T&)`                    | `const T& scalar` constant of the same type                 | -           |
| `public` | `virtual` | `Vector<T, dimension>&` | `operator-=(const Vector<T, dimension>&)` | `const Vector<T, dimension>& other` another vector of same type and dimension | -           |
| `public` | `virtual` | `Vector<T, dimension>&` | `operator-=(const T&)`                    | `const T& scalar` constant of the same type                 | -           |
| `public` | `virtual` | `Vector<T, dimension>&` | `operator*=(const Vector<T, dimension>&)` | `const Vector<T, dimension>& other` another vector of same type and dimension | -           |
| `public` | `virtual` | `Vector<T, dimension>&` | `operator*=(const T&)`                    | `const T& scalar` constant of the same type                 | -           |
| `public` | `virtual` | `Vector<T, dimension>&` | `operator/=(const Vector<T, dimension>&)` | `const Vector<T, dimension>& other` another vector of same type and dimension | -           |
| `public` | `virtual` | `Vector<T, dimension>&` | `operator/=(const T&)`                    | `const T& scalar` constant of the same type                 | -           |

```C++
EM::Vector<float,2> vectorA{ 1, 2};

vectorA += 1;
std::cout << vectorA << std::endl;

vectorA += EM::Vector<float,2>{1,2};
std::cout << vectorA << std::endl;

vectorA -= 1;
std::cout << vectorA << std::endl;
	
vectorA -= EM::Vector<float,2>{1,2};
std::cout << vectorA << std::endl;

vectorA *= 3;
std::cout << vectorA << std::endl;
	
vectorA *= EM::Vector<float,2>{2,4};
std::cout << vectorA << std::endl;

vectorA /= 3;
std::cout << vectorA << std::endl;
	
vectorA /= EM::Vector<float,2>{2,4};
std::cout << vectorA << std::endl;
```

### Common operators
| Access             | Modifier | Return                 | Name                                                                  | Parameter                                                                       | Description |
| ------------------ | -------- | ---------------------- | --------------------------------------------------------------------- | ------------------------------------------------------------------------------- | ----------- |
| `private` `friend` | -        | `Vector<T, dimension>` | `operator+(const Vector<T, dimension>&, const Vector<T, dimension>&)` | `const Vector<T, dimension>& lhs` vector,<br>`const Vector<T, dimension>& rhs` vector | -           |
| `private` `friend` | -        | `Vector<T, dimension>` | `operator+(const Vector<T, dimension>&, const T&)`                    | `const Vector<T, dimension>& vec` vector,<br>`const T& scalar` scalar                 | -           |
| `private` `friend` | -        | `Vector<T, dimension>` | `operator+(const T&，const Vector<T, dimension>&)`                    | `const T& scalar` scalar，<br>`const Vector<T, dimension>& vec` vector                | -           |
| `private` `friend` | -        | `Vector<T, dimension>` | `operator-(const Vector<T, dimension>&, const Vector<T, dimension>&)` | `const Vector<T, dimension>& lhs` vector,<br>`const Vector<T, dimension>& rhs` vector | -           |
| `private` `friend` | -        | `Vector<T, dimension>` | `operator-(const Vector<T, dimension>&, const T&)`                    | `const Vector<T, dimension>& vec` vector,<br>`const T& scalar` scalar                 | -           |
| `private` `friend` | -        | `Vector<T, dimension>` | `operator*(const Vector<T, dimension>&, const Vector<T, dimension>&)` | `const Vector<T, dimension>& lhs` vector,<br>`const Vector<T, dimension>& rhs` vector | -           |
| `private` `friend` | -        | `Vector<T, dimension>` | `operator*(const Vector<T, dimension>&, const T&)`                    | `const Vector<T, dimension>& vec` vector,<br>`const T& scalar` scalar                 | -           |
| `private` `friend` | -        | `Vector<T, dimension>` | `operator*(const T&, const Vector<T, dimension>&)`                    | `const T& scalar` scalar,<br>`const Vector<T, dimension>& vec` vector                 | -           |
| `private` `friend` | -        | `Vector<T, dimension>` | `operator/(const Vector<T, dimension>&, const Vector<T, dimension>&)` | `const Vector<T, dimension>& lhs` vector,<br>`const Vector<T, dimension>& rhs` vector | -           |
| `private` `friend` | -        | `Vector<T, dimension>` | `operator/(const Vector<T, dimension>&, const T&)`                    | `const Vector<T, dimension>& vec` vector,<br>`const T& scalar` scalar                 | -           |

```C++
EM::Vector<float,2> vectorA{ 1, 2};
EM::Vector<float,2> vectorB{ 3, 4};

std::cout << vectorA+vectorB << std::endl;
std::cout << vectorA-vectorB << std::endl;
std::cout << vectorA*vectorB << std::endl;
std::cout << vectorA/vectorB << std::endl;
```

### Unary
| Access             | Modifier | Return                 | Name                                     | Parameter                             | Description      |
| ------------------ | -------- | ---------------------- | ---------------------------------------- | ------------------------------------- | ---------------- |
| `private` `friend` | -        | `Vector<T, dimension>` | `operator-(const Vector<T, dimension>&)` | `const Vector<T, dimension>& vec`     | Vector supports direct negation |

```C++
EM::Vector<float,2> vectorA{ 1, 2};
std::cout << -vectorA << std::endl;
```

### Comparison
| Access             | Modifier | Return | Name                                                                   | Parameter                                                                       | Description |
| ------------------ | -------- | ------ | ---------------------------------------------------------------------- | ------------------------------------------------------------------------------- | ----------- |
| `private` `friend` | -        | `bool` | `operator==(const Vector<T, dimension>&, const Vector<T, dimension>&)` | `const Vector<T, dimension>& lhs` vector,<br>`const Vector<T, dimension>& rhs` vector | -           |
| `private` `friend` | -        | `bool` | `operator!=(const Vector<T, dimension>&, const Vector<T, dimension>&)` | `const Vector<T, dimension>& lhs` vector,<br>`const Vector<T, dimension>& rhs` vector | -           |


## Output
| Access   | Modifier                                 | Return          | Name                                                     | Parameter                 | Description                        |
| -------- | ---------------------------------------- | --------------- | -------------------------------------------------------- | ------------------------- | ---------------------------------- |
| `public` | `virtual`                                | `std::ostream&` | `print(std::ostream&)`                                   | `std::ostream& out` output | Output to console (this is a member function) |
| -        | `template<typename T, size_t dimension>` | `std::ostream&` | `operator<<(std::ostream&, const Vector<T,dimension>&)<` | `std::ostream& out` output | Output to console (this is a global function)  |

```C++
EM::Vector<float,2> vectorA{8,9};
std::cout << vectorA << std::endl;
```

## Member methods
| Access   | Modifier                                                | Return                            | Name                                   | Parameter                                                                                                                               | Description                                          |
| -------- | ------------------------------------------------------- | --------------------------------- | -------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------- |
| `public` | `virtual`                                               | `Matrix<T, dimension, 1>`         | `toColMatrix()`                        | -                                                                                                                                       | Represent the vector as a column matrix             |
| `public` | `virtual`                                               | `Matrix<T, 1, dimension>`         | `toRowMatrix()`                        | -                                                                                                                                       | Represent the vector as a row matrix                |
| `public` | `virtual` `constexpr`                                   | `T`                               | `lengthSquared(bool)`                  | `bool dimensionalityReduction = false` Whether to perform dimensionality reduction<br>Only the first 3 dimensions are computed           | Compute the squared length of the vector            |
| `public` | `virtual` `const` `noexcept`                            | `T`                               | `length(bool)`                         | `bool dimensionalityReduction = false` Whether to perform dimensionality reduction<br>Only the first 3 dimensions are computed           | Compute the length of the vector                    |
| `public` | `virtual` `const`                                       | `Vector<T, dimension>`            | `normalized(bool)`                     | `bool dimensionalityReduction = false` Whether to perform dimensionality reduction<br>Only the first 3 dimensions are computed           | Compute the normalized result of the vector         |
| `public` | `virtual`                                               | `Vector<T, dimension>`                            | `normalize(bool)`                      | `bool dimensionalityReduction = false` Whether to perform dimensionality reduction<br>Only the first 3 dimensions are computed           | Normalize the vector                                |
| `public` | `virtual` `const` `noexcept`                            | `bool`                            | `isNormalized(T,bool)`                 | `bool dimensionalityReduction = false` Whether to perform dimensionality reduction<br>Only the first 3 dimensions are computed<br>`T epsilon = T{NEARZERO_THRESHOLD}` Floating-point comparison epsilon | Whether it is a unit vector                         |
| `public` | `virtual` `const` `noexcept`                            | `bool`                            | `isZero(T,bool)`                       | `T epsilon = T{NEARZERO_THRESHOLD}` Floating-point comparison epsilon                                                                     | Whether it is a zero vector                         |
| `public` | `virtual` `const` `noexcept`                            | `Vector<T, dimension>`            | `lerp(const Vector<T, dimension>&, T)` | `const Vector<T, dimension>& other` interpolation vector,<br> `T t` interpolation weight                                                  | Interpolation function                              |
| `public` | `virtual` `const`                                       | `Vector<T, dimension>`            | `reflect(const Vector<T, dimension>&)` | `const Vector<T, dimension>& normal` Normal vector (must be normalized)                                                                  | Reflection vector based on a normal                 |
| `public` | `virtual` `const`                                       | `Vector`                          | `project(const Vector&)`               | `const Vector& onto` Projection target                                                                                                   | Projection vector onto a certain vector             |
| `public` | `virtual` `const`                                       | `Vector`                          | `project(const Vector&, bool)`         | `const Vector& onto` Projection target,<br>`bool dimensionalityReduction = false` Whether to perform dimensionality reduction<br>Only the first 3 dimensions are computed | Projection vector onto a certain vector             |
| `public` | `template<size_t newDim = dimension + 1>` `const`       | `Vector<T, newDim>`               | `toHomogeneous(T)`                     | `T w = T{ 1 }` Homogeneous coordinate w component (w=1 for points, w=0 for directions)                                                   | Convert current coordinates to homogeneous          |
| `public` | `template<size_t newDim = dimension - 1>` `const`       | `Vector<T, newDim>`               | `fromHomogeneous()`                    | -                                                                                                                                       | Convert the current vector from homogeneous         |
| `public` | `template<typename = std::enable_if_t<dimension == 3>>` | `Matrix<T, 4, 4>`                 | `toTranslationMatrix(bool)`            | `bool usedWithOrient = false` Whether used for moving direction vectors (can be optimized away)                                          | Convert the current vector to a translation matrix  |
| `public` | `const`                                                 | `Matrix<T, dimension, dimension>` | `skewSymmetric()`                      | -                                                                                                                                       | Obtain the skew-symmetric matrix of the vector      |
| `public` | `template<size_t D = dimension>` `const`                | `Matrix<T, 3, 3>`                 | `skewSymmetric_2D()`                   | -                                                                                                                                       | Obtain the skew-symmetric matrix of current 2D vector (not accessible for non-2D) |


```C++
EM::Vector<float,2> vectorA{ 1, 2};
EM::Vector<float,3> vectorB{ 3, 4, 5};
EM::Vector<float,4> vectorC{ 6, 7, 8, 9};

std::cout << vectorB.toColMatrix() << std::endl;
std::cout << vectorB.toRowMatrix() << std::endl;
std::cout << vectorB.length() << std::endl;
std::cout << vectorB.normalized() << std::endl;
std::cout << vectorB.isNormalized() << std::endl;
std::cout << vectorB.isZero() << std::endl;
std::cout << vectorB.lerp(EM::Vector<float,3>{1,1,1},0.5f) << std::endl;
std::cout << vectorB.reflect(EM::Vector<float,3>{1,1,1}.normalize()) << std::endl;
std::cout << vectorB.project(EM::Vector<float,3>{1,1,1}.normalize()) << std::endl;
std::cout << vectorB.toHomogeneous(1) << std::endl;
std::cout << vectorB.toTranslationMatrix() << std::endl;
std::cout << vectorB.skewSymmetric() << std::endl;

std::cout << vectorA.skewSymmetric_2D() << std::endl;
std::cout << vectorC.fromHomogeneous() << std::endl;
```


## Global methods
| Access | Modifier                         | Return         | Name                                                              | Parameter                                                                                                                                          | Description  |
| ------ | -------------------------------- | -------------- | ----------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------- | ------------ |
| -      | `template<typename T, size_t D>` | `T`            | `dot(const Vector<T, D>&, const Vector<T, D>&)`                   | `const Vector<T, D>& a` Vector A,<br> `const Vector<T, D>& b` Vector B                                                                             | Dot product   |
| -      | `template<typename T>`           | `T`            | `cross(const Vector<T, 3>&, const Vector<T, 3>&)`                 | `const Vector<T, 3>& a` 3D Vector A,<br> `const Vector<T, 3>& b` 3D Vector B                                                                       | 3D cross      |
| -      | `template<typename T>`           | `T`            | `cross(const Vector<T, 2>&, const Vector<T, 2>&)`                 | `const Vector<T, 2>& a` 2D Vector A,<br> `const Vector<T, 2>& b` 2D Vector B                                                                       | 2D cross      |
| -      | `template<typename T, size_t D>` | `T`            | `distance(const Vector<T, D>&, const Vector<T, D>&, bool)`        | `const Vector<T, D>& a` Vector A,<br> `const Vector<T, D>& b` Vector B,<br> `bool dimensionalityReduction=true` Whether to perform dimensionality reduction<br>Only the first 3 dimensions are computed | Distance      |
| -      | `template<typename T, size_t D>` | `T`            | `distanceSquared(const Vector<T, D>&, const Vector<T, D>&, bool)` | `const Vector<T, D>& a` Vector A,<br> `const Vector<T, D>& b` Vector B,<br> `bool dimensionalityReduction=true` Whether to perform dimensionality reduction<br>Only the first 3 dimensions are computed | Distance squared |
| -      | `template<typename T, size_t D>` | `Vector<T, D>` | `lerp(const Vector<T, D>&, const Vector<T, D>&, T)`               | `const Vector<T, D>& a` Vector A,<br> `const Vector<T, D>& b` Vector B,<br> `T t` interpolation weight                                              | Linear interpolation |
| -      | `template<typename T, size_t D>` | `Vector<T, D>` | `slerp(const Vector<T, D>&, const Vector<T, D>&, T)`              | `const Vector<T, D>& a` Vector A,<br> `const Vector<T, D>& b` Vector B,<br> `T t` interpolation weight                                              | Spherical linear interpolation |

```C++
EM::Vector<float,3> vectorA{ 1, 2, 3};
EM::Vector<float,3> vectorB{ 4, 5, 6};

std::cout << dot(vectorA,vectorB) << std::endl;
std::cout << cross(vectorA,vectorB) << std::endl;
std::cout << distance(vectorA,vectorB) << std::endl;
std::cout << lerp(vectorA,vectorB,0.5f) << std::endl;
std::cout << slerp(vectorA,vectorB,0.5f) << std::endl;
std::cout << lerp(vectorA.normalize(),vectorB.normalize(),0.5f) << std::endl;
std::cout << slerp(vectorA.normalize(),vectorB.normalize(),0.5f) << std::endl;
```

--- 

[![Static Badge](https://img.shields.io/badge/Back_to_Top-gray)](#class-vector "点击返回顶部") | [![Static Badge](https://img.shields.io/badge/Back_to_README-gray)](../README.md "点击返回README")


