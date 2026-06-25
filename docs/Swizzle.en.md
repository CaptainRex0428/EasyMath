# Struct `Swizzle` 

[![Static Badge](https://img.shields.io/badge/Back_to_README-gray)](../README.md "Click to return to README") | [![Static Badge](https://img.shields.io/badge/Language-en-blue)](Swizzle.md "点击查看中文版本")

---

![Static Badge](https://img.shields.io/badge/Type-struct-green)
![Static Badge](https://img.shields.io/badge/State-completed-blue)

Swizzle is a template struct that implements a vector component access mechanism similar to GLSL/HLSL, providing flexible vector element access and reorganization functionality

| Attribute     | Value                                                              |
| ------------- | ------------------------------------------------------------------ |
| **Namespace** | **`EasyMath`**                                                     |
| **File**      | [**`include\Vector.h`**](../include/Vector.h)                      |
| **Parent**    | -                                                                  |
| **Feature**   | ![Static Badge](https://img.shields.io/badge/Feature-template-red) |

## Declaration

```C++
template <typename T, size_t dimension, int... Indices>
struct Swizzle
```

### Template

| Type        | Name        | Description                                                        |
| ----------- | ----------- | ------------------------------------------------------------------ |
| typename    | `T`         | Element type, constrained by `std::is_arithmetic_v` (arithmetic type) |
| size_t      | `dimension` | Dimension of the original vector                                   |
| int...      | `Indices`   | Variadic template, specifies the sequence of element indices to access |

## Core Features

### 1D Swizzle Specialization

| Access   | Modifier                           | Return              | Name                    | Parameter           | Description                              |
| -------- | ---------------------------------- | ------------------- | ----------------------- | ------------------- | ---------------------------------------- |
| `public` | `template<size_t D = swizzle_dim>` | `Swizzle&`          | `operator=(const T&)`   | `const T& value`    | Assign a single scalar to 1D Swizzle     |
| `public` | `template<size_t D = swizzle_dim>` | `T`                 | `operator T()`          | -                   | Implicit conversion of 1D Swizzle to scalar |

### Multi-dimensional Swizzle Functions

| Access   | Modifier                           | Return                        | Name                                  | Parameter                              | Description                        |
| -------- | ---------------------------------- | ----------------------------- | ------------------------------------- | -------------------------------------- | ---------------------------------- |
| `public` | `template<size_t D = swizzle_dim>` | `Swizzle&`                    | `operator=(const Vector<T, D>&)`     | `const Vector<T, swizzle_dim>& vec`   | Assign vector to Swizzle           |
| `public` | `template<size_t D = swizzle_dim>` | `Vector<T, swizzle_dim>`      | `operator Vector<T, D>()`            | -                                      | Convert Swizzle to vector          |
| `public` | `static` `constexpr`               | `size_t`                      | `swizzle_dim`                        | -                                      | Return the dimension of Swizzle    |
| `public` | `static` `constexpr`               | `bool`                        | `check_indices()`                    | -                                      | Compile-time check if all indices are within range |

## Supported Access Patterns

### Coordinate System Access (xyzw)

```C++
EasyMath::Vector<float, 4> vec{1.0f, 2.0f, 3.0f, 4.0f};

// Single component access
float x = vec.x;    // 1.0f
float y = vec.y;    // 2.0f
float z = vec.z;    // 3.0f
float w = vec.w;    // 4.0f

// Two-component access
EasyMath::Vector<float, 2> xy = vec.xy;   // {1.0f, 2.0f}
EasyMath::Vector<float, 2> yx = vec.yx;   // {2.0f, 1.0f}
EasyMath::Vector<float, 2> xz = vec.xz;   // {1.0f, 3.0f}

// Three-component access
EasyMath::Vector<float, 3> xyz = vec.xyz; // {1.0f, 2.0f, 3.0f}
EasyMath::Vector<float, 3> zyx = vec.zyx; // {3.0f, 2.0f, 1.0f}

// Four-component access
EasyMath::Vector<float, 4> wzyx = vec.wzyx; // {4.0f, 3.0f, 2.0f, 1.0f}
```

### Color System Access (rgba)

```C++
EasyMath::Vector<float, 4> color{0.5f, 0.7f, 0.9f, 1.0f};

// Single component access
float r = color.r;    // 0.5f
float g = color.g;    // 0.7f
float b = color.b;    // 0.9f
float a = color.a;    // 1.0f

// Multi-component access
EasyMath::Vector<float, 3> rgb = color.rgb;   // {0.5f, 0.7f, 0.9f}
EasyMath::Vector<float, 3> bgr = color.bgr;   // {0.9f, 0.7f, 0.5f}
EasyMath::Vector<float, 4> argb = color.argb; // {1.0f, 0.5f, 0.7f, 0.9f}
```

## Operator Overloading

### Swizzle-to-Swizzle Operations

| Operator | Description                                | Example                         |
| -------- | ------------------------------------------ | ------------------------------- |
| `+`      | Addition between Swizzles                  | `vec.xy + vec.zw`               |
| `-`      | Subtraction between Swizzles               | `vec.xyz - vec.zyx`             |
| `*`      | Multiplication between Swizzles (element-wise) | `vec.xy * vec.yx`               |
| `/`      | Division between Swizzles (element-wise)   | `vec.xz / vec.zx`               |

### Swizzle-to-Vector Operations

| Operator | Description                                | Example                              |
| -------- | ------------------------------------------ | ------------------------------------ |
| `+`      | Addition of Swizzle and Vector            | `vec.xy + Vector2{1.0f, 2.0f}`      |
| `-`      | Subtraction of Swizzle and Vector         | `vec.xyz - Vector3{0.5f, 0.5f, 0.5f}` |
| `*`      | Multiplication of Swizzle and Vector (element-wise) | `vec.xy * Vector2{2.0f, 3.0f}`      |
| `/`      | Division of Swizzle and Vector (element-wise) | `vec.xyz / Vector3{2.0f, 2.0f, 2.0f}` |

### Swizzle-to-Scalar Operations

| Operator | Description                                | Example                         |
| -------- | ------------------------------------------ | ------------------------------- |
| `+`      | Addition of Swizzle and scalar            | `vec.xy + 1.0f`                 |
| `-`      | Subtraction of Swizzle and scalar         | `vec.xyz - 0.5f`                |
| `*`      | Multiplication of Swizzle and scalar      | `vec.xy * 2.0f`                 |
| `/`      | Division of Swizzle and scalar            | `vec.xyz / 2.0f`                |

## Assignment Operations

```C++
EasyMath::Vector<float, 4> vec{1.0f, 2.0f, 3.0f, 4.0f};

// Single component assignment
vec.x = 5.0f;        // vec = {5.0f, 2.0f, 3.0f, 4.0f}
vec.w = 10.0f;       // vec = {5.0f, 2.0f, 3.0f, 10.0f}

// Multi-component assignment
vec.xy = EasyMath::Vector<float, 2>{7.0f, 8.0f};  // vec = {7.0f, 8.0f, 3.0f, 10.0f}
vec.zw = vec.xy;     // vec = {7.0f, 8.0f, 7.0f, 8.0f}

// Rearrangement assignment
vec.wzyx = vec.xyzw; // vec = {8.0f, 7.0f, 8.0f, 7.0f}
```

## Swizzle Application Areas

### 1. Shader-Style Vector Operations

Swizzle allows developers to manipulate vectors using syntax similar to GLSL/HLSL shader languages, making code more intuitive and concise.

```C++
// Traditional approach
EasyMath::Vector<float, 3> position{10.0f, 20.0f, 30.0f};
EasyMath::Vector<float, 2> screenPos{position[0], position[1]};

// Swizzle approach
EasyMath::Vector<float, 2> screenPos = position.xy;  // More concise and intuitive
```

### 2. Color Channel Operations

In graphics programming, Swizzle is particularly suitable for handling color data channel operations.

```C++
// Color processing example
EasyMath::Vector<float, 4> color{1.0f, 0.5f, 0.3f, 0.8f};  // RGBA

// Extract RGB components (ignore Alpha)
EasyMath::Vector<float, 3> rgb = color.rgb;

// Swap red and green channels
color.rg = color.gr;  // color = {0.5f, 1.0f, 0.3f, 0.8f}

// Create grayscale value
float gray = (color.r + color.g + color.b) / 3.0f;
color.rgb = EasyMath::Vector<float, 3>{gray, gray, gray};

// Premultiply Alpha
color.rgb = color.rgb * color.a;
```

### 3. Vector Normalization and Projection

```C++
EasyMath::Vector<float, 4> vec4{3.0f, 4.0f, 5.0f, 1.0f};

// Normalize only xyz components (keep w unchanged)
EasyMath::Vector<float, 3> xyz = vec4.xyz;
xyz = xyz.normalized();
vec4.xyz = xyz;

// 2D projection (ignore z component)
EasyMath::Vector<float, 2> projected2D = vec4.xy / vec4.z;

// Homogeneous coordinate conversion
EasyMath::Vector<float, 3> cartesian = vec4.xyz / vec4.w;
```

### 4. Fast Vector Reorganization

Swizzle provides an efficient way to reorganize vector elements, avoiding manual creation of temporary vectors.

```C++
EasyMath::Vector<float, 3> velocity{10.0f, 20.0f, 30.0f};

// Create different permutations
EasyMath::Vector<float, 3> reversed = velocity.zyx;     // {30.0f, 20.0f, 10.0f}
EasyMath::Vector<float, 3> repeated = velocity.xxx;     // {10.0f, 10.0f, 10.0f}
EasyMath::Vector<float, 3> mixed = velocity.xzy;        // {10.0f, 30.0f, 20.0f}

// Used for cross product optimization
EasyMath::Vector<float, 3> a{1.0f, 2.0f, 3.0f};
EasyMath::Vector<float, 3> b{4.0f, 5.0f, 6.0f};
EasyMath::Vector<float, 3> cross_result = a.yzx * b.zxy - a.zxy * b.yzx;
```

### 5. Texture Coordinate Operations

In graphics rendering, Swizzle simplifies texture coordinate handling.

```C++
// UV coordinate processing
EasyMath::Vector<float, 4> texCoord{0.5f, 0.7f, 0.0f, 1.0f};

// Extract UV components
EasyMath::Vector<float, 2> uv = texCoord.xy;

// UV flip
texCoord.xy = EasyMath::Vector<float, 2>{texCoord.x, 1.0f - texCoord.y};

// UV scaling
texCoord.xy = texCoord.xy * 2.0f;

// Create cube map coordinates
EasyMath::Vector<float, 3> cubeMapCoord = texCoord.xyz;
```

### 6. Applications in Physics Simulation

```C++
// Particle system
struct Particle {
    EasyMath::Vector<float, 4> positionAndMass;  // xyz=position, w=mass
    EasyMath::Vector<float, 4> velocityAndLife;  // xyz=velocity, w=lifetime
};

Particle particle;
particle.positionAndMass = {10.0f, 20.0f, 30.0f, 1.5f};
particle.velocityAndLife = {1.0f, 2.0f, 3.0f, 100.0f};

// Update position (only modify xyz, keep mass unchanged)
particle.positionAndMass.xyz = particle.positionAndMass.xyz + 
                               particle.velocityAndLife.xyz * deltaTime;

// Apply gravity (only affects y-direction velocity)
particle.velocityAndLife.y -= 9.8f * deltaTime;

// Get kinetic energy
float kineticEnergy = 0.5f * particle.positionAndMass.w * 
                     particle.velocityAndLife.xyz.lengthSquared();
```

### 7. Matrix Column/Row Access Optimization

When used with the Matrix class, Swizzle can simplify matrix column and row operations.

```C++
// Extract information from transformation matrix
EasyMath::Matrix<float, 4, 4> transform = GetTransformMatrix();
EasyMath::Vector<float, 4> column0 = GetColumn(transform, 0);

// Extract scale factor
float scaleX = column0.xyz.length();

// Extract direction vector (normalized first three components)
EasyMath::Vector<float, 3> forward = column0.xyz.normalized();

// Construct new transformation
EasyMath::Vector<float, 4> newColumn;
newColumn.xyz = forward * scaleX;
newColumn.w = 0.0f;  // w component of direction vector is 0
```

### 8. Data Packing and Unpacking

Swizzle is very useful in data packing and unpacking scenarios, especially when communicating with the GPU.

```C++
// Pack normal and tangent information
struct CompressedVertex {
    EasyMath::Vector<float, 4> normalAndTangentSign;  // xyz=normal, w=tangent sign
};

CompressedVertex vertex;
EasyMath::Vector<float, 3> normal{0.0f, 1.0f, 0.0f};
float tangentSign = 1.0f;

// Pack
vertex.normalAndTangentSign.xyz = normal;
vertex.normalAndTangentSign.w = tangentSign;

// Unpack
EasyMath::Vector<float, 3> unpackedNormal = vertex.normalAndTangentSign.xyz;
float unpackedSign = vertex.normalAndTangentSign.w;

// Reconstruct tangent
EasyMath::Vector<float, 3> tangent = CalculateTangent(unpackedNormal) * unpackedSign;
```

## Performance Optimization Recommendations

1. **Compile-time optimization**: All index checks of Swizzle are completed at compile time, with no runtime overhead
2. **Inline optimization**: Swizzle operations are typically inlined by the compiler, with performance comparable to direct array access
3. **Avoid overuse**: Although Swizzle syntax is concise, excessive chaining may affect code readability
4. **Cache-friendly**: When using Swizzle for batch operations, pay attention to data locality

```C++
// Recommended: Process multiple components at once
vec.xyz = vec.xyz * 2.0f;

// Not recommended: Multiple individual accesses
vec.x *= 2.0f;
vec.y *= 2.0f;
vec.z *= 2.0f;
```

## Notes

1. **Index out of bounds**: Accessing components beyond the vector dimension will result in a compile-time error
```C++
EasyMath::Vector<float, 2> vec2{1.0f, 2.0f};
// float z = vec2.z;  // Compile error
```

2. **Duplicate index assignment**: Certain Swizzle patterns do not support assignment operations
```C++
EasyMath::Vector<float, 3> vec{1.0f, 2.0f, 3.0f};
// vec.xxx = Vector3{4.0f, 5.0f, 6.0f}; // Compile error
```

3. **Type compatibility**: The dimension of the Swizzle result must match the target type
```C++
EasyMath::Vector<float, 4> vec4{1.0f, 2.0f, 3.0f, 4.0f};
// EasyMath::Vector<float, 3> vec3 = vec4.xy; // Compile error
EasyMath::Vector<float, 2> vec2 = vec4.xy;
```


--- 

[![Static Badge](https://img.shields.io/badge/Back_to_Top-gray)](#struct-swizzle "Click to return to top") | [![Static Badge](https://img.shields.io/badge/Back_to_README-gray)](../README.md "Click to return to README")
