# `Swizzle`结构体

[![Static Badge](https://img.shields.io/badge/Back_to_README-gray)](../README.md "点击返回README") | [![Static Badge](https://img.shields.io/badge/Language-zh-red)](Swizzle.en.md "click to english version")

---

![Static Badge](https://img.shields.io/badge/Type-struct-green)
![Static Badge](https://img.shields.io/badge/State-completed-blue)

Swizzle是一个模板结构体，实现了类似GLSL/HLSL中的向量分量访问机制，提供灵活的向量元素访问和重组功能

| Attribute     | Value                                                              |
| ------------- | ------------------------------------------------------------------ |
| **Namespace** | **`EM`**                                                           |
| **File**      | [**`include\Vector.h`**](../include/Vector.h)                      |
| **Parent**    | -                                                                  |
| **Feature**   | ![Static Badge](https://img.shields.io/badge/Feature-template-red) |

## 声明

```C++
template <typename T, size_t dimension, int... Indices>
struct Swizzle
```

### 模板

| Type        | Name        | Description                                               |
| ----------- | ----------- | --------------------------------------------------------- |
| typename    | `T`         | 元素类型，受到 `std::is_arithmetic_v` （算数类型）限制    |
| size_t      | `dimension` | 原始向量的维度                                           |
| int...      | `Indices`   | 可变参数模板，指定要访问的元素索引序列                    |

## 核心特性

### 1D Swizzle特化

| Access   | Modifier                           | Return              | Name                    | Parameter           | Description                |
| -------- | ---------------------------------- | ------------------- | ----------------------- | ------------------- | -------------------------- |
| `public` | `template<size_t D = swizzle_dim>` | `Swizzle&`          | `operator=(const T&)`   | `const T& value`    | 为1D Swizzle赋值单个标量   |
| `public` | `template<size_t D = swizzle_dim>` | `T`                 | `operator T()`          | -                   | 1D Swizzle隐式转换为标量   |

### 多维Swizzle功能

| Access   | Modifier                           | Return                        | Name                                  | Parameter                              | Description                        |
| -------- | ---------------------------------- | ----------------------------- | ------------------------------------- | -------------------------------------- | ---------------------------------- |
| `public` | `template<size_t D = swizzle_dim>` | `Swizzle&`                    | `operator=(const Vector<T, D>&)`     | `const Vector<T, swizzle_dim>& vec`   | 从向量赋值给Swizzle                |
| `public` | `template<size_t D = swizzle_dim>` | `Vector<T, swizzle_dim>`      | `operator Vector<T, D>()`            | -                                      | Swizzle转换为向量                  |
| `public` | `static` `constexpr`               | `size_t`                      | `swizzle_dim`                        | -                                      | 返回Swizzle的维度                  |
| `public` | `static` `constexpr`               | `bool`                        | `check_indices()`                    | -                                      | 编译期检查所有索引是否在范围内      |

## 支持的访问模式

### 坐标系访问 (xyzw)

```C++
EM::Vector<float, 4> vec{1.0f, 2.0f, 3.0f, 4.0f};

// 单分量访问
float x = vec.x;    // 1.0f
float y = vec.y;    // 2.0f
float z = vec.z;    // 3.0f
float w = vec.w;    // 4.0f

// 双分量访问
EM::Vector<float, 2> xy = vec.xy;   // {1.0f, 2.0f}
EM::Vector<float, 2> yx = vec.yx;   // {2.0f, 1.0f}
EM::Vector<float, 2> xz = vec.xz;   // {1.0f, 3.0f}

// 三分量访问
EM::Vector<float, 3> xyz = vec.xyz; // {1.0f, 2.0f, 3.0f}
EM::Vector<float, 3> zyx = vec.zyx; // {3.0f, 2.0f, 1.0f}

// 四分量访问
EM::Vector<float, 4> wzyx = vec.wzyx; // {4.0f, 3.0f, 2.0f, 1.0f}
```

### 颜色系访问 (rgba)

```C++
EM::Vector<float, 4> color{0.5f, 0.7f, 0.9f, 1.0f};

// 单分量访问
float r = color.r;    // 0.5f
float g = color.g;    // 0.7f
float b = color.b;    // 0.9f
float a = color.a;    // 1.0f

// 多分量访问
EM::Vector<float, 3> rgb = color.rgb;   // {0.5f, 0.7f, 0.9f}
EM::Vector<float, 3> bgr = color.bgr;   // {0.9f, 0.7f, 0.5f}
EM::Vector<float, 4> argb = color.argb; // {1.0f, 0.5f, 0.7f, 0.9f}
```

## 运算符重载

### Swizzle与Swizzle运算

| Operator | Description                          | Example                         |
| -------- | ------------------------------------ | ------------------------------- |
| `+`      | Swizzle之间的加法                    | `vec.xy + vec.zw`               |
| `-`      | Swizzle之间的减法                    | `vec.xyz - vec.zyx`             |
| `*`      | Swizzle之间的乘法（逐元素）           | `vec.xy * vec.yx`               |
| `/`      | Swizzle之间的除法（逐元素）           | `vec.xz / vec.zx`               |

### Swizzle与Vector运算

| Operator | Description                          | Example                              |
| -------- | ------------------------------------ | ------------------------------------ |
| `+`      | Swizzle与Vector加法                  | `vec.xy + Vector2{1.0f, 2.0f}`      |
| `-`      | Swizzle与Vector减法                  | `vec.xyz - Vector3{0.5f, 0.5f, 0.5f}` |
| `*`      | Swizzle与Vector乘法（逐元素）         | `vec.xy * Vector2{2.0f, 3.0f}`      |
| `/`      | Swizzle与Vector除法（逐元素）         | `vec.xyz / Vector3{2.0f, 2.0f, 2.0f}` |

### Swizzle与标量运算

| Operator | Description                          | Example                         |
| -------- | ------------------------------------ | ------------------------------- |
| `+`      | Swizzle与标量加法                    | `vec.xy + 1.0f`                 |
| `-`      | Swizzle与标量减法                    | `vec.xyz - 0.5f`                |
| `*`      | Swizzle与标量乘法                    | `vec.xy * 2.0f`                 |
| `/`      | Swizzle与标量除法                    | `vec.xyz / 2.0f`                |

## 赋值操作

```C++
EM::Vector<float, 4> vec{1.0f, 2.0f, 3.0f, 4.0f};

// 单分量赋值
vec.x = 5.0f;        // vec = {5.0f, 2.0f, 3.0f, 4.0f}
vec.w = 10.0f;       // vec = {5.0f, 2.0f, 3.0f, 10.0f}

// 多分量赋值
vec.xy = EM::Vector<float, 2>{7.0f, 8.0f};  // vec = {7.0f, 8.0f, 3.0f, 10.0f}
vec.zw = vec.xy;     // vec = {7.0f, 8.0f, 7.0f, 8.0f}

// 重排赋值
vec.wzyx = vec.xyzw; // vec = {8.0f, 7.0f, 8.0f, 7.0f}
```

## Swizzle应用板块

### 1. 着色器风格的向量操作

Swizzle允许开发者使用类似GLSL/HLSL着色器语言的语法来操作向量，使代码更加直观和简洁。

```C++
// 传统方式
EM::Vector<float, 3> position{10.0f, 20.0f, 30.0f};
EM::Vector<float, 2> screenPos{position[0], position[1]};

// Swizzle方式
EM::Vector<float, 2> screenPos = position.xy;  // 更简洁直观
```

### 2. 颜色通道操作

在图形编程中，Swizzle特别适合处理颜色数据的通道操作。

```C++
// 颜色处理示例
EM::Vector<float, 4> color{1.0f, 0.5f, 0.3f, 0.8f};  // RGBA

// 提取RGB分量（忽略Alpha）
EM::Vector<float, 3> rgb = color.rgb;

// 交换红绿通道
color.rg = color.gr;  // color = {0.5f, 1.0f, 0.3f, 0.8f}

// 创建灰度值
float gray = (color.r + color.g + color.b) / 3.0f;
color.rgb = EM::Vector<float, 3>{gray, gray, gray};

// 预乘Alpha
color.rgb = color.rgb * color.a;
```

### 3. 向量归一化和投影

```C++
EM::Vector<float, 4> vec4{3.0f, 4.0f, 5.0f, 1.0f};

// 只归一化xyz分量（保持w不变）
EM::Vector<float, 3> xyz = vec4.xyz;
xyz = xyz.normalized();
vec4.xyz = xyz;

// 2D投影（忽略z分量）
EM::Vector<float, 2> projected2D = vec4.xy / vec4.z;

// 齐次坐标转换
EM::Vector<float, 3> cartesian = vec4.xyz / vec4.w;
```

### 4. 快速向量重组

Swizzle提供了一种高效的向量元素重组方式，避免了手动创建临时向量。

```C++
EM::Vector<float, 3> velocity{10.0f, 20.0f, 30.0f};

// 创建不同的排列组合
EM::Vector<float, 3> reversed = velocity.zyx;     // {30.0f, 20.0f, 10.0f}
EM::Vector<float, 3> repeated = velocity.xxx;     // {10.0f, 10.0f, 10.0f}
EM::Vector<float, 3> mixed = velocity.xzy;        // {10.0f, 30.0f, 20.0f}

// 用于叉乘优化
EM::Vector<float, 3> a{1.0f, 2.0f, 3.0f};
EM::Vector<float, 3> b{4.0f, 5.0f, 6.0f};
EM::Vector<float, 3> cross_result = a.yzx * b.zxy - a.zxy * b.yzx;
```

### 5. 纹理坐标操作

在图形渲染中，Swizzle简化了纹理坐标的处理。

```C++
// UV坐标处理
EM::Vector<float, 4> texCoord{0.5f, 0.7f, 0.0f, 1.0f};

// 提取UV分量
EM::Vector<float, 2> uv = texCoord.xy;

// UV翻转
texCoord.xy = EM::Vector<float, 2>{texCoord.x, 1.0f - texCoord.y};

// UV缩放
texCoord.xy = texCoord.xy * 2.0f;

// 创建立方体贴图坐标
EM::Vector<float, 3> cubeMapCoord = texCoord.xyz;
```

### 6. 物理模拟中的应用

```C++
// 粒子系统
struct Particle {
    EM::Vector<float, 4> positionAndMass;  // xyz=位置, w=质量
    EM::Vector<float, 4> velocityAndLife;  // xyz=速度, w=生命值
};

Particle particle;
particle.positionAndMass = {10.0f, 20.0f, 30.0f, 1.5f};
particle.velocityAndLife = {1.0f, 2.0f, 3.0f, 100.0f};

// 更新位置（只修改xyz，保持质量不变）
particle.positionAndMass.xyz = particle.positionAndMass.xyz + 
                               particle.velocityAndLife.xyz * deltaTime;

// 应用重力（只影响y方向速度）
particle.velocityAndLife.y -= 9.8f * deltaTime;

// 获取动能
float kineticEnergy = 0.5f * particle.positionAndMass.w * 
                     particle.velocityAndLife.xyz.lengthSquared();
```

### 7. 矩阵列/行访问优化

与Matrix类配合使用时，Swizzle可以简化矩阵的列和行操作。

```C++
// 从变换矩阵中提取信息
EM::Matrix<float, 4, 4> transform = GetTransformMatrix();
EM::Vector<float, 4> column0 = GetColumn(transform, 0);

// 提取缩放因子
float scaleX = column0.xyz.length();

// 提取方向向量（归一化后的前三个分量）
EM::Vector<float, 3> forward = column0.xyz.normalized();

// 构建新的变换
EM::Vector<float, 4> newColumn;
newColumn.xyz = forward * scaleX;
newColumn.w = 0.0f;  // 方向向量的w分量为0
```

### 8. 数据打包与解包

Swizzle在数据打包和解包场景中非常有用，特别是在与GPU通信时。

```C++
// 打包法线和切线信息
struct CompressedVertex {
    EM::Vector<float, 4> normalAndTangentSign;  // xyz=法线, w=切线符号
};

CompressedVertex vertex;
EM::Vector<float, 3> normal{0.0f, 1.0f, 0.0f};
float tangentSign = 1.0f;

// 打包
vertex.normalAndTangentSign.xyz = normal;
vertex.normalAndTangentSign.w = tangentSign;

// 解包
EM::Vector<float, 3> unpackedNormal = vertex.normalAndTangentSign.xyz;
float unpackedSign = vertex.normalAndTangentSign.w;

// 重建切线
EM::Vector<float, 3> tangent = CalculateTangent(unpackedNormal) * unpackedSign;
```

## 性能优化建议

1. **编译期优化**：Swizzle的所有索引检查都在编译期完成，运行时无额外开销
2. **内联优化**：Swizzle操作通常被编译器内联，性能与直接数组访问相当
3. **避免过度使用**：虽然Swizzle语法简洁，但过度链式调用可能影响代码可读性
4. **缓存友好**：使用Swizzle进行批量操作时，注意数据局部性

```C++
// 推荐：一次性处理多个分量
vec.xyz = vec.xyz * 2.0f;

// 不推荐：多次单独访问
vec.x *= 2.0f;
vec.y *= 2.0f;
vec.z *= 2.0f;
```

## 注意事项

1. **索引越界**：访问超出向量维度的分量会在编译期报错
```C++
EM::Vector<float, 2> vec2{1.0f, 2.0f};
// float z = vec2.z;  // 编译错误
```

2. **重复索引赋值**：某些Swizzle模式不支持赋值操作
```C++
EM::Vector<float, 3> vec{1.0f, 2.0f, 3.0f};
// vec.xxx = Vector3{4.0f, 5.0f, 6.0f}; // 编译错误
```

3. **类型兼容性**：Swizzle结果的维度必须与目标类型匹配
```C++
EM::Vector<float, 4> vec4{1.0f, 2.0f, 3.0f, 4.0f};
// EM::Vector<float, 3> vec3 = vec4.xy; // 编译错误
EM::Vector<float, 2> vec2 = vec4.xy;
```


--- 

[![Static Badge](https://img.shields.io/badge/Back_to_Top-gray)](#swizzle结构体 "点击返回顶部") | [![Static Badge](https://img.shields.io/badge/Back_to_README-gray)](../README.md "点击返回README")
