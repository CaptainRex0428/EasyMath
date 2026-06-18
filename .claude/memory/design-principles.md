# 设计原则与核心机制

## 核心设计哲学

EasyMath 的设计围绕三个核心原则：
1. **编译期类型安全** - 用类型系统防止错误
2. **零运行时开销** - 模板元编程 + constexpr
3. **直观的 API** - 图形程序员熟悉的接口

---

## 1. 泛型设计模式

### Vector - 维度泛型

```cpp
template<typename T, size_t dimension, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
class Vector
```

**设计要点**：
- `T` - 元素类型（float, double, int 等）
- `size_t dimension` - 编译期常量维度
- `std::enable_if_t<std::is_arithmetic_v<T>>` - SFINAE 类型约束

**底层存储**：
```cpp
std::array<T, dimension> data;
```
- ✅ 编译期已知大小 → 零运行时开销
- ✅ 栈上分配 → 高速缓存友好
- ✅ 类型安全 → 编译器验证大小

**支持的维度**：
```cpp
Vector<float, 2>    // 2D 向量
Vector<float, 3>    // 3D 向量
Vector<float, 4>    // 4D 向量（齐次坐标）
Vector<float, 100>  // 甚至任意维度！
```

### Matrix - 行列泛型

```cpp
template <typename T, size_t rows, size_t cols,
    typename = std::enable_if_t<std::is_arithmetic_v<T>>>
class Matrix
```

**存储布局 - 行主序（Row-major）**：
```cpp
std::array<T, rows * cols> data;

// 访问公式
data[row * cols + col]  // 第 row 行，第 col 列
```

**为什么选择行主序？**
- ✅ C/C++ 多维数组的自然顺序
- ✅ 与 DirectX 一致
- ⚠️ 注意：OpenGL 使用列主序

**示例**：
```
3×3 矩阵内存布局：
逻辑视图:        内存布局:
[0, 1, 2]       [0, 1, 2, 3, 4, 5, 6, 7, 8]
[3, 4, 5]
[6, 7, 8]
```

---

## 2. Swizzle 机制 - Union 魔法

### 实现原理

```cpp
union {
    std::array<T, dimension> data;  // 实际存储
    Swizzle<T, dimension, 0> x;     // 索引 0
    Swizzle<T, dimension, 1> y;     // 索引 1
    Swizzle<T, dimension, 0, 1> xy; // 索引 0, 1
    Swizzle<T, dimension, 1, 0> yx; // 索引 1, 0（反向）
    // ... 更多组合
}
```

**Union 共享内存**：
```
内存布局（共享）：
[data[0]] [data[1]] [data[2]] [data[3]]
    ↓         ↓         ↓         ↓
    x         y         z         w     (几何坐标)
    r         g         b         a     (颜色通道)
```

**Swizzle 模板参数 = 索引序列**：
```cpp
Swizzle<T, dimension, 0, 1, 2>
//                    ↑  ↑  ↑
//                    访问索引 [0], [1], [2]

vec.xyz  → 访问 data[0], data[1], data[2]
vec.zyx  → 访问 data[2], data[1], data[0]  // 反向！
vec.xxyy → 访问 data[0], data[0], data[1], data[1]  // 重复！
```

**零开销抽象**：
- 编译期完成所有索引计算
- 运行时等价于直接内存访问
- 类型安全（编译期检查）

---

## 3. Color 类型安全系统

### 设计核心：模板参数编码类型状态

```cpp
template <typename T, bool bIsRGB, size_t modeIndex, typename>
class Color : public Vector<T, 4>
```

**两个关键参数**：
- `bool bIsRGB` - 是否是 RGB 家族（vs HSV/HSL/HSI 家族）
- `size_t modeIndex` - 具体模式索引

### 类型别名系统

```cpp
// RGB 家族 (bIsRGB = true)
using sRGBColor<T>   = Color<T, true, 1>;  // 显示颜色空间
using LinearColor<T> = Color<T, true, 0>;  // 物理光强度

// HSx 家族 (bIsRGB = false)
using HSV<T> = Color<T, false, 0>;  // 色相-饱和度-明度
using HSL<T> = Color<T, false, 1>;  // 色相-饱和度-亮度
using HSI<T> = Color<T, false, 2>;  // 色相-饱和度-强度
```

### 类型安全的转换系统

**关键设计**：用 SFINAE 控制方法可见性

```cpp
// 只有 RGB 家族才能调用 toHSV
template <bool rgbmode = bIsRGB>
typename std::enable_if_t<rgbmode, HSV<T>>
toHSV() const

// 只有 sRGB 才能调用 toLinear
template <bool rgbmode = bIsRGB>
typename std::enable_if_t<rgbmode && modeIndex == 1, LinearColor<T>>
toLinear() const
```

**调用权限矩阵**：

| 类型 | bIsRGB | modeIndex | toHSV() | toLinear() | toSRGB() |
|------|--------|-----------|---------|------------|----------|
| sRGBColor | true | 1 | ✅ | ✅ | ❌ |
| LinearColor | true | 0 | ✅ | ❌ | ✅ |
| HSV | false | 0 | ❌ | ❌ | ❌ |
| HSL | false | 1 | ❌ | ❌ | ❌ |
| HSI | false | 2 | ❌ | ❌ | ❌ |

**防止的错误**：
```cpp
// ❌ 防止：HSV 再转 HSV（无意义）
HSV<float> hsv;
// hsv.toHSV();  // 编译错误！

// ❌ 防止：LinearColor 转 Linear（已经是 Linear）
LinearColor<float> linear;
// linear.toLinear();  // 编译错误！

// ❌ 防止：sRGB 转 sRGB（已经是 sRGB）
sRGBColor<float> srgb;
// srgb.toSRGB();  // 编译错误！
```

**强制的正确转换路径**：
```cpp
// ✅ RGB → HSx（实例方法）
sRGBColor<float> srgb;
HSV<float> hsv = srgb.toHSV();

// ✅ HSx → RGB（静态方法）
sRGBColor<float> back = sRGBColor<float>::fromHSV(hsv);

// ✅ sRGB ↔ Linear（往返）
LinearColor<float> linear = srgb.toLinear();
sRGBColor<float> back = linear.toSRGB();
```

---

## 4. 颜色空间语义

### sRGB vs Linear - 何时使用

**sRGB（标准 RGB）**：
```cpp
sRGBColor<float> displayColor{0.5f, 0.3f, 0.8f, 1.0f};
```
- 📺 **用途**：屏幕显示、UI、纹理采样
- 🎨 **特性**：伽马校正（γ ≈ 2.2）、人眼感知均匀
- 📊 **值范围**：[0, 1] 对应显示器的像素值

**Linear RGB**：
```cpp
LinearColor<float> physicalLight{0.214f, 0.073f, 0.604f, 1.0f};
```
- 💡 **用途**：光照计算、HDR、物理正确渲染
- ⚡ **特性**：物理线性、没有伽马校正
- 🔬 **值范围**：[0, ∞) 对应实际光强度（支持 HDR）

**关键规则**：
```
显示器显示 → sRGB
光照计算   → Linear
颜色混合   → Linear（物理正确）
后处理     → 根据效果选择
```

### 伽马校正公式

**sRGB → Linear**（解码）：
```cpp
if (srgb <= 0.04045)
    linear = srgb / 12.92           // 线性段
else
    linear = pow((srgb + 0.055) / 1.055, 2.4)  // 幂函数段
```

**Linear → sRGB**（编码）：
```cpp
if (linear <= 0.0031308)
    srgb = linear * 12.92           // 线性段
else
    srgb = 1.055 * pow(linear, 1/2.4) - 0.055  // 幂函数段
```

### HSV/HSL/HSI - 三种色彩模型

**HSV（Hue, Saturation, Value）**：
- **V = max(R, G, B)** - 明度 = 最大分量
- 🎨 艺术家最常用（Photoshop、Blender）
- 适合调整颜色亮度和鲜艳度

**HSL（Hue, Saturation, Lightness）**：
- **L = (max + min) / 2** - 亮度 = 中点
- 🌐 Web 开发常用（CSS `hsl()`）
- 亮度更"感知均匀"

**HSI（Hue, Saturation, Intensity）**：
- **I = (R + G + B) / 3** - 强度 = 算术平均
- 🔬 图像处理/计算机视觉
- 数学上最简单

### RGB → HSV/HSL/HSI 的语义

**sRGB → HSV**：有意义，艺术家直观
```cpp
sRGBColor<float> purple{0.5f, 0.3f, 0.8f};
HSV<float> hsv = purple.toHSV();
// H ≈ 270° (紫色调)
// S ≈ 0.625 (中等饱和)
// V ≈ 0.8 (较亮)
// ✅ 艺术家可理解："一个亮紫色"
```

**Linear → HSV**：数学可行，语义混乱
```cpp
LinearColor<float> light{0.214f, 0.073f, 0.604f};
// 如果允许：hsv = light.toHSV();
// ⚠️ 这个 HSV 代表"光强度的 HSV 表示"
// ❌ 不是感知上的色相/饱和度
// 建议：先转 sRGB 再转 HSV
```

---

## 5. 亮度计算标准

### 三种国际标准

```cpp
enum class LuminanceStandard {
    Rec601,   // NTSC 标准
    Rec709,   // HDTV 标准（最常用）
    Rec2020   // UHD 4K/8K 标准
};
```

**公式**：
```cpp
Rec601:  L = 0.299*R  + 0.587*G  + 0.114*B   // 老电视
Rec709:  L = 0.2126*R + 0.7152*G + 0.0722*B  // 高清电视 ⭐
Rec2020: L = 0.2627*R + 0.6780*G + 0.0593*B  // 超高清
```

**为什么绿色权重最高？**
- 人眼对绿色（550nm）最敏感
- 视网膜绿色感光细胞最多
- 绿色贡献约 70% 的感知亮度

### 三种亮度类型

**1. 物理正确亮度（Luminance）**：
```cpp
T luminance(LuminanceStandard ls = Rec709)
```
- 在 **Linear 空间**计算
- sRGB 会先转 Linear 再计算
- 用于：光照计算、HDR

**2. 感知亮度（Perceived Luminance）**：
```cpp
T perceivedLuminance(LuminanceStandard ls = Rec709) const
```
- 直接在当前颜色空间计算
- 不转换到 Linear
- 用于：UI 显示、视觉效果

**3. 相对亮度（Relative Luminance - WCAG）**：
```cpp
T relativeLuminance() const
```
- 始终用 Rec709 标准
- 始终在 Linear 空间
- 用于：对比度计算、无障碍访问

**对比度计算（WCAG 标准）**：
```cpp
T contrastRatio(const Color& other) const
{
    T L1 = relativeLuminance();
    T L2 = other.relativeLuminance();
    if (L2 > L1) std::swap(L1, L2);
    return (L1 + 0.05) / (L2 + 0.05);  // 范围 [1, 21]
}

// WCAG 标准：
// AA 级别：对比度 ≥ 4.5:1（正文）
// AAA 级别：对比度 ≥ 7:1（正文）
```

---

## 6. Matrix 线性代数设计

### 行列式计算 - 递归算法

```cpp
template<size_t R = rows, size_t C = cols>
typename std::enable_if_t<R == C, T> determinant() const
{
    if constexpr (rows == 1) {
        return data[0];  // 基础情况
    }
    else if constexpr (rows == 2) {
        return data[0] * data[3] - data[1] * data[2];  // 优化 2×2
    }
    else {
        // 余子式展开（第一行）
        T det = 0;
        for (size_t i = 0; i < cols; ++i) {
            auto submat = submatrix(0, i);
            det += (i % 2 == 0 ? 1 : -1) * data[i] * submat.determinant();
        }
        return det;
    }
}
```

**特点**：
- `if constexpr` - 编译期分支，零运行时开销
- 2×2 特化 - 避免不必要的递归
- 递归深度 = 矩阵大小（3×3 递归3层，4×4 递归4层）

### 逆矩阵 - 伴随矩阵法

```
数学公式：
A⁻¹ = adj(A) / det(A)

其中：
- det(A) = 行列式
- adj(A) = 伴随矩阵 = cofactor(A)ᵀ
```

**实现步骤**：
```cpp
Matrix inverse() const
{
    T det = determinant();  // 1. 计算行列式
    assert(det != 0);       // 2. 检查可逆性

    auto adj = adjugate();  // 3. 计算伴随矩阵
    return (1.0 / det) * adj;  // 4. 缩放
}
```

---

## 7. 类型安全最佳实践

### 使用 SFINAE 防止错误

**模板参数约束**：
```cpp
// 只允许算术类型
template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
class Vector;

// 只允许方阵求行列式
template<size_t R = rows, size_t C = cols>
typename std::enable_if_t<R == C, T> determinant() const;

// 只允许 RGB 家族转 HSV
template <bool rgbmode = bIsRGB>
typename std::enable_if_t<rgbmode, HSV<T>>
toHSV() const;
```

### 使用 static_assert 提供友好错误

```cpp
static_assert(std::is_arithmetic_v<T>,
    "Vector element type must be arithmetic");

static_assert(dimension > 0,
    "vector dimension must be positive.");

static_assert(rows > 0 && cols > 0,
    "Matrix row & column must be positive.");
```

### 使用 constexpr 进行编译期计算

```cpp
static constexpr size_t Size() noexcept { return dimension; }
static constexpr size_t Rows() { return rows; }
static constexpr size_t Cols() { return cols; }

virtual constexpr T lengthSquared() const noexcept
{
    T sum = T(0);
    for (size_t idx = 0; idx < dimension; ++idx) {
        sum += data[idx] * data[idx];
    }
    return sum;
}
```

---

## 8. 继承策略

### Color 继承 Vector - 合理的 IS-A 关系

```cpp
class Color : public Vector<T, 4>
```

**优势**：
- ✅ 复用 Vector 的所有功能（算术、Swizzle、迭代器）
- ✅ Color 确实"是一个"4D 向量（RGBA）
- ✅ 多态行为（可以当 Vector 使用）

**注意事项**：
- 析构函数为 `virtual` - 支持多态删除
- 重写 `print()` - 自定义输出格式

### Quaternion 继承 Vector - 复用实现

```cpp
class Quaternion : public Vector<T, 4>
```

**设计考量**：
- ✅ 复用存储和基本操作
- ⚠️ 语义上 Quaternion 不是"向量"（是旋转表示）
- ✅ 重写关键方法（`normalized()`, `print()`）

---

## 9. 命名约定总结

### 类型名
```cpp
Vector<T, N>        // PascalCase，泛型参数
Matrix<T, R, C>
QuaternionF         // 类型别名用完整单词
sRGBColor<T>        // 首字母小写的缩写
```

### 函数名
```cpp
// 成员函数：camelCase
normalized()
toLinear()
lengthSquared()

// 全局函数：camelCase
dot(a, b)
cross(a, b)
distance(a, b)

// 矩阵构造函数：MTX 前缀
MTXIdentity()
MTXRotationX()
MTXTranslation()
```

### 常量
```cpp
// 数学常量：UPPER_SNAKE_CASE
PI
E
NEARZERO_THRESHOLD

// 枚举值：PascalCase
enum class LuminanceStandard {
    Rec601,
    Rec709,
    Rec2020
};
```

---

## 10. 性能优化策略

### 编译期优化
- `constexpr` - 编译期常量折叠
- `if constexpr` - 编译期分支消除
- 模板实例化 - 为每个具体类型生成专门代码

### 内存布局优化
- `std::array` - 连续内存，缓存友好
- 栈分配 - 避免堆分配开销
- 行主序 - 符合 C++ 内存模型

### 零抽象开销
- Swizzle - 编译期索引计算
- 类型别名 - 零运行时成本
- SFINAE - 编译期方法过滤

---

## 设计哲学总结

EasyMath 的核心设计理念：

1. **让编译器工作** - 用类型系统在编译期捕获错误
2. **零运行时开销** - 所有抽象在编译期展开
3. **物理正确优先** - 遵循图形学标准（sRGB、WCAG、Rec709）
4. **API 自文档化** - 类型名即文档，错误用法无法编译
5. **继承有度** - 只在真正 IS-A 关系时使用继承
