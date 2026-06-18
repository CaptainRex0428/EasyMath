# 代码规范

## 命名约定

### 类型名称
- **类/结构体**: PascalCase
  ```cpp
  class Vector
  class Matrix
  class Quaternion
  struct LuminanceStandard
  ```

- **模板参数**: PascalCase 或单字母
  ```cpp
  template<typename T, size_t N>  // 推荐
  template<typename T, size_t dimension>  // 描述性名称
  ```

### 函数名称
- **成员函数**: camelCase
  ```cpp
  void normalize()
  Vector normalized() const
  T length() const
  Matrix transpose() const
  ```

- **全局/静态函数**: camelCase
  ```cpp
  T dot(const Vector& a, const Vector& b)
  Vector cross(const Vector& a, const Vector& b)
  T distance(const Vector& a, const Vector& b)
  ```

- **矩阵构造函数**: MTX 前缀 + PascalCase
  ```cpp
  Matrix<T,4,4> MTXIdentity()
  Matrix<T,4,4> MTXRotationX(T radians)
  Matrix<T,4,4> MTXTranslation(T x, T y, T z)
  Matrix<T,4,4> MTXScale(T x, T y, T z)
  ```

### 常量
- **数学常量**: UPPER_SNAKE_CASE
  ```cpp
  constexpr double PI = 3.14159265358979323846;
  constexpr double E = 2.71828182845904523536;
  constexpr double NEARZERO_THRESHOLD = 1e-6;
  ```

- **枚举值**: PascalCase
  ```cpp
  enum class LuminanceStandard {
      Rec601,
      Rec709,
      Rec2020
  };
  ```

### 变量名称
- **成员变量**: camelCase（无前缀）
  ```cpp
  std::array<T, dimension> data;
  ```

- **局部变量**: camelCase
  ```cpp
  T result = 0;
  size_t rowIndex = 0;
  ```

- **模板参数变量**: camelCase 或描述性
  ```cpp
  template<size_t rows, size_t cols>  // 推荐
  template<size_t R, size_t C>        // 简洁版本
  ```

---

## 代码风格

### 缩进与格式
```cpp
// 使用 Tab 缩进（项目当前风格）
class Vector
{
public:
    Vector()
        : data({})  // 初始化列表对齐
    {
    }

    void function()
    {
        if (condition)
        {
            // 代码
        }
    }
};
```

### 括号风格
- **类/函数**: K&R 风格（左括号换行）
```cpp
class MyClass
{
    void function()
    {
        // 代码
    }
};
```

- **控制流**: 同行
```cpp
if (condition) {
    // 代码
} else {
    // 代码
}
```

### 空格约定
```cpp
// 运算符周围加空格
result = a + b * c;
vec = vec1 + vec2;

// 模板尖括号内不加空格
Vector<float, 3> vec;

// 函数调用括号不加空格
dot(a, b);

// 逗号后加空格
Vector{1.0f, 2.0f, 3.0f};
```

---

## 设计原则

### 1. 类型安全
```cpp
// 使用 SFINAE 限制模板参数
template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
class Vector { };

// 使用 static_assert 提供友好错误
static_assert(std::is_arithmetic_v<T>, "Vector element type must be arithmetic");
static_assert(dimension > 0, "vector dimension must be positive.");
```

### 2. Header-Only 优先
```cpp
// 大部分实现直接在头文件中
template<typename T, size_t N>
class Vector
{
public:
    T length() const noexcept
    {
        return std::sqrt(lengthSquared());
    }
};

// 仅在必要时使用 .cpp 实现（如 ByteSize）
```

### 3. constexpr 和 noexcept
```cpp
// 能在编译期计算的函数标记 constexpr
static constexpr size_t Size() noexcept { return dimension; }

// 不抛异常的函数标记 noexcept
T* Data() noexcept { return data.data(); }
const T* Data() const noexcept { return data.data(); }
```

### 4. 明确的接口
```cpp
// 使用 virtual 标记可重写函数
virtual T& operator[](size_t idx)
{
    assert(idx < dimension && "Index out of bounds");
    return data[idx];
}

// 使用 override 标记重写函数
virtual std::ostream& print(std::ostream& out) const override
{
    // 实现
}
```

### 5. 避免深层嵌套
```cpp
// ❌ 不推荐
if (condition1) {
    if (condition2) {
        if (condition3) {
            // 代码
        }
    }
}

// ✅ 推荐：早返回
if (!condition1) return;
if (!condition2) return;
if (!condition3) return;
// 代码
```

### 6. 清晰的错误处理
```cpp
// 使用 assert 处理编程错误
assert(idx < dimension && "Index out of bounds");

// 使用异常处理运行时错误
if (idx >= dimension) {
    throw std::out_of_range("Vector index out of range");
}

// 使用 static_assert 处理编译期错误
static_assert(rows > 0 && cols > 0, "Matrix dimensions must be positive");
```

---

## 文档规范

### 注释风格
```cpp
/**
 * 计算向量长度
 * @param dimensionalityReduction 是否降维计算（仅用于齐次坐标）
 * @return 向量的欧几里得长度
 */
virtual T length(bool dimensionalityReduction = false) const noexcept;

/*
* 数据访问
* 提供多种访问方式以适应不同场景
*/
virtual T& operator[](size_t idx);
```

### 代码分组
```cpp
class Vector
{
public:
    // 构造与析构
    Vector();
    ~Vector();

/*
* 数据访问
*/
    T& operator[](size_t idx);
    const T& operator[](size_t idx) const;

/*
* 向量运算
*/
    Vector operator+(const Vector& other) const;
    Vector operator-(const Vector& other) const;

/*
* 几何运算
*/
    T length() const;
    void normalize();
};
```

---

## 模板最佳实践

### SFINAE 使用
```cpp
// 方法1: 模板参数默认值
template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
class Vector { };

// 方法2: 返回类型 SFINAE
template<size_t R = rows, size_t C = cols>
typename std::enable_if_t<R == C, T> determinant() const;

// 方法3: 函数参数 SFINAE
template<bool rgbmode = isRGBMode>
typename std::enable_if_t<rgbmode, T> luminance();
```

### 模板特化
```cpp
// 偏特化用于优化特定情况
template<typename T>
T determinant() const
{
    if constexpr (rows == 1) {
        return data[0];
    } else if constexpr (rows == 2) {
        return data[0] * data[3] - data[1] * data[2];
    } else {
        // 通用实现
    }
}
```

---

## Unreal Engine 兼容性

### 条件编译
```cpp
#ifdef ENGINE_API
    // Unreal Engine 插件模式
    // 跳过 sandbox，使用 UE 导出宏
#else
    // 独立库模式
    #define EASYMATH_API  // 导出宏
#endif
```

### 宏命名
```cpp
// 避免与 UE 宏冲突
#define EASYMATH_API  // 而不是简单的 API
#define EM_INLINE inline  // 明确的命名空间前缀
```

---

## 版本控制

### Commit Message 格式
```
#Update <简短描述>
#Add <新增功能>
#Fix <修复问题>
#Refactor <重构>
#Docs <文档更新>
```

### 分支策略
- `master`: 主分支（稳定版本）
- `develop`: 开发分支
- `feature/*`: 功能分支
- `fix/*`: 修复分支
