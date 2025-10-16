# `Vector`类

[![Static Badge](https://img.shields.io/badge/Back_to_README-gray)](../README.md "点击返回README") | [![Static Badge](https://img.shields.io/badge/Language-zh-red)](Vector.en.md "click to english version")

---

![Static Badge](https://img.shields.io/badge/Type-class-green)
![Static Badge](https://img.shields.io/badge/State-completed-blue)



Vector是一个长度可变、类型可变的泛型向量

| Attribute     | Value                                                              |
| ------------- | ------------------------------------------------------------------ |
| **Namespace** | **`EM`**                                                           |
| **File**      | [**`include\Vector.h`**](../include/Vector.h)                      |
| **Parent**    | -                                                                  |
| **Feature**   | ![Static Badge](https://img.shields.io/badge/Feature-template-red) |

## 声明

```C++
template<typename T, size_t dimension, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
class Vector
```

### 模板

| Type     | Name        | Description                                            |
| -------- | ----------- | ------------------------------------------------------ |
| typename | `T`         | 元素类型，受到 `std::is_arithmetic_v` （算数类型）限制 |
| size_t   | `dimension` | 向量维度/向量长度                                      |

### 构建函数

| Access   | Modifier                                                                       | Return | Name                                   | Parameter                                                       | Description                                                      |
| -------- | ------------------------------------------------------------------------------ | ------ | -------------------------------------- | --------------------------------------------------------------- | ---------------------------------------------------------------- |
| `public` | -                                                                              |        | `Vector()`                             | -                                                               | 默认构建，全部元素初始化为0                                      |
| `public` | `explicit`                                                                     |        | `Vector(const T& value)`               | `const T& value`-类型为 `T`的单一值                             | 使用单一输入的值初始化向量的所有元素                             |
| `public` | -                                                                              |        | `Vector(std::initializer_list<T>)`     | `std::initializer_list<T> InitializeList`-类型为T的初始化数列   | 使用数列初始化向量元素如果数列元素数量不等于 `dimension`则会报错 |
| `public` | `template<typename T2, typename = std::enable_if_t<std::is_arithmetic_v<T2>>>` |        | `Vector(const Vector<T2, dimension>&)` | `const Vector<T2, dimension>& other`-另一不同类型但同维度的向量 | -                                                                |

```C++
EM::Vector<float,2> vectorA;
EM::Vector<float,3> vectorB(1);
EM::Vector<float,4> vectorC{ 1, 2, 3, 4};
```

## 数据访问

### 一般数据访问

| Access   | Modifier                                 | Return     | Name                 | Parameter                 | Description                                  |
| -------- | ---------------------------------------- | ---------- | -------------------- | ------------------------- | -------------------------------------------- |
| `public` | `virtual`                                | `T&`       | `Operator[](size_t)` | `size_t idx`-元素索引位置 | -                                            |
| `public` | `virtual` `const`                        | `const T&` | `Operator[](size_t)` | `size_t idx`-元素索引位置 | -                                            |
| `public` | `virtual`                                | `T&`       | `at(size_t)`         | `size_t idx`-元素索引位置 | -                                            |
| `public` | `virtual` `const`                        | `const T&` | `at(size_t)`         | `size_t idx`-元素索引位置 | -                                            |
| `public` | `virtual` `noexcept`                     | `T*`       | `Data()`             | -                         | 直接返回data数列的指针                       |
| `public` | `virtual` `const` `noexcept`             | `const T*` | `Data()`             | -                         | 直接返回data数列的指针                       |
| `public` | `virtual` `noexcept`                     | `T*`       | `begin()`            | -                         | 返回data数列第一个元素的指针                 |
| `public` | `virtual` `const` `noexcept`             | `const T*` | `begin()`            | -                         | 返回data数列第一个元素的指针                 |
| `public` | `virtual` `noexcept`                     | `T*`       | `end()`              | -                         | 返回data数列的末尾指针（最后一个元素索引+1） |
| `public` | `virtual` `const` `noexcept`             | `const T*` | `end()`              | -                         | 返回data数列的末尾指针（最后一个元素索引+1） |
| `public` | `virtual` `constexpr` `const` `noexcept` | `size_t`   | `getDimension()`     | -                         | 返回向量的维度（数组的维度）                 |
| `public` | `static` `constexpr` `noexcept`          | `size_t`   | `Size()`             | -                         | 返回向量的维度（数组的维度）                 |
| `public` | `static` `constexpr` `noexcept`          | `size_t`   | `Dimension()`        | -                         | 返回向量的维度（数组的维度）                 |


### Swizzle
向量元素的访问支持像glsl、hlsl等语言中使用xyzw或者rgba访问和修改元素
```C++
EM::Vector<float,4> vectorA{1,2,3,4};
std::cout << vectorA << std::endl;

vectorA.x = 5;
vectorA.zy = EM::Vector<float,2> {6,7};
std::cout << vectorA << std::endl;

std::cout << vectorA.xyz << std::endl;
```
```
// 输出结果
(1, 2, 3, 4) (V mode: xyzw)
(5, 7, 6, 4) (V mode: xyzw)
(5, 7, 6) (V mode: xyz)
```

## 运算支持

### 复合赋值运算
| Access   | Modifier  | Return                  | Name                                      | Parameter                                                    | Description |
| -------- | --------- | ----------------------- | ----------------------------------------- | ------------------------------------------------------------ | ----------- |
| `public` | `virtual` | `Vector<T, dimension>&` | `operator=(const Vector<T, dimension>&)`  | `const Vector<T, dimension>& other`-同类型和同维度的另一向量 | -           |
| `public` | `virtual` | `Vector<T, dimension>&` | `operator=(Vector<T, dimension>&)`        | `Vector<T, dimension>& other`-同类型和同维度的另一向量       | -           |
| `public` | `virtual` | `Vector<T, dimension>&` | `operator+=(const Vector<T, dimension>&)` | `const Vector<T, dimension>& other`-同类型和同维度的另一向量 | -           |
| `public` | `virtual` | `Vector<T, dimension>&` | `operator+=(const T&)`                    | `const T& scalar`-同类型的常量                               | -           |
| `public` | `virtual` | `Vector<T, dimension>&` | `operator-=(const Vector<T, dimension>&)` | `const Vector<T, dimension>& other`-同类型和同维度的另一向量 | -           |
| `public` | `virtual` | `Vector<T, dimension>&` | `operator-=(const T&)`                    | `const T& scalar`-同类型的常量                               | -           |
| `public` | `virtual` | `Vector<T, dimension>&` | `operator*=(const Vector<T, dimension>&)` | `const Vector<T, dimension>& other`-同类型和同维度的另一向量 | -           |
| `public` | `virtual` | `Vector<T, dimension>&` | `operator*=(const T&)`                    | `const T& scalar`-同类型的常量                               | -           |
| `public` | `virtual` | `Vector<T, dimension>&` | `operator/=(const Vector<T, dimension>&)` | `const Vector<T, dimension>& other`-同类型和同维度的另一向量 | -           |
| `public` | `virtual` | `Vector<T, dimension>&` | `operator/=(const T&)`                    | `const T& scalar`-同类型的常量                               | -           |

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

### 二元运算
| Access             | Modifier | Return                 | Name                                                                  | Parameter                                                                         | Description |
| ------------------ | -------- | ---------------------- | --------------------------------------------------------------------- | --------------------------------------------------------------------------------- | ----------- |
| `private` `friend` | -        | `Vector<T, dimension>` | `operator+(const Vector<T, dimension>&, const Vector<T, dimension>&)` | `const Vector<T, dimension>& lhs`-向量,<br>`const Vector<T, dimension>& rhs`-向量 | -           |
| `private` `friend` | -        | `Vector<T, dimension>` | `operator+(const Vector<T, dimension>&, const T&)`                    | `const Vector<T, dimension>& vec`-向量,<br>`const T& scalar`-标量                 | -           |
| `private` `friend` | -        | `Vector<T, dimension>` | `operator+(const T&，const Vector<T, dimension>&)`                    | `const T& scalar`-标量，<br>`const Vector<T, dimension>& vec`-向量                | -           |
| `private` `friend` | -        | `Vector<T, dimension>` | `operator-(const Vector<T, dimension>&, const Vector<T, dimension>&)` | `const Vector<T, dimension>& lhs`-向量,<br>`const Vector<T, dimension>& rhs`-向量 | -           |
| `private` `friend` | -        | `Vector<T, dimension>` | `operator-(const Vector<T, dimension>&, const T&)`                    | `const Vector<T, dimension>& vec`-向量,<br>`const T& scalar`-标量                 | -           |
| `private` `friend` | -        | `Vector<T, dimension>` | `operator*(const Vector<T, dimension>&, const Vector<T, dimension>&)` | `const Vector<T, dimension>& lhs`-向量,<br>`const Vector<T, dimension>& rhs`-向量 | -           |
| `private` `friend` | -        | `Vector<T, dimension>` | `operator*(const Vector<T, dimension>&, const T&)`                    | `const Vector<T, dimension>& vec`-向量,<br>`const T& scalar`-标量                 | -           |
| `private` `friend` | -        | `Vector<T, dimension>` | `operator*(const T&, const Vector<T, dimension>&)`                    | `const T& scalar`-标量,<br>`const Vector<T, dimension>& vec`-向量                 | -           |
| `private` `friend` | -        | `Vector<T, dimension>` | `operator/(const Vector<T, dimension>&, const Vector<T, dimension>&)` | `const Vector<T, dimension>& lhs`-向量,<br>`const Vector<T, dimension>& rhs`-向量 | -           |
| `private` `friend` | -        | `Vector<T, dimension>` | `operator/(const Vector<T, dimension>&, const T&)`                    | `const Vector<T, dimension>& vec`-向量,<br>`const T& scalar`-标量                 | -           |

```C++
EM::Vector<float,2> vectorA{ 1, 2};
EM::Vector<float,2> vectorB{ 3, 4};

std::cout << vectorA+vectorB << std::endl;
std::cout << vectorA-vectorB << std::endl;
std::cout << vectorA*vectorB << std::endl;
std::cout << vectorA/vectorB << std::endl;
```


### 取反
| Access             | Modifier | Return                 | Name                                     | Parameter                              | Description      |
| ------------------ | -------- | ---------------------- | ---------------------------------------- | -------------------------------------- | ---------------- |
| `private` `friend` | -        | `Vector<T, dimension>` | `operator-(const Vector<T, dimension>&)` | `const Vector<T, dimension>& vec`-向量 | 向量支持直接取反 |

```C++
EM::Vector<float,2> vectorA{ 1, 2};
std::cout << -vectorA << std::endl;
```

### 比较运算
| Access             | Modifier | Return | Name                                                                   | Parameter                                                                         | Description |
| ------------------ | -------- | ------ | ---------------------------------------------------------------------- | --------------------------------------------------------------------------------- | ----------- |
| `private` `friend` | -        | `bool` | `operator==(const Vector<T, dimension>&, const Vector<T, dimension>&)` | `const Vector<T, dimension>& lhs`-向量,<br>`const Vector<T, dimension>& rhs`-向量 | -           |
| `private` `friend` | -        | `bool` | `operator!=(const Vector<T, dimension>&, const Vector<T, dimension>&)` | `const Vector<T, dimension>& lhs`-向量,<br>`const Vector<T, dimension>& rhs`-向量 | -           |


## 输出
| Access   | Modifier                                 | Return          | Name                                                    | Parameter                                                                | Description                        |
| -------- | ---------------------------------------- | --------------- | ------------------------------------------------------- | ------------------------------------------------------------------------ | ---------------------------------- |
| `public` | `virtual` `const`                        | `std::ostream&` | `print(std::ostream&)`                                  | `std::ostream& out`-输出流                                               | 在控制台输出（这个函数为成员函数） |
| -        | `template<typename T, size_t dimension>` | `std::ostream&` | `operator<<(std::ostream&, const Vector<T,dimension>&)` | `std::ostream& out`-输出流,<br>`const Vector<T,dimension>& vec`-输出向量 | 在控制台输出（这个函数为全局函数） |

```C++
EM::Vector<float,2> vectorA{8,9};
std::cout << vectorA << std::endl;
```

## 成员方法
| Access   | Modifier                                                | Return                            | Name                                   | Parameter                                                                                                                                | Description                                          |
| -------- | ------------------------------------------------------- | --------------------------------- | -------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------- |
| `public` | `virtual`                                               | `Matrix<T, dimension, 1>`         | `toColMatrix()`                        | -                                                                                                                                        | 用列矩阵表示向量                                     |
| `public` | `virtual`                                               | `Matrix<T, 1, dimension>`         | `toRowMatrix()`                        | -                                                                                                                                        | 用行矩阵表示向量                                     |
| `public` | `virtual` `constexpr`                                   | `T`                               | `lengthSquared(bool)`                  | `bool dimensionalityReduction = false`-是否要进行降维运算<br>降维将只计算前3个维度                                                       | 计算向量长度的平方                                   |
| `public` | `virtual` `const` `noexcept`                            | `T`                               | `length(bool)`                         | `bool dimensionalityReduction = false`-是否要进行降维运算<br>降维将只计算前3个维度                                                       | 计算向量长度                                         |
| `public` | `virtual` `const`                                       | `Vector<T, dimension>`            | `normalized(bool)`                     | `bool dimensionalityReduction = false`-是否要进行降维运算<br>降维将只计算前3个维度                                                       | 计算向量归一化的结果                                 |
| `public` | `virtual`                                               | `Vector<T, dimension>`            | `normalize(bool)`                      | `bool dimensionalityReduction = false`-是否要进行降维运算<br>降维将只计算前3个维度                                                       | 向量归一化                                           |
| `public` | `virtual` `const` `noexcept`                            | `bool`                            | `isNormalized(T,bool)`                 | `bool dimensionalityReduction = false`-是否要进行降维运算<br>降维将只计算前3个维度<br>`T epsilon = T{NEARZERO_THRESHOLD}`-浮点数判断精度 | 是否为单位向量                                       |
| `public` | `virtual` `const` `noexcept`                            | `bool`                            | `isZero(T,bool)`                       | `T epsilon = T{NEARZERO_THRESHOLD}`-浮点数判断精度                                                                                       | 是否为0向量                                          |
| `public` | `virtual` `const` `noexcept`                            | `Vector<T, dimension>`            | `lerp(const Vector<T, dimension>&, T)` | `const Vector<T, dimension>& other`-插值向量,<br> `T t`-插值权重                                                                         | 插值函数                                             |
| `public` | `virtual` `const`                                       | `Vector<T, dimension>`            | `reflect(const Vector<T, dimension>&)` | `const Vector<T, dimension>& normal`-法向量（必须归一化）                                                                                | 基于法线的反射向量                                   |
| `public` | `virtual` `const`                                       | `Vector`                          | `project(const Vector&)`               | `const Vector& onto`-投影向量                                                                                                            | 在某一向量上的投影向量                               |
| `public` | `virtual` `const`                                       | `Vector`                          | `project(const Vector&, bool)`         | `const Vector& onto`-投影向量,<br>`bool dimensionalityReduction = false`-是否要进行降维运算<br>降维将只计算前3个维度<br>                 | 在某一向量上的投影向量                               |
| `public` | `template<size_t newDim = dimension + 1>` `const`       | `Vector<T, newDim>`               | `toHomogeneous(T)`                     | `T w = T{ 1 }`-齐次坐标的w分量（点向量 w=1, 方向向量w=0）                                                                                | 转换当前坐标为齐次坐标                               |
| `public` | `template<size_t newDim = dimension - 1>` `const`       | `Vector<T, newDim>`               | `fromHomogeneous()`                    | -                                                                                                                                        | 将当前向量作为齐次坐标转换为对应的向量               |
| `public` | `template<typename = std::enable_if_t<dimension == 3>>` | `Matrix<T, 4, 4>`                 | `toTranslationMatrix(bool)`            | `bool usedWithOrient = false`-是否用于移动方向向量（这个条件可以被优化掉）                                                               | 将当前向量转换为平移矩阵                             |
| `public` | `const`                                                 | `Matrix<T, dimension, dimension>` | `skewSymmetric()`                      | -                                                                                                                                        | 获得当前向量的反对称矩阵                             |
| `public` | `template<size_t D = dimension>` `const`                | `Matrix<T, 3, 3>`                 | `skewSymmetric_2D()`                   | -                                                                                                                                        | 获得当前2D向量的反对称矩阵(非2D向量不可访问这个函数) |

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


## 全局方法
| Access | Modifier                         | Return         | Name                                                              | Parameter                                                                                                                                             | Description  |
| ------ | -------------------------------- | -------------- | ----------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- | ------------ |
| -      | `template<typename T, size_t D>` | `T`            | `dot(const Vector<T, D>&, const Vector<T, D>&)`                   | `const Vector<T, D>& a`-向量A,<br> `const Vector<T, D>& b`-向量B                                                                                      | 点乘         |
| -      | `template<typename T>`           | `T`            | `cross(const Vector<T, 3>&, const Vector<T, 3>&)`                 | `const Vector<T, 3>& a`-3D向量A,<br> `const Vector<T, 3>& b`-3D向量B                                                                                  | 3D叉乘       |
| -      | `template<typename T>`           | `T`            | `cross(const Vector<T, 2>&, const Vector<T, 2>&)`                 | `const Vector<T, 2>& a`-2D向量A,<br> `const Vector<T, 2>& b`-2D向量B                                                                                  | 2D叉乘       |
| -      | `template<typename T, size_t D>` | `T`            | `distance(const Vector<T, D>&, const Vector<T, D>&, bool)`        | `const Vector<T, D>& a`-向量A,<br> `const Vector<T, D>& b`-向量B,<br> `bool dimensionalityReduction=true`-是否要进行降维运算<br>降维将只计算前3个维度 | 距离         |
| -      | `template<typename T, size_t D>` | `T`            | `distanceSquared(const Vector<T, D>&, const Vector<T, D>&, bool)` | `const Vector<T, D>& a`-向量A,<br> `const Vector<T, D>& b`-向量B,<br> `bool dimensionalityReduction=true`-是否要进行降维运算<br>降维将只计算前3个维度 | 距离的平方   |
| -      | `template<typename T, size_t D>` | `Vector<T, D>` | `lerp(const Vector<T, D>&, const Vector<T, D>&, T)`               | `const Vector<T, D>& a`-向量A,<br> `const Vector<T, D>& b`-向量B,<br> `T t`-插值权重                                                                  | 线性插值     |
| -      | `template<typename T, size_t D>` | `Vector<T, D>` | `slerp(const Vector<T, D>&, const Vector<T, D>&, T)`              | `const Vector<T, D>& a`-向量A,<br> `const Vector<T, D>& b`-向量B,<br> `T t`-插值权重                                                                  | 球面线性插值 |

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

[![Static Badge](https://img.shields.io/badge/Back_to_Top-gray)](#vector类 "点击返回顶部") | [![Static Badge](https://img.shields.io/badge/Back_to_README-gray)](../README.md "点击返回README")