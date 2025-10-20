# EasyMath

[![Static Badge](https://img.shields.io/badge/Language-English-purple)](#english) | [![Static Badge](https://img.shields.io/badge/Language-中文-red)](#中文)

---

![Static Badge](https://img.shields.io/badge/Version-1.0.10-orange) ![Static Badge](https://img.shields.io/badge/Code->C%2B%2B17-red) ![Static Badge](https://img.shields.io/badge/Export-dll-brightgreen) ![Static Badge](https://img.shields.io/badge/Sandbox-console-brown) ![Static Badge](https://img.shields.io/badge/State-developing-blue)

--- 

## English

EasyMath is a modern C++ math library for real-time graphics and rendering. It provides dimension-generic vectors and matrices, quaternions, colors, and practical utilities like byte-size conversion. The API emphasizes type safety, readability, and ease of integration.

- Namespace: `EM`
- Core types are header-first; a small amount of implementation lives in `src/`
- Recommended standard: C++17 or newer

### Directory layout
```text
E:.\EasyMath
├─include
│  ├─ByteSize.h
│  ├─Color.h
│  ├─Common.h
│  ├─EasyConversion.h
│  ├─EasyMath.h
│  ├─EasyMathAPI.h
│  ├─EMConst.h
│  ├─Matrix.h
│  ├─Quaternion.h
│  ├─Vector.h
│  └─Version.h
├─src
│  └─ByteSize.cpp
├─sandbox
│  └─sandbox.cpp
├─Scripts
│  ├─GenerateWIN.bat
│  └─premake
│      └─premake5.exe
├─Build.lua
├─Config.lua
├─Dependencies.lua
├─Directory.lua
├─EasyMath.lua
├─Sandbox.lua
├─LICENSE
└─README.md
```

### Build & Integrate (WIN)
- Run `Scripts/GenerateWIN.bat` to generate/configure the Visual Studio solution on Windows
- Open `EasyMath.sln`, build and run `Sandbox` to see the demo
- Generated DLL and import library are under `bin/windows-x86_64/Debug/EasyMath/` (e.g., `EasyMath.dll`, `.lib`) for linking from other projects

### Features
- Namespace: `EM`
- Vectors: `Vector<T, N>` with arithmetic, swizzle, length/normalize, dot/cross, distance/angle, lerp/slerp, homogeneous conversions, projection/reflection, skew-symmetric matrix, and common aliases (`Vector2/3/4`, `Vector2f/d/i`, ...)
- Matrices: `Matrix<T, R, C>` with element access, transpose, determinant/cofactor/adjugate/inverse, arithmetic (+, unary -, scalar ×), matrix×matrix, matrix×vector, and helpers (identity, rotations X/Y/Z, translation, scale)
- Quaternions: `Quaternion<T>` based on `Vector<T,4>`, with length, normalization and output
- Color: `Color<T>` based on `Vector<T,4>`, with luminance helpers
- Byte size utilities: `ByteSizeTo`, `ByteSizeConvert`

### Quick Start
```cpp
#include <iostream>
#include "EasyMath.h"
#include "Vector.h"
#include "Matrix.h"
#include "Quaternion.h"
#include "Color.h"

int main() {
    using namespace EM;

    Vector3 a{1.0f, 2.0f, 3.0f};
    Vector3 b{2.0f, 0.0f, 1.0f};

    float d = dot(a, b);
    auto c = cross(a, b);
    auto an = a.normalized();

    Matrix<float, 4, 4> rx = MTXRotationX<float, 4>(3.14159f * 0.5f);
    auto v4 = a.toHomogeneous();
    auto rotated = rx * v4;

    QuaternionF q{0, 0, 0, 1};
    auto qn = q.normalized();

    Color<float> col{0.2f, 0.7f, 0.1f, 1.0f};
    float L = col.luminance();

    std::cout << "dot=" << d << " cross=" << c << " L=" << L << "\n";
}
```

### Documentaion
- [Vector](/docs/Vector.en.md)
- [Swizzle](/docs/Swizzle.en.md)
- [Matrix](/docs/Matrix.en.md)

### Roadmap
![Static Badge](https://img.shields.io/badge/next_version-1.1-yellow)


- [ ] Optimize matrix and vector algorithms
- [ ] Enrich quaternion APIs
- [ ] Add Euler angles APIs
- [ ] Add conversions among Euler angles, quaternions and axis-angle
- [x] Optimize algorithm of vector swizzle operation
- [x] Sync algorithms of vector to swizzle

### Version & Requirements
- Compilers: MSVC v143+, Clang 12+, GCC 10+ (recommended)
- C++17 or newer

### Contributing
- Issues and PRs are welcome. Keep naming and semantics consistent. Prefer clear naming and early returns; avoid deep nesting.

### License
Please see the English license `LICENSE`. A Chinese reference translation is available in `LICENSE.zh-CN`.

## 中文

EasyMath 是一个为实时图形与渲染场景打造的现代 C++ 数学库。它提供维度泛化的向量与矩阵、四元数、颜色类型，以及实用的字节大小转换工具。API 强调类型安全、可读性与高可用性。

 - 命名空间：`EM`
 - 核心类型以头文件为主，少量实现位于 `src/`
 - 建议使用 C++17 及以上

### 目录结构
```text
E:.\EasyMath
├─include
│  ├─ByteSize.h
│  ├─Color.h
│  ├─Common.h
│  ├─EasyConversion.h
│  ├─EasyMath.h
│  ├─EasyMathAPI.h
│  ├─EMConst.h
│  ├─Matrix.h
│  ├─Quaternion.h
│  ├─Vector.h
│  └─Version.h
├─src
│  └─ByteSize.cpp
├─sandbox
│  └─sandbox.cpp
├─Scripts
│  ├─GenerateWIN.bat
│  └─premake
│      └─premake5.exe
├─Build.lua
├─Config.lua
├─Dependencies.lua
├─Directory.lua
├─EasyMath.lua
├─Sandbox.lua
├─LICENSE
└─README.md
```

### 构建与集成（WIN）
- 运行 `Scripts/GenerateWIN.bat` 即可在 Windows 上生成/配置 Visual Studio 工程
- 打开 `EasyMath.sln`，构建并运行 `Sandbox` 查看演示
- 生成的 DLL 与导入库位于 `bin/windows-x86_64/Debug/EasyMath/`（如 `EasyMath.dll`、`.lib`），可直接供其他项目链接使用

### 特性
- 命名空间：`EM`
- 向量：`Vector<T, N>`，支持算术运算、swizzle、长度与归一化、点积/叉积、距离/角度、线性/球面插值、齐次坐标转换、投影/反射、反对称矩阵、常用别名（`Vector2/3/4`、`Vector2f/d/i` 等）
- 矩阵：`Matrix<T, R, C>`，支持元素访问、转置、行列式/代数余子式/伴随/逆，算术运算（加减、取负、标量乘），矩阵×矩阵、矩阵×向量，常用构造（单位矩阵、绕 X/Y/Z 旋转、平移、缩放）
- 四元数：`Quaternion<T>` 基于 `Vector<T,4>`，提供长度、归一化与输出等功能
- 颜色：`Color<T>` 基于 `Vector<T,4>`，提供亮度计算等功能
- 字节大小：`ByteSizeTo`、`ByteSizeConvert`



### 快速开始
```cpp
#include <iostream>
#include "EasyMath.h"
#include "Vector.h"
#include "Matrix.h"
#include "Quaternion.h"
#include "Color.h"

int main() {
    using namespace EM;

    Vector3 a{1.0f, 2.0f, 3.0f};
    Vector3 b{2.0f, 0.0f, 1.0f};

    float d = dot(a, b);
    auto c = cross(a, b);
    auto an = a.normalized();

    Matrix<float, 4, 4> rx = MTXRotationX<float, 4>(3.14159f * 0.5f);
    auto v4 = a.toHomogeneous();
    auto rotated = rx * v4;

    QuaternionF q{0, 0, 0, 1};
    auto qn = q.normalized();

    Color<float> col{0.2f, 0.7f, 0.1f, 1.0f};
    float L = col.luminance();

    std::cout << "dot=" << d << " cross=" << c << " L=" << L << "\n";
}
```
### 文档
- [Vector](/docs/Vector.md)
- [Swizzle](/docs/Swizzle.md)
- [Matrix](/docs/Matrix.md)

### 版本开发计划 
![Static Badge](https://img.shields.io/badge/Next_Version-1__1__0-yellow)
- [ ] 优化矩阵和向量算法
- [ ] 完善四元数计算API
- [ ] 添加欧拉角API
- [ ] 添加欧拉角、四元数、轴角转换方法
- [x] 优化vector swizzle算法
- [x] 为swizzle同步vector算法


### 版本与要求
- 编译器：MSVC v143+、Clang 12+、GCC 10+（推荐）
- C++17 或更高

### 参与贡献
- 欢迎 Issue 与 PR。请保持 API 命名与语义一致，优先清晰命名与早返回，避免深层嵌套。

### 许可协议
请参阅英文许可 `LICENSE`
