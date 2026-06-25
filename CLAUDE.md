# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## 项目概述

EasyMath 是一个为实时图形与渲染场景打造的现代 C++ 数学库（命名空间 `EasyMath`），提供维度泛化的向量、矩阵、四元数、颜色类型，以及字节大小转换工具。强调编译期类型安全、零运行时开销（模板元编程 + `constexpr`）和直观的 API。

- 当前版本：1.0.11，正在开发 1.1.0
- 编译要求：C++17 或更高（Premake 脚本实际使用 C++20），MSVC v143+ / Clang 12+ / GCC 10+
- 库以 Header-Only 为主；仅 `src/ByteSize.cpp` 包含非模板实现

## 构建与运行

### 生成 Visual Studio 工程（Windows）
```bash
Scripts/GenerateWIN.bat
```
该脚本调用 `Scripts/premake/premake5.exe --file=Build.lua vs2022` 在仓库根目录生成 `EasyMath.sln`。

### 运行 Sandbox（运行时验证）
唯一运行时验证入口是 `Sandbox` 控制台程序：
```bash
# 1. 用 Visual Studio 打开 EasyMath.sln，选择 "Sandbox" 项目 → 生成 → 运行
#    或在仓库根目录：
msbuild EasyMath.sln /p:Configuration=Debug /p:Platform=x64

# 2. 直接执行
./bin/windows-x86_64/Debug/Sandbox.exe
```
输出 MVP（Model-View-Projection）管线结果到控制台，可作为 API 行为的活文档。

### 生成物位置
- 输出根目录：`bin/{system}-{arch}/{Debug|Release|Dist}/`
  - 库 DLL/LIB：`bin/windows-x86_64/Debug/EasyMath/`
  - Sandbox 可执行文件：`bin/windows-x86_64/Debug/Sandbox.exe`
- 中间对象文件：`bin/obj/{system}-{arch}/{Config}/{ProjectName}/`
- Visual Studio 配置：`settings/`

### Premake 配置结构
- `Build.lua` - 顶层 workspace 配置（x64，三种配置：Debug / Release / Dist）
- `EasyMath.lua` - 定义 `EasyMath` 项目（SharedLib，定义 `EASYMATH_DLL`）
- `Sandbox.lua` - 定义 `Sandbox` 控制台测试程序（C++20，链接 `EasyMath`）
- `Directory.lua` - 输出路径变量（`LocationDir` / `EXEDir` / `TargetDir` / `ObjectDir`）
- `Dependencies.lua` - 包含目录与库目录（`ProjIncludeDir["EasyMath"]` 指向 `include/`）
- `Config.lua` - Rider `.run` 与 Fork 自定义命令生成

### 重新生成项目文件
修改任何 `*.lua` 后必须重新运行 `GenerateWIN.bat`，否则 Visual Studio 工程不会更新。

### 开发文档站（demo/）
`demo/` 是一个独立的 Astro（Node.js）站点，用于文档展示与可视化：
```bash
cd demo
npm install        # 一次性
npm run dev        # 本地开发服务器
npm run build      # 静态站点构建到 demo/dist/
```

## 目录结构

```
E:\EasyMath
├─include/         公开头文件（核心库）
│  ├─EasyMath.h        头文件聚合入口，末尾 `using namespace EasyMath;`
│  ├─EasyMathAPI.h     DLL 导入/导出宏 (EASYMATH_API)
│  ├─EMConst.h         数学常量（PI）
│  ├─Common.h          Max/Min 变长模板 + NEARZERO_THRESHOLD
│  ├─Version.h         C++ 标准版本宏（CPP11/14/17/20_OR_LATER）
│  ├─Vector.h          Vector<T,N> + Swizzle<T,N,Indices...>
│  ├─Matrix.h          Matrix<T,R,C> + 投影/视图矩阵工厂
│  ├─Quaternion.h      Quaternion<T> : Vector<T,4>
│  ├─Color.h           Color<T,bIsRGB,modeIndex> + sRGB/Linear/HSV/HSL/HSI
│  ├─ByteSize.h        字节单位转换函数声明
│  └─EasyConversion.h  Vector ↔ Matrix 互转 + 提取/外积/Rodrigues
├─src/              非模板实现
│  └─ByteSize.cpp
├─sandbox/          独立测试程序源码
│  └─sandbox.cpp       MVP 管线测试（MTXLookAt/Perspective/Ortho）
├─Scripts/
│  ├─GenerateWIN.bat   Premake 一键生成 VS 工程
│  └─premake/premake5.exe
├─*.lua             Premake 构建脚本
├─EasyMath.sln      生成的 VS 解决方案
├─demo/             Astro (Node.js) 文档/演示站点
└─docs/             API 文档 Markdown（中英文）
   ├─Vector.md / Vector.en.md
   ├─Matrix.md / Matrix.en.md
   └─Swizzle.md / Swizzle.en.md
```

## 核心架构

### 1. 维度泛型 `Vector<T, N>`
- 存储：`std::array<T, N> data;`，栈分配、缓存友好
- 通过 `union` 将 `data` 与多个 `Swizzle` 成员共享内存，实现 **零开销** swizzle：
  ```cpp
  union {
      std::array<T, dimension> data;
      Swizzle<T, dimension, 0> x;
      Swizzle<T, dimension, 0, 1> xy;
      Swizzle<T, dimension, 2, 1, 0> zyx;  // 反向
  };
  ```
- 别名：`Vector2<T>` / `Vector3<T>` / `Vector4<T>`（在 `Vector.h` 末尾）
- 元素类型约束：`std::enable_if_t<std::is_arithmetic_v<T>>`
- 几何命名（位置 x/y/z/w 或颜色 r/g/b/a）通过 Swizzle 暴露

### 2. 行主序 `Matrix<T, R, C>`
- 存储布局：`data[row * cols + col]`（与 DirectX 一致，**OpenGL 是列主序**——交互时需注意）
- 平方操作支持 SFINAE：行列式 / 余子式 / 伴随 / 逆矩阵仅在 `R == C` 时可用
- 反演算法：伴随矩阵法（`A⁻¹ = adj(A) / det(A)`），递归求行列式
- 自由函数工厂：`MTXIdentity` / `MTXRotationX/Y/Z` / `MTXTranslation` / `MTXScale` / `MTXLookAt` / `MTXPerspective` / `MTXOrtho`
- 矩阵 × 向量返回 `Matrix<T, rows, 1>`（列向量）—— **`Vector3 × Matrix4x4` 等不直接支持**，需显式 `toHomogeneous()`

### 3. `Quaternion<T> : Vector<T, 4>`
- 当前完成度约 30%，是 v1.1.0 重点完善对象
- 暂只实现 length/normalize/print；缺失乘法、共轭、逆、轴角/欧拉角/矩阵互转、slerp、nlerp
- 输出顺序：`[w, x, y, z]`（注意！非 xyzw）

### 4. `Color<T, bIsRGB, modeIndex> : Vector<T, 4>`
- 用模板参数编码颜色空间状态，配合 SFINAE 实现**类型安全**的转换：
  | 别名 | bIsRGB | modeIndex | 可用方法 |
  |------|--------|-----------|----------|
  | `sRGBColor<T>` | true | 1 | `toLinear()`, `toHSV/HSL/HSI()` |
  | `LinearColor<T>` | true | 0 | `toSRGB()`, `toHSV/HSL/HSI()` |
  | `HSV<T>` / `HSL<T>` / `HSI<T>` | false | 0/1/2 | `fromHSV/H/H()` 静态构造 |
- 不允许的转换在编译期被阻止（如 HSV 调用 `toHSV()` 编译错误）
- RGB → HSx 是实例方法，HSx → RGB 是静态工厂 `Color::fromHSV(h,s,v,a)`

### 5. Swizzle 机制
- 模板签名 `Swizzle<T, dimension, int... Indices>`，参数是**索引序列**
- 通过 `reinterpret_cast` 直接访问 `union` 中的 `data` 数组（`elem(i)` 辅助函数）
- 在 `Vector` 和 `Swizzle` 之间的所有算术运算都已同步（参见最近的 commit "Sync algorithms of vector to swizzle"）

### 6. 投影管线（`sandbox/sandbox.cpp` 演示）
完整 MVP 管线：
```cpp
auto model = MTXTranslation(...) * MTXRotationY(...);
auto view  = MTXLookAt(eye, target, up);
auto proj  = MTXPerspective<float>(fov, aspect, near, far);
auto MVP   = proj * view * model;
```

## 关键约定

- **命名空间**：所有公开符号在 `EasyMath`（无 `EM` 简写别名；头文件末尾 `using namespace EasyMath;`）
- **命名**：
  - 类/类型：`PascalCase`（`Vector`, `Matrix`, `Quaternion`）
  - 成员函数 / 全局函数：`camelCase`（`normalized`, `dot`, `cross`）
  - 矩阵工厂：`MTX` 前缀 + `PascalCase`（`MTXIdentity`, `MTXRotationX`）
  - 颜色空间别名：首字母小写（`sRGBColor`, `linearColor`）
  - 常量：`UPPER_SNAKE_CASE`（`PI`, `NEARZERO_THRESHOLD`）
  - 枚举值：`PascalCase`（`LuminanceStandard::Rec709`）
- **缩进**：Tab（与仓库现有风格保持一致）
- **大括号**：K&R 风格（类/函数左括号换行，控制流同行）
- **错误处理**：
  - 编程错误用 `assert`
  - 运行时错误用 `throw std::out_of_range`（`at()`）
  - 编译期错误用 `static_assert` 给出友好消息
- **宏约定**：避免与 Unreal Engine 冲突——使用 `EASYMATH_API` 而非 `API`
- **Unreal 集成**：通过 `#ifdef ENGINE_API` 条件编译切换为 UE 插件模式
- **Commit Message 前缀**：`#Add` / `#Update` / `#Fix` / `#Refactor` / `#Docs` / `#Delete`

## 重要参考资料

- `docs/Vector.md` / `docs/Vector.en.md` - Vector 与 Swizzle 详细 API
- `docs/Matrix.md` / `docs/Matrix.en.md` - Matrix API 与投影矩阵推导
- `docs/Swizzle.md` / `docs/Swizzle.en.md` - Swizzle 运算符与边界
- `.claude/memory/design-principles.md` - Swizzle/Color 类型安全机制详解
- `.claude/memory/code-conventions.md` - 完整代码规范
- `.claude/memory/development-status.md` - 各模块完成度（Vector/Matrix/Color 95%、Quaternion 30%）
- `.claude/memory/roadmap.md` - v1.1.0 待办（Quaternion 完善 / 欧拉角 / 旋转互转 / Color 运算）
- `.claude/memory/optimization-summary.md` - 投影矩阵性能优化细节
- `.claude/memory/matrix-optimization-analysis.md` - 矩阵优化分析
- `.claude/memory/project-context.md` - 项目背景与定位
- `.claude/DEMO_INTEGRATION.md` - 文档站（Astro）与库源码的协作约定
- `.claude/easymath-design-system.md` - 文档站视觉/设计规范
- `.claude/agents/` `.claude/skills/` `.claude/mcp/` - Claude Code 配套资源（agents、slash 命令、MCP 配置）
- `README.md` - 项目说明、构建步骤、特性列表（中英文）

## 验证

- **沙盒测试**：`sandbox/sandbox.cpp` 是当前唯一的运行时测试入口，编译运行后输出 MVP 管线结果到控制台
- **类型测试**：通过 SFINAE 在编译期捕获错误（如错误颜色空间调用、错误矩阵维度）
- **没有自动化单元测试框架**（v1.1.0 路线图中计划添加）
- **没有 Lint 配置**（`warnings "off"` 已在 `EasyMath.lua` 中）

## 已知陷阱

- **Matrix × Vector**：矩阵乘向量返回 `Matrix<T, rows, 1>`（列向量），不是 `Vector`；要得到 Vector3 结果需显式 `.template get<>()` 或提取分量
- **`Vector3 × Matrix4x4`** 不直接支持：先把 Vector3 通过 `toHomogeneous()` 提升到 Vector4
- **Matrix 存储顺序**：`data[row * cols + col]` 行主序（DirectX 风格），与 GLSL 列主序相反
- **Quaternion 输出顺序**：`[w, x, y, z]`，不是 xyzw
