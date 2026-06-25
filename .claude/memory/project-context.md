# EasyMath 项目上下文

## 项目身份
- **项目名称**: EasyMath
- **当前版本**: 1.0.11
- **目标版本**: 1.1.0
- **开发者**: CaptainRex
- **项目类型**: C++17 实时渲染数学库
- **命名空间**: `EasyMath`

## 核心定位
为实时图形与渲染场景打造的现代 C++ 数学库，支持：
1. 独立库使用
2. Unreal Engine 插件集成

## 项目结构
```
E:\Private\EasyMath
├─ include/           # 头文件（核心实现，header-only）
│  ├─ Vector.h        # 泛型向量 Vector<T,N>
│  ├─ Matrix.h        # 泛型矩阵 Matrix<T,R,C>
│  ├─ Quaternion.h    # 四元数 Quaternion<T>
│  ├─ Color.h         # 颜色类 Color<T>
│  ├─ EasyMath.h      # 主入口头文件
│  ├─ EasyMathAPI.h   # 导出宏与 UE 兼容
│  ├─ Common.h        # 通用定义
│  ├─ EMConst.h       # 数学常量
│  ├─ ByteSize.h      # 字节工具
│  ├─ EasyConversion.h # 转换工具
│  └─ Version.h       # 版本信息
├─ src/               # 源文件（少量实现）
│  └─ ByteSize.cpp
├─ sandbox/           # 测试沙盒
│  └─ sandbox.cpp
├─ docs/              # API 文档（中英文）
├─ Scripts/           # 构建脚本
└─ .claude/           # AI 辅助开发配置
```

## Unreal Engine 集成机制
通过 `ENGINE_API` 宏实现双模式：
```cpp
#ifdef ENGINE_API
    // UE 插件模式：跳过 sandbox，使用 UE 模块系统
#else
    // 独立模式：编译 sandbox 测试
#endif
```

## CPU / GPU 共用算法（2026-06-24 讨论结论）

**核心结论**：算法语义可以共用，实现不能。

| 层面 | 是否共用 |
|--|--|
| 数学语义（公式、定义） | ✅ 共用 |
| API 命名风格（`MTXLookAt(eye, target, up)`） | ✅ 共用 |
| 算法文档 / 教学概念 | ✅ 共用 |
| CPU 字节码 / GPU 指令 | ❌ 不共用 |
| C++ 模板实例 vs GLSL shader 代码 | ❌ 不共用 |

**为什么不能"一份代码两边跑"**：
- C++ 模板在编译期实例化，类型 `Matrix<float,4,4>` 是 C++ 对象
- GLSL 的 `mat4` 是另一套类型系统，编译进 GPU 指令
- CPU 计算结果通过 `glUniformMatrix4fv` / UBO 上传给 GPU，**中间是数据交换，不是代码共享**

**实际工程中的三种模式**：
1. **手写两份，靠纪律同步** — C++ 一份，GLSL 一份，公式必须一致。Unreal / Unity 大部分代码这么做。
2. **代码生成** — 写一份 DSL，工具链吐 C++ + GLSL + HLSL + MSL。工业级方案。
3. **CPU SIMD stub + GPU 精确实现** — CPU 端用 intrinsics，GPU 端用 GLSL 精确版。

**EasyMath 现状**：模式 1 的起点。已命名 `EasyMath::MTXLookAt` / `EasyMath::MTXPerspective` 等 CPU 实现。

**未来若加 GLSL 端**（最低成本路径）：
- 新建 `shaders/em_math.glsl`
- 从 `include/Matrix.h` 复制 `MTXLookAt / MTXPerspective / MTXRotationX/Y/Z / MTXTranslation / MTXScale` 的 GLSL 实现
- 命名一致（`emLookAt` / `emPerspective`）
- 文件头注释："This file mirrors include/Matrix.h, keep formulas in sync"
- 不进 Premake 编译（着色器是字符串资源），但可加 Python 脚本做两边公式对照检查

**反面清单（明确不共用的）**：
- GPU 着色器里的特殊优化（分支避免、向量合并）
- CPU 端的 SIMD 加速
- GLSL 精度限定（`mediump` / `highp`）
- 模板元编程（GLSL 不支持）

## 技术栈
- **语言标准**: C++17+
- **构建系统**: Premake5 (Lua)
- **编译器**: MSVC v143+ / Clang 12+ / GCC 10+
- **依赖**: 无外部依赖（仅 C++ 标准库）

## 设计原则
1. 类型安全优先（SFINAE with `std::enable_if_t`）
2. Header-only 优先（性能与易用性）
3. 编译期计算（constexpr when possible）
4. 清晰的 API 命名
5. 泛型编程（模板元编程）
