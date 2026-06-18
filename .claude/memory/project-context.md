# EasyMath 项目上下文

## 项目身份
- **项目名称**: EasyMath
- **当前版本**: 1.0.11
- **目标版本**: 1.1.0
- **开发者**: CaptainRex
- **项目类型**: C++17 实时渲染数学库
- **命名空间**: `EM`

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
