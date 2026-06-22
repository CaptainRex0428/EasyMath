# EasyMath 项目记忆系统

## 📁 文件说明

### 1. project-context.md
**项目身份与核心信息**
- 项目名称、版本、开发者
- 项目定位与目标
- 技术栈与设计原则
- 项目结构概览

### 2. development-status.md
**当前开发状态**
- 各模块完成度（Vector 95%, Matrix 90%, Quaternion 30%, Color 85%）
- 已实现功能清单
- 待实现功能列表
- 最近提交记录
- 集成状态与文档状态

### 3. roadmap.md
**开发路线图**
- v1.1.0 详细计划（按优先级分类）
- 高优先级：Quaternion 完善、EulerAngles、旋转转换
- 中优先级：Color 运算、Matrix 投影
- 低优先级：SIMD 优化、测试、文档
- 未来版本规划（v1.2.0, v2.0.0）

### 4. code-conventions.md
**代码规范与最佳实践**
- 命名约定（类型、函数、常量、变量）
- 代码风格（缩进、括号、空格）
- 设计原则（类型安全、header-only、constexpr）
- 模板最佳实践（SFINAE、特化）
- UE 兼容性规范
- 版本控制规范

### 5. design-principles.md
**设计原则与核心机制**（⭐ 重要）
- 泛型设计模式（Vector、Matrix、Color 的模板参数设计）
- Swizzle 机制详解（Union 魔法 + 模板参数包）
- Color 类型安全系统（编译期防止错误转换）
- 颜色空间语义（sRGB vs Linear，何时使用）
- 伽马校正公式（标准 sRGB 转换）
- HSV/HSL/HSI 三种色彩模型对比
- 亮度计算标准（Rec601/709/2020）
- Matrix 线性代数实现（递归行列式、伴随矩阵法）
- 类型安全最佳实践（SFINAE、static_assert、constexpr）
- 性能优化策略（编译期优化、内存布局、零抽象开销）

---

## 🔄 使用方式

### 对于 AI 助手
这些文件帮助 AI 理解：
1. **项目是什么** - project-context.md
2. **进展到哪里** - development-status.md
3. **接下来做什么** - roadmap.md
4. **如何编写代码** - code-conventions.md
5. **为什么这样设计** - design-principles.md ⭐

### 对于开发者（CaptainRex）
1. 快速回顾项目状态
2. 规划下一步开发
3. 保持代码风格一致
4. 跟踪功能实现进度

---

## 📝 更新频率

### 每次开发后更新
- **development-status.md**: 更新完成度、提交记录
- **roadmap.md**: 标记已完成项，调整优先级

### 重大变更时更新
- **project-context.md**: 版本号、项目结构变化
- **code-conventions.md**: 新增规范、修改设计原则

---

## 🎯 快速索引

| 想要了解... | 查看文件 |
|------------|---------|
| 项目是做什么的 | project-context.md |
| 当前完成了哪些功能 | development-status.md |
| 下一步要做什么 | roadmap.md |
| 如何写代码 | code-conventions.md |
| **为什么这样设计** ⭐ | **design-principles.md** |
| Swizzle 是什么 | design-principles.md → Swizzle 机制 |
| Color 类型安全 | design-principles.md → Color 类型安全系统 |
| sRGB vs Linear | design-principles.md → 颜色空间语义 |
| Matrix 如何求逆 | design-principles.md → Matrix 线性代数 |
| Quaternion 进度 | development-status.md → Quaternion 部分 |
| 命名规范 | code-conventions.md → 命名约定 |
| UE 集成方式 | project-context.md → UE 集成机制 |
| **Demo 开发** | **demo-overview.md** |
| Matrix Demo 设计 | demo-matrix-design.md |
| Demo 开发技能 | ../skills/demo-dev.md |
