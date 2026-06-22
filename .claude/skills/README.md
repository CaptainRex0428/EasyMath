# Skills 目录

此目录用于存放 EasyMath 项目的自定义 Claude Code 技能。

## 可用技能

### 1. /demo-dev
**文件**: `demo-dev.md`
**功能**: EasyMath Demo 开发助手
- 创建新 Demo 页面
- 构建和预览
- 代码生成和检查
- 遵循 EasyMath 代码规范

## 规划中的技能

### 2. /math-test
**功能**: 快速运行数学库测试
- 编译 Sandbox 项目
- 运行测试用例
- 显示测试结果

### 3. /ue-build
**功能**: 构建 Unreal Engine 插件版本
- 切换到 UE 插件模式
- 生成 UE 插件项目文件
- 编译插件

### 4. /doc-gen
**功能**: 生成 API 文档
- 扫描头文件
- 提取函数签名和注释
- 生成 Markdown 文档

## 如何添加技能

技能文件应遵循 Claude Code 技能格式，包含：
- 技能名称
- 描述
- 参数
- 执行逻辑
- 模板配置

参考 `demo-dev.md` 的结构创建新技能。
