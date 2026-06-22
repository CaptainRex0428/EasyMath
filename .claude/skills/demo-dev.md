# Demo 开发技能

**技能名称**: `/demo-dev`
**描述**: EasyMath Demo 开发助手

## 功能

此技能帮助开发者快速创建、测试和部署 EasyMath Demo 页面。

## 使用方法

### 创建新 Demo

```bash
/demo-dev create --type <matrix|vector|color|quaternion> --name <demo-name>
```

**示例**:
```bash
/demo-dev create --type vector --name dot-product
# 创建 src/pages/vector/dot-product.astro
```

### 构建和预览

```bash
/demo-dev build        # 构建静态站点
/demo-dev preview     # 预览构建结果
/demo-dev dev         # 启动开发服务器
```

### 测试 Demo

```bash
/demo-dev test --name <demo-name>
# 在浏览器中打开指定 demo
```

## 技能模板

### Matrix Demo 模板

当创建 Matrix 类型的 Demo 时，技能会自动生成：

```astro
---
import Layout from '../../layouts/Layout.astro';
---
<Layout title="<Demo Name> - Matrix Demo" back={{ href: '/', label: '返回主页' }}>
<div id="app">
  <div id="controls">
    <!-- 控制面板 -->
  </div>
  <canvas id="glcanvas"></canvas>
  <div id="matrix-panel">
    <!-- 矩阵显示 -->
  </div>
</div>
</Layout>

<style>
/* 复用 demo 样式 */
</style>

<script is:inline>
// 自实现矩阵数学
class Matrix {
  constructor(rows, cols) { /* ... */ }
  // ... 对应 EasyMath C++ 实现
}

// WebGL 渲染
// ... 渲染代码
</script>
```

### Vector Demo 模板

```astro
---
import Layout from '../../layouts/Layout.astro';
---
<Layout title="<Demo Name> - Vector Demo" back={{ href: '/', label: '返回主页' }}>
<div id="app">
  <div id="controls">
    <!-- 向量控制 -->
  </div>
  <canvas id="canvas"></canvas>
  <div id="vector-panel">
    <!-- 向量显示 -->
  </div>
</div>
</Layout>
```

## 自动化任务

### 1. 添加主页入口

创建新 Demo 后，技能会自动更新 `src/pages/index.astro`，添加新 demo 的入口。

### 2. 生成文档

自动生成 demo 说明文档，包括：
- 功能描述
- 对应的 EasyMath C++ 函数
- 数学原理
- 使用说明

### 3. 代码检查

检查新代码是否符合：
- EasyMath 代码规范
- Demo UI 设计原则
- 响应式设计要求

## 最佳实践

### 1. 对应 C++ 源码

所有数学函数应与 EasyMath C++ 一一对应：

```javascript
// C++ 版本
template<typename T>
Matrix<T, 4, 4> MTXTranslation(T x, T y, T z)
{
    return {1,0,0,x, 0,1,0,y, 0,0,1,z, 0,0,0,1};
}

// JavaScript 版本（完全对应）
function MTXTranslation(x, y, z) {
    return new Matrix([
        1, 0, 0, x,
        0, 1, 0, y,
        0, 0, 1, z,
        0, 0, 0, 1
    ]);
}
```

### 2. 矩阵显示规范

```javascript
function displayMatrix(name, matrix, functionName) {
    // 格式：名称
    //      函数名
    //      数值（保留2位小数）
    console.log(`${name}`);
    console.log(`${functionName}`);
    console.log(fmtMatrix(matrix));
}
```

### 3. 管道可视化

对于涉及变换管道的 Demo，应显示：

```
[Local] → [World] → [Camera] → [Clip] → [NDC]
   M         Model      View        Proj
```

### 4. 响应式设计

确保在所有屏幕尺寸下可用：

```css
body { overflow: auto; }
#app { min-height: 100vh; }
canvas { min-height: 200px; }
.info-panel { overflow-x: auto; }
```

## 故障排除

### 问题：开发服务器无法启动

```bash
# 清理缓存
rm -rf .astro/
rm -rf dist/

# 重新安装依赖
rm -rf node_modules/
npm install

# 重新启动
npx astro dev
```

### 问题：矩阵显示不正确

检查：
1. 矩阵存储顺序（行主序 vs 列主序）
2. WebGL uniform 是否正确转置
3. 矩阵乘法顺序是否正确

### 问题：UI 在小屏幕上被隐藏

检查：
1. 是否使用了 `overflow: hidden`
2. 是否使用了固定 `height: 100vh`
3. 是否给滚动容器添加了 `overflow-x: auto`

## 技能配置

在 `.claude/skills/demo-dev.md` 中配置：

```yaml
name: demo-dev
version: 1.0.0
author: CaptainRex
description: EasyMath Demo 开发助手

templates:
  matrix: templates/matrix.astro
  vector: templates/vector.astro
  color: templates/color.astro
  quaternion: templates/quaternion.astro

rules:
  - 遵循 EasyMath 代码规范
  - 对应 C++ 源码实现
  - 深色主题
  - 响应式设计
```

## 示例工作流

```bash
# 1. 创建新的矩阵 Demo
/demo-dev create --type matrix --name "transform-chain"

# 2. 技能自动生成文件
# Created: src/pages/matrix/transform-chain.astro

# 3. 编辑 Demo 内容
# 添加控制面板、渲染逻辑、可视化

# 4. 测试 Demo
/demo-dev test --name transform-chain

# 5. 构建静态站点
/demo-dev build

# 6. 预览构建结果
/demo-dev preview
```

## 注意事项

1. **自实现数学函数** - 不要依赖外部数学库
2. **行主序存储** - 与 C++ 保持一致
3. **WebGL 转置** - 传递矩阵给 shader 时需要转置
4. **内联脚本** - 使用 `<script is:inline>` 避免打包问题
5. **深色主题** - 使用指定的颜色变量
6. **响应式** - 确保在所有屏幕尺寸下可用

---

## 技能状态

**当前版本**: 1.0.0
**最后更新**: 2026-06-22
**维护状态**: 活跃开发中
