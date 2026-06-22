# EasyMath Demo 开发规范总结

## 概述

本文档总结了 EasyMath Demo 部分的开发规范、设计原则和组织结构，作为整个项目 Demo 开发的权威指南。

## 📋 文档结构

### Memory 文档（记录知识）

1. **demo-overview.md** - Demo 项目总览
   - 技术架构选择（Astro 框架）
   - 项目结构和组织
   - 已实现 Demo 清单
   - 开发规范和最佳实践

2. **demo-matrix-design.md** - Matrix Demo 详细设计
   - 变换管道架构
   - 矩阵实现细节
   - WebGL 实现要点
   - UI 和交互设计
   - 故障排除指南

### Skills 文档（定义技能）

1. **demo-dev.md** - Demo 开发技能
   - 技能定义和接口
   - 自动化模板
   - 最佳实践指南
   - 故障排除

## 🎯 核心原则

### 1. 与 C++ 源码一一对应

所有数学函数必须与 EasyMath C++ 实现对应：

```javascript
// C++ 版本
template<typename T>
Matrix<T, 4, 4> MTXTranslation(T x, T y, T z);

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

### 2. 自实现而非依赖库

- ❌ 不使用 gl-matrix、Three.js 内部矩阵
- ✅ 完全自实现所有数学运算
- ✅ 可实时显示计算过程
- ✅ 便于理解和学习

### 3. 统一的 UI 设计

**深色主题**（与 EasyMath 主项目一致）:
```css
--bg-color: #0d1117;
--primary-color: #58a6ff;
--border-color: #21262d;
--card-bg: #161b22;
--text-color: #e6edf3;
```

**响应式布局**:
```css
body { overflow: auto; }
#app { min-height: 100vh; }
canvas { min-height: 200px; }
.info-panel { overflow-x: auto; }
```

### 4. 清晰的信息展示

每个 Demo 应显示：
- **使用的函数** - 例如 `MTXTranslation(tx, ty, tz)`
- **实时数值** - 矩阵/向量当前值
- **变换管道** - 空间变换链
- **标注说明** - 关键概念标注

## 🔧 技术栈

### Astro 框架

**为什么选择 Astro**:
- ✅ 静态站点编译，部署简单
- ✅ 组件化设计，易于维护
- ✅ 零配置路由管理
- ✅ 优秀的开发体验（热重载）
- ✅ 易于扩展

### 为什么不是纯 HTML？

| 特性 | 纯 HTML | Astro |
|-----|---------|-------|
| 模块化 | ❌ | ✅ |
| 路由管理 | ❌ | ✅ |
| 代码复用 | ❌ | ✅ |
| 构建优化 | ❌ | ✅ |
| 开发体验 | ⚠️ | ✅ |

## 📁 目录结构

```
demo/
├── src/
│   ├── layouts/
│   │   └── Layout.astro          # 基础布局（支持 head slot）
│   └── pages/
│       ├── index.astro          # 主页（Demo 画廊）
│       └── matrix/
│           ├── webgl.astro      # Pure WebGL 版本
│           └── threejs.astro    # Three.js 版本
├── astro.config.mjs              # Astro 配置
├── package.json                  # 依赖管理
└── dist/                         # 构建输出
```

## 🚀 开发工作流

### 1. 创建新 Demo

```bash
# 方式 1: 手动创建
cd demo/src/pages
mkdir -p new-demo
touch new-demo/index.astro

# 方式 2: 使用技能（规划中）
/demo-dev create --type matrix --name "new-demo"
```

### 2. Demo 页面模板

```astro
---
import Layout from '../../layouts/Layout.astro';
---
<Layout title="<Demo Name>" back={{ href: '/', label: '返回主页' }}>
<div id="app">
  <div id="controls">
    <!-- 控制面板 -->
  </div>
  <canvas id="glcanvas"></canvas>
  <div id="info-panel">
    <!-- 信息显示 -->
  </div>
</div>
</Layout>

<style>
/* 样式 */
</style>

<script is:inline>
// 内联脚本（避免打包）
</script>
```

### 3. 开发和测试

```bash
# 开发模式
npx astro dev --port 4321

# 构建
npx astro build

# 预览
npx astro preview
```

## 📐 矩阵显示规范

### 格式化

```javascript
function fmtMatrix(m) {
  let out = '';
  for (let i = 0; i < 4; i++) {
    out += '| ';
    for (let j = 0; j < 4; j++) {
      const val = m[i * 4 + j];
      out += (val >= 0 ? ' ' : '') + val.toFixed(2) + ' ';
    }
    out += '|\n';
  }
  return out;
}
```

### 标注

```
┌─ Model Matrix ─────────────┐
│ MTXTranslation(tx, ty, tz)  │
│ ┌────────────────────────┐ │
│ │ 1.00  0.00  0.00  0.00│ │
│ │ 0.00  1.00  0.00  0.00│ │
│ │ 0.00  0.00  1.00  0.00│ │
│ │ 0.00  0.00  0.00  1.00│ │
│ └────────────────────────┘ │
└────────────────────────────┘
```

## 🎨 UI 组件规范

### 控制面板

```html
<div class="section">
  <div class="section-title">Model — 平移 Translation</div>
  <div class="row">
    <label>tx</label>
    <input type="range" id="tx" min="-3" max="3" step="0.1" value="0">
    <span id="tx-v">0.0</span>
  </div>
</div>
```

### 管道可视化

```
[Local] → [World] → [Camera] → [Clip] → [NDC]
   M         Model      View        Proj
```

### 切换按钮

```html
<div class="toggle-row">
  <button class="toggle-btn active" id="btn-persp" onclick="setProj('persp')">
    Perspective
  </button>
  <button class="toggle-btn" id="btn-ortho" onclick="setProj('ortho')">
    Orthographic
  </button>
</div>
```

## ⚠️ 常见陷阱

### 1. WebGL 矩阵传递

```javascript
// ❌ 错误：未转置
gl.uniformMatrix4fv(uMVP, false, MVP.data);

// ✅ 正确：转置
gl.uniformMatrix4fv(uMVP, false, transpose(MVP.data));
```

### 2. 矩阵乘法顺序

```javascript
// ❌ 错误：从左到右
const MVP = M.multiply(V).multiply(P);

// ✅ 正确：从右到左
const MVP = P.multiply(V).multiply(M);
```

### 3. 网格渲染空间

```javascript
// ❌ 错误：网格跟着立方体转
drawGrid(MVP);

// ✅ 正确：网格固定在世界空间
drawGrid(VP);
```

### 4. 响应式设计

```css
/* ❌ 错误：小屏幕隐藏内容 */
body { overflow: hidden; }
#app { height: 100vh; }

/* ✅ 正确：保持可见性 */
body { overflow: auto; }
#app { min-height: 100vh; }
```

## 📊 性能优化

### 1. 避免内存分配

```javascript
// ❌ 每帧分配
function render() {
  const MVP = new Matrix(16);
}

// ✅ 复用数组
const MVP = new Matrix(16);
function render() {
  MVP.update();
}
```

### 2. 缓存计算结果

```javascript
// 缓存不变的矩阵
let cachedP = null;
function getProjection() {
  if (!cachedP || projChanged) {
    cachedP = MTXPerspective(fov, aspect, near, far);
  }
  return cachedP;
}
```

### 3. 减少矩阵乘法

```javascript
// 复用转置结果
const MVP = P.multiply(V).multiply(M);
const transposed = transpose(MVP);
gl.uniformMatrix4fv(uMVP, false, transposed);
gl.uniformMatrix4fv(uMVP2, false, transposed);
```

## 🎯 质量检查清单

在提交新 Demo 前，检查：

- [ ] 所有数学函数与 C++ 源码对应
- [ ] 使用深色主题颜色变量
- [ ] 响应式布局正常工作
- [ ] 矩阵显示格式正确
- [ ] 管道可视化清晰
- [ ] 函数名标注完整
- [ ] 无控制台错误
- [ ] 在不同浏览器测试通过
- [ ] 性能可接受（60 FPS）

## 🔗 相关资源

### C++ 源码
- `include/Matrix.h` - 矩阵实现
- `include/Vector.h` - 向量运算
- `sandbox/sandbox.cpp` - C++ 测试

### Demo 文件
- `demo/src/pages/matrix/webgl.astro` - Pure WebGL 版本
- `demo/src/pages/matrix/threejs.astro` - Three.js 版本

### 文档
- `docs/Matrix_API.md` - Matrix API 文档
- `docs/Matrix_CN.md` - Matrix 中文文档

## 📈 未来规划

### 短期（进行中）
- [ ] Quaternion Demo
- [ ] Vector Demo
- [ ] Color Demo

### 中期（规划中）
- [ ] 光照计算 Demo
- [ ] 贝塞尔曲线 Demo
- [ ] 性能基准测试

### 长期（设想）
- [ ] 交互式教程模式
- [ ] 代码生成器
- [ ] 在线编辑器

---

**维护者**: CaptainRex
**最后更新**: 2026-06-22
**版本**: 1.0.0
