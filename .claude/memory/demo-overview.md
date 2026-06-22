# EasyMath Demo 概述

**最后更新**: 2026-06-22
**维护者**: CaptainRex

## Demo 定位

EasyMath Demo 是 EasyMath C++ 数学库的**交互式可视化演示平台**，旨在帮助开发者理解数学概念在实际渲染中的应用。

### 核心目标

1. **可视化数学原理** - 将抽象的数学概念通过实时渲染展示
2. **对应 C++ 实现** - Web 版本与 EasyMath C++ 源码一一对应
3. **教育价值** - 帮助学习者理解图形学数学基础
4. **技术展示** - 展示不同技术栈的实现方式（Pure WebGL vs Three.js）

---

## 技术架构

### 框架选择

**Astro** - 静态站点生成器
- ✅ 零配置静态编译
- ✅ 原子化组件设计
- ✅ 支持 React/Vue/Svelte 等框架集成
- ✅ 优秀的 SEO 和性能
- ✅ 易于部署（可部署到 GitHub Pages、Netlify 等）

### 为什么使用 Astro 而不是纯 HTML？

1. **模块化** - 组件可复用，易于维护
2. **路由管理** - 自动生成路由结构
3. **构建优化** - 自动压缩、优化资源
4. **扩展性** - 未来可轻松添加新的 demo 案例
5. **开发体验** - 热重载、TypeScript 支持、开发服务器

---

## 项目结构

```
demo/
├── src/
│   ├── components/               # 可复用组件
│   │   ├── CollapsibleSection.astro  # 可折叠面板
│   │   ├── ResizableSidebar.astro    # 可调整侧栏
│   │   └── TabSwitcher.astro         # 标签页切换
│   ├── layouts/
│   │   └── Layout.astro          # 基础布局（引入全局 CSS）
│   ├── pages/
│   │   ├── index.astro           # 主页（Demo 画廊）
│   │   └── matrix/
│   │       ├── index.astro       # 统一 Matrix Demo（WebGL + Three.js）
│   │       ├── webgl.astro       # 重定向 → /matrix/?tab=webgl
│   │       └── threejs.astro     # 重定向 → /matrix/?tab=threejs
│   ├── scripts/
│   │   ├── resizable-sidebar.js  # 侧栏拖拽逻辑
│   │   └── storage.js            # localStorage 工具
│   └── styles/
│       └── global.css            # 全局 CSS 变量（GitHub Dark 主题）
├── dist/                         # 编译输出
├── astro.config.mjs
└── package.json
```

### 历史遗留

> ⚠️ `demo/matrix-visualizer/`（早期独立 HTML 版本）已删除，请使用 Astro 版本进行开发

---

## 已实现 Demo

### 1. Matrix 矩阵演示

**统一 Matrix Demo** (`/matrix/`) — 2026-06-22 重构为单页面架构
- ✅ WebGL / Three.js 标签页切换（URL 参数 `?tab=webgl` / `?tab=threejs`）
- ✅ 可调整侧栏（200px–600px，localStorage 持久化）
- ✅ 可折叠控制面板（平移/旋转/缩放/相机/投影）
- ✅ 完全自实现矩阵数学（与 EasyMath C++ 一一对应）
- ✅ 管道可视化（Local → World → Camera → Clip → NDC）
- ✅ 实时矩阵数值显示（4 个矩阵：Model / View / Proj / MVP）
- ✅ 响应式布局（min-height 而非固定 height）
- ✅ 移动端适配（汉堡菜单 + 侧栏叠加层）
- ✅ 所有状态 localStorage 持久化

**旧链接重定向**（向后兼容）:
- `/matrix/webgl/` → `/matrix/?tab=webgl`
- `/matrix/threejs/` → `/matrix/?tab=threejs`

**对应 EasyMath C++ 源码**：
- `Matrix.h` - MTXTranslation, MTXRotationX/Y/Z, MTXScale, MTXLookAt
- `Matrix.h` - MTXPerspective, MTXOrtho
- `Vector.h` - 向量运算（点积、叉积、归一化）

---

## 设计规范

### 1. 页面结构

每个 Demo 页面应包含：

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
      <!-- 实时信息显示 -->
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

### 2. UI 设计原则

**深色主题** — GitHub Dark 调色板，通过 CSS 变量统一管理 (`src/styles/global.css`)
```css
--bg-primary: #0d1117;    --bg-secondary: #161b22;  --bg-tertiary: #21262d;
--border-color: #30363d;  --text-primary: #e6edf3;  --text-secondary: #8b949e;
--accent-blue: #58a6ff;   --accent-green: #3fb950;  --accent-orange: #f0883e;
--sidebar-width: 280px;   /* 动态，由 JS 更新 */
```

**响应式布局**
```css
#app { min-height: 100vh; }   /* 不用 height: 100vh */
```

**可调整侧栏**: 200px–600px，鼠标拖拽 `#sidebar-resizer`，宽度持久化到 localStorage

**可折叠面板**: `<details>/<summary>` 原生实现，状态持久化到 localStorage

**标签页切换**: URL 参数 `?tab=webgl|threejs`，状态持久化到 localStorage

### 3. 矩阵显示规范

**格式**
```javascript
function fmtMatrix(m) {
  // 4x4 矩阵，每行 4 个元素
  // 保留 2 位小数
  // 对齐显示
}
```

**标注函数名**
- 明确标注每个矩阵使用的 EasyMath 函数
- 例如：`MTXTranslation(tx, ty, tz)`
- 例如：`MTXLookAt(eye, target, up)`

**管道可视化**
```
[Local] → [World] → [Camera] → [Clip] → [NDC]
   M         Model      View        Proj
```

### 4. 代码风格

**自实现数学函数**
- 不使用外部数学库（gl-matrix、Three.js 内部）
- 完全对应 EasyMath C++ 实现
- 行主序存储（Row-major）
- 与 C++ 一致的函数签名

**WebGL 适配**
- WebGL 使用列主序，需要转置
- 使用 `gl.uniformMatrix4fv(loc, false, transpose(m))`
- `false` = 不转置（我们手动转置）

---

## 命令速查

### 开发

```bash
cd demo
npm install                 # 安装依赖
npx astro dev --port 4321   # 启动开发服务器
```

### 构建

```bash
npx astro build             # 构建静态站点
npx astro preview           # 预览构建结果
```

### 清理

```bash
rm -rf dist/                # 删除构建输出
rm -rf node_modules/        # 删除依赖
```

---

## 未来规划

### v1.0 Demo 路线图

**短期**（已完成）:
- ✅ Matrix 基础变换演示
- ✅ Pure WebGL 版本
- ✅ Three.js 版本
- ✅ Astro 框架集成

**中期**（进行中）:
- [ ] Quaternion 四元数演示
- [ ] Color 颜色空间演示
- [ ] Vector 向量运算演示
- [ ] 投影矩阵对比演示

**长期**（规划中）:
- [ ] 光照计算可视化
- [ ] 贝塞尔曲线演示
- [ ] 性能基准测试
- [ ] 交互式教程模式

---

## 开发规范

### 新增 Demo 步骤

1. **规划 Demo 内容**
   - 确定要演示的数学概念
   - 设计交互方式
   - 规划 UI 布局

2. **创建 Astro 页面**
   - 在 `src/pages/` 下创建新目录和 `.astro` 文件
   - 使用 Layout 组件
   - 遵循命名规范

3. **实现核心逻辑**
   - 自实现数学函数（对应 C++）
   - WebGL/Canvas 渲染代码
   - 交互控制逻辑

4. **添加可视化**
   - 实时数值显示
   - 矩阵/向量显示
   - 管道/图解可视化

5. **测试和优化**
   - 测试不同浏览器
   - 测试不同屏幕尺寸
   - 性能优化

### 文档要求

每个 Demo 应包含：
- 功能说明
- 对应的 EasyMath C++ 函数
- 数学原理说明
- 使用说明

---

## 常见问题

### Q: 为什么不直接使用 Three.js 的矩阵？

A: Demo 的目的是展示 EasyMath C++ 的实现原理。Three.js 的矩阵是内部实现，无法看到计算过程。我们自实现矩阵可以：
1. 与 C++ 源码一一对应
2. 实时显示中间结果
3. 理解矩阵组合顺序
4. 学习底层数学原理

### Q: WebGL 和 Three.js 版本的区别？

A:
- **Pure WebGL**: 使用原生 WebGL API，完全手写 shader 和矩阵运算
- **Three.js**: 使用 Three.js 渲染，但矩阵运算仍是自实现（显示用）

### Q: 为什么使用 Astro？

A:
- 静态站点，部署简单
- 组件化，易于维护
- 扩展性强，未来可添加更多 demo
- 开发体验好，热重载快

### Q: 如何添加新的 Demo？

A:
1. 在 `src/pages/` 下创建新文件
2. 复用 Layout 组件
3. 实现渲染和交互逻辑
4. 在主页添加入口

---

## 贡献指南

如果你有新的 Demo 想法：
1. 提出 Issue 描述你的想法
2. 设计 Demo 的交互方式
3. 实现并提交 PR
4. 确保遵循现有代码风格和 UI 规范

---

## 许可

遵循 EasyMath 主项目的许可证。
