# Matrix Demo UI 重构记录

**日期**: 2026-06-22
**维护者**: CaptainRex
**状态**: ✅ 已完成

---

## 背景与动机

旧版 Matrix Demo 存在以下问题：
1. 两个独立页面（`webgl.astro` / `threejs.astro`），切换需要返回主页
2. 侧栏宽度固定 280px，无法调整
3. 使用 `height: 100vh` 导致小屏幕内容裁切
4. 无法折叠控制面板，侧栏内容拥挤
5. 缺少现代交互特性（状态持久化等）

---

## 改动清单

### 新增文件

| 文件 | 用途 |
|------|------|
| `demo/src/pages/matrix/index.astro` | 统一 Matrix Demo 主页面 |
| `demo/src/styles/global.css` | 全局 CSS 变量系统 |
| `demo/src/scripts/storage.js` | localStorage 安全封装 |
| `demo/src/scripts/resizable-sidebar.js` | 侧栏拖拽逻辑（独立模块） |
| `demo/src/components/ResizableSidebar.astro` | 可调整侧栏组件 |
| `demo/src/components/TabSwitcher.astro` | 标签页切换组件 |
| `demo/src/components/CollapsibleSection.astro` | 可折叠面板组件 |

### 修改文件

| 文件 | 改动 |
|------|------|
| `demo/src/layouts/Layout.astro` | 引入全局 CSS，更新 topbar 样式 |
| `demo/src/pages/index.astro` | 更新链接指向 `/matrix/` |
| `demo/src/pages/matrix/webgl.astro` | 改为重定向到 `/matrix/?tab=webgl` |
| `demo/src/pages/matrix/threejs.astro` | 改为重定向到 `/matrix/?tab=threejs` |
| `.claude/memory/demo-overview.md` | 更新项目结构和 UI 规范 |
| `.claude/DEMO_INTEGRATION.md` | 修正绝对路径为相对路径 |
| `.claude/settings.local.json` | 替换硬编码绝对路径权限为通用规则 |

---

## 核心架构决策

### 单页面 + 标签页切换

将 WebGL 和 Three.js 合并到 `/matrix/index.astro`，通过 `?tab=webgl|threejs` URL 参数切换。两个渲染引擎共享：
- 同一套控制面板滑块
- 同一套矩阵数学函数
- 同一个矩阵显示面板
- 同一个管道可视化条

### CSS Grid 布局

```
grid-template-columns: var(--sidebar-width) 4px 1fr
grid-template-rows: 44px minmax(300px,1fr) 58px 190px
```

第二列（4px）是拖拽手柄 `#sidebar-resizer`。

### localStorage 键名约定

所有键名使用 `em-` 前缀：

| 键 | 值类型 | 默认 |
|----|--------|------|
| `em-active-tab` | `'webgl'` \| `'threejs'` | `'webgl'` |
| `em-sidebar-width` | number (px) | 280 |
| `em-proj-mode` | `'persp'` \| `'ortho'` | `'persp'` |
| `em-section-model-trans` | boolean | true |
| `em-section-model-rot` | boolean | true |
| `em-section-model-scale` | boolean | true |
| `em-section-camera` | boolean | true |
| `em-section-projection` | boolean | true |

### 向后兼容

旧路由通过 Astro `Astro.redirect()` 301 重定向：
- `/matrix/webgl/` → `/matrix/?tab=webgl`
- `/matrix/threejs/` → `/matrix/?tab=threejs`

---

## 验证

构建输出（`npx astro build`）：
```
4 page(s) built in 457ms — Complete!
```

- `/index.html` ✅
- `/matrix/index.html` ✅
- `/matrix/webgl/index.html` ✅ (redirect)
- `/matrix/threejs/index.html` ✅ (redirect)

---

## 已知限制

- Three.js 通过 CDN 加载（`cdn.jsdelivr.net/npm/three@0.158.0`），离线环境无法使用 Three.js 标签页，但 WebGL 标签页不受影响
- `localStorage` 在隐私模式下不可用，所有持久化降级为默认值（try-catch 保护）
- 侧栏拖拽仅支持鼠标（桌面端）；移动端侧栏改为汉堡菜单叠加层模式
