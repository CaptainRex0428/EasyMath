# 引用资料索引（References）

> 本文件列出项目可用的外部参考资料的 memory 索引。当任务需要相关知识时，**先查阅这里**——memory 中已抽取的概要足够回答绝大多数问题，**不要回读原始源文件**（如 epub、PDF 等），它们通常体积很大。

---

## 1. RTR4：实时渲染第四版（**核心 reference**）

**位置**：`.claude/memory/reference/RTR4/`
**源文件**：`reference/Real-Time Rendering, Fourth Edition (Tomas Akenine-Moller, Eric Haines, Naty Hoffman etc.).epub`
**作者**：Tomas Akenine-Möller, Eric Haines, Naty Hoffman, Angelo Pesce, Michal Iwanicki, Sébastien Hillaire
**中文版**：实时渲染第四版（pandoc 翻译）
**总章节数**：26 章

### 何时查阅
当任务涉及**计算机图形学 / 实时渲染**相关概念时，必须先查 RTR4：
- 矩阵 / 变换 / 四元数 API 设计 → Ch4
- 投影 / 视图矩阵推导 → Ch4 §4.7
- 颜色空间 / 光度学 / 亮度 → Ch8
- 物理着色（BRDF、菲涅尔、微表面） → Ch9
- 光照模型（Phong / Blinn-Phong / PBR Cook-Torrance） → Ch5 + Ch9 + Ch10
- 阴影算法 → Ch7
- 纹理过滤（mipmap、各向异性） → Ch6
- 加速结构（BVH、k-d 树、八叉树、LOD） → Ch19
- 碰撞检测（GJK、SAT、BVH CD） → Ch25
- 光线追踪（BVH 遍历、一致性） → Ch26
- 性能分析 → Ch18
- GPU 架构 → Ch23
- 相交测试（射线 / 包围体 / 平面） → Ch22
- 体积渲染 → Ch14

### 索引文件
- **`RTR4-overview.md`**：全书 26 章索引 + 与 EasyMath 项目代码的关联映射
- **`formulas.md`**：跨章节关键公式（LaTeX + 来源章节 + C++ 实现片段）
- **`principles.md`**：核心原理 / 算法思想（含 20 条原理）
- **`chapters/ch01-...md` ~ `ch26-...md`**：每章详细概要（1500-3000 字）

### 重点章节（与 EasyMath 强相关）
| 章节 | 标题 | 关联 EasyMath 头文件 |
|------|------|---------------------|
| Ch2 | 图形渲染管线 | `sandbox/sandbox.cpp` (MVP) |
| Ch4 | 变换 | `Matrix.h` + `Quaternion.h` + `EasyConversion.h` |
| Ch8 | 光与颜色 | `Color.h` |
| Ch9 | 基于物理的着色 | `Color.h`（未来 PBR） |
| Ch19 | 加速算法 | 未来 AABB / 剔除 API |
| Ch22 | 相交测试 | 未来求交 API |
| Ch25 | 碰撞检测 | 未来碰撞检测 API |
| Ch26 | 实时光线追踪 | 未来 ray tracing API |

---

## 2. 未来可添加的 reference

未来添加新 reference 时，请遵循以下结构：

```
.claude/memory/reference/<reference-name>/
├── <name>-overview.md      # 主概要（含章节索引）
├── chapters/               # 每章详细概要
│   └── chXX-name.md
├── formulas.md             # 关键公式（LaTeX + C++）
└── principles.md           # 核心原理
```

并在 `.claude/memory/references.md`（本文件）中追加新条目。

**建议格式**（每个 reference 条目）：
1. 位置（路径）
2. 源文件（项目内路径）
3. 作者/版本
4. 何时查阅（按任务领域）
5. 索引文件清单
6. 重点章节（如适用）
