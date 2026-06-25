# Reference 提取标准流程（SOP）

> **目的**：未来添加新的 reference（教材、论文集、规范）时，遵循统一流程，产出可复用的结构化 memory。
> **核心原则**：**memory 是快速入口，源文件是真相**——内容不足或冲突时**必须回读源文件**。

---

## 0. 适用场景

当项目有新的外部参考（epub、pdf、网页、规范）需要抽取到 memory 时。

**当前已抽取的 reference**：
- `reference/RTR4/`（Real-Time Rendering 4th Edition，26 章）

**未来可能添加**：
- Physically Based Rendering (PBRT)
- Computer Graphics: Principles and Practice
- Game Engine Architecture
- OpenGL/Vulkan/DirectX 规范
- IEEE 浮点标准
- 各类 SIGGRAPH 论文集

---

## 1. 整体流程

```
[识别需求] → [创建目录骨架] → [预处理源文件] → [抽取概要] → [写公式/原理] → [建索引] → [验证]
```

每一步的细节见后续章节。

---

## 2. 步骤 1：识别需求

**问题清单**（开工前问自己）：
- Q1. 这本 reference 在项目中的角色是什么？（核心？补充？规范？）
- Q2. 哪些章节与本项目当前/未来代码强相关？
- Q3. 是否需要"详尽"（每章 1500-3000 字）还是"骨架"（每章 200-400 字）？

**决策依据**：
- 核心参考（如 RTR4 这种"百科全书"）→ 详尽
- 补充参考（如单一主题论文集）→ 骨架
- 规范（如 OpenGL 规范）→ 只抽与本项目相关的部分

---

## 3. 步骤 2：创建目录骨架

每个 reference 使用统一结构：

```
.claude/memory/reference/<short-name>/
├── <name>-overview.md          # 主概要（章节索引 + 与项目代码的关联映射）
├── formulas.md                  # 跨章关键公式（LaTeX + 来源 + C++ 片段）
├── principles.md                # 核心原理/算法思想
├── chapters/                    # 每章详细概要
│   └── chXX-name.md
└── _extraction-log.md           # 提取日志（本次提取做了什么、跳过了什么、为什么）
```

**注意**：**不保留 `_raw/` 预处理文件**。预处理（如去噪的 .txt）会丢失 LaTeX 公式，而公式是核心内容。**未来回读直接读原始 .xhtml/.html**（用 offset/limit 减少 token）。

**`_extraction-log.md` 很重要**——记录哪些章节因为"与项目无关"被跳过、哪些是"低优先级待补"。

**模板**：

```bash
mkdir -p .claude/memory/reference/<name>/chapters
touch .claude/memory/reference/<name>/_extraction-log.md
```

---

## 4. 步骤 3：预处理源文件（**保留公式**）

**目的**：去掉 EPUB/PDF 的样式噪音，**必须保留 LaTeX 公式**。

**核心原则**：**不存储预处理后的 .txt/.md 文件**到 memory——预处理过程要么丢公式，要么文件大。未来回读源文件时，**直接读原始 .xhtml/.html**（用 offset/limit 控制读取量）。

### 4.1 抽取时使用 pandoc（**推荐**）

```bash
# 单文件转换（保留 MathML，pandoc 默认支持）
pandoc -f epub -t markdown --wrap=none \
       -o /tmp/book_full.md \
       "path/to/book.epub"

# 拆分章节
csplit -f /tmp/ch_ -b "%02d.md" /tmp/book_full.md '/^# /' '{*}' 2>/dev/null
```

**优点**：保留 `$...$` 和 `$$...$$` 公式。
**缺点**：在 /tmp 临时存放，**不要 commit 到 memory**。

### 4.2 不存储预处理文件

**为什么**：
- 预处理文件通常 50-200KB/章
- 抽取完成后没用了
- memory 应保持精简

**回读策略**：未来遇到需要源文件时，**直接 Read 原始 .xhtml/.html**：
- 用 `Grep` 精确定位关键词
- 用 `Read offset=N limit=M` 控制读取段

### 4.3 EPUB 直接处理

EPUB 章节文件 (.xhtml/.html) 通常包含：
- 大量 MathJax SVG 噪音（占 80%+ 内容）
- 关键正文段落
- 标题、列表

**读取技巧**：
- `Read` 工具带 `offset` 和 `limit` 跳过 SVG 块
- 优先读 `nav.xhtml` 找章节定位
- 用 `Grep` 找公式关键词（如 "Fresnel"、"Schlick"）

### 4.4 PDF → Markdown

```bash
# 选项 1：marker（开源，质量好）
pip install marker-pdf
marker_single "path/to/book.pdf" --output_dir /tmp/

# 选项 2：mathpix（付费，质量最好）
mathpix --format markdown "path/to/book.pdf" > /tmp/book.md

# 选项 3：pdftotext（无公式）
pdftotext -layout "path/to/book.pdf" /tmp/book.txt
```

**警告**：所有 PDF 工具输出文件**不存储**到 memory，仅作临时读取。

---

## 5. 步骤 4：抽取概要

### 5.1 抽取顺序

**不要一次性抽取所有章节**。推荐顺序：

1. **先读完整目录**（如 `nav.xhtml`）→ 建立全局图
2. **先写主概要**（`<name>-overview.md`）→ 锁定哪些章节相关
3. **按优先级抽取**：
   - **P0（必抽）**：核心参考章节、与项目代码强相关
   - **P1（应抽）**：重要概念章节
   - **P2（可跳）**：与项目无关的章节（如历史回顾、未来展望）

### 5.2 每章概要模板

每个 `chXX-name.md` 应包含以下结构：

```markdown
# ChXX <标题>

> **与项目关联**：<一句话说明这个章节与项目代码的关系>
> **核心问题**：<这个章节要回答的核心问题>
> **完成度**：<P0/P1/P2 + 是否完整抽取>

## 核心概念
<用一段话概述章节核心>

## 关键算法/思想
<用列表列出章节的主要算法/思想，每条 1-2 句话>

## 关键公式
<LaTeX 块 + 来源小节>
**重要**：公式来自 `_raw/chXX.md` 的具体小节（不是凭印象写）

## 与项目代码的关联
<哪些头文件/类/方法与本章相关，列出已实现 vs 待实现>

## 与其他章节的关联
<本章引用了哪些其他章节，被哪些章节引用>

## 关键术语索引
| 术语 | 释义 |

## 待补
<本章中尚未抽取的细节，列出优先级>

## 抽取历史
- 日期、来源、置信度
```

### 5.3 公式约定（**必须严格**）

- **LaTeX 块**：`$$...$$` 或 `$...$`
- **来源标注**：每条公式后注明 `来源：ChX §X.X` 或 `来源：本节`
- **C++ 片段**：用 EasyMath 风格的 `namespace EasyMath`（即使代码不一定能编译）
- **不复制源码**——给出关键行 + 解释
- **公式必须有数学正确性**——遇到模糊公式必须**回读源文件**

### 5.4 原理/算法抽取

- **不复制原书**——压缩为关键步骤 + 直觉
- **包含对比表格**（如 X vs Y 的优缺点）
- **包含决策树**（什么时候用哪个）

---

## 6. 步骤 5：跨章汇总

`formulas.md` 和 `principles.md` 是**索引式**文件：

- 公式按主题分组（变换/着色/纹理/碰撞/光线追踪）
- 每条公式标注来源章节
- 原理按"问题→方案→权衡"组织

**目的**：未来查找时，**不必打开每个 chapter 文件**，先看汇总。

---

## 7. 步骤 6：建索引

更新 `.claude/memory/references.md`，添加新条目。

**条目格式**：
1. 位置（路径）
2. 源文件（项目内路径）
3. 作者/版本
4. 何时查阅（按任务领域）
5. 索引文件清单
6. 重点章节

---

## 8. 步骤 7：验证与日志

### 8.1 写 `_extraction-log.md`

记录：
- 抽取日期、抽取人（AI/人类）
- 抽取范围（P0/P1/P2 各抽了多少章）
- 跳过的章节及原因
- 已知不完整/低置信度的章节
- 未来回读源文件的优先顺序

### 8.2 质量自检清单

- [ ] 公式有来源标注？
- [ ] C++ 片段有 namespace？
- [ ] 章节概要包含"与项目关联"段落？
- [ ] 关键术语有索引？
- [ ] `references.md` 已更新？
- [ ] `_extraction-log.md` 已写？
- [ ] **公式正确性**（关键）—— 不确定时**必须回读源文件**

---

## 9. 关键原则

### 9.1 Memory 是快速入口，源文件是真相

**不要**试图把 reference 的所有内容复制到 memory。

**memory 应承担**：
- 概念解释、定义、关键术语
- **关键公式（保留 LaTeX）**、算法直觉
- 章节地图、"看哪一章"
- 与项目代码的关联

**memory 不应承担**：
- 完整长公式推导（保留骨架，详细回读源文件）
- 完整代码示例
- 大段文字照搬
- 习题、边角案例

### 9.2 何时回读源文件

**触发条件**（任一即应回读）：
1. ✅ memory 中没有该概念
2. ✅ memory 中有该概念但与用户问题**不完全匹配**
3. ✅ memory 中内容与**其他参考资料**（论文、规范）**冲突**
4. ✅ 用户问的是**精确数值**（如某个常数、某个 RGB 值）
5. ✅ 用户问的是**完整推导**或**实现细节**
6. ✅ memory 标注"待补"或"低置信度"
7. ✅ 涉及**数学公式细节**（用户问参数、系数、推导步骤）

**回读策略**（按优先级）：
1. 打开 `<name>-overview.md` → 定位章节
2. 打开对应 `chapters/chXX-name.md` → 检查是否有现成概要
3. 如果内容不足/冲突 → 读**原始 EPUB 的 .xhtml/.html**（用 `Read offset=N limit=M`）
4. 用 `Grep` 在源文件中精确定位关键词（如"Schlick"、"fresnel"）
5. 仅读必要段落（不读整章）

**重要**：
- **memory 是快速入口，不是源文件副本**
- 公式若不准确 → 读源��件原始 LaTeX
- 源文件 MathJax SVG 多 → 用 offset/limit 分段

### 9.3 不要"过度抽取"

- 与项目**无关**的章节（P2）→ 只在 overview 里写一行话，不写详细 md
- 与项目**间接相关**的章节 → 在 overview 里写"与 ChX 关联"
- 重复内容 → 在主概要里 cross-reference，不重复写

### 9.4 复用而非重写

未来添加新 reference 时：
- **先看 `RTR4/` 的结构**——它是最完整的样板
- **复用其模板**（章节 md 结构、公式格式、原理分类）
- **复用其约定**（术语索引、与项目映射的写法）

---

## 10. 工具与脚本

| 工具 | 用途 |
|------|------|
| `pandoc` | EPUB/MOBI → Markdown（**保留公式**，临时使用） |
| `pdftotext` | PDF → 纯文本（无公式） |
| `marker-pdf` | PDF → Markdown（开源、含公式） |
| `mathpix` | PDF → Markdown（付费、质量最好） |
| `csplit` | 按一级标题拆分大文件 |
| `Grep` | 在原始 .xhtml/.html 中精确定位 |
| `Read` + offset/limit | 读取源文件局部段落 |

**重要**：所有这些工具的输出**仅作临时读取**，**不存储到 memory**。

### 脚本

- `_scripts/extract_epub_chapter.py`：EPUB XHTML → 干净文本（保留 LaTeX）

---

## 11. 失败模式与避免

| 失败模式 | 避免方法 |
|---------|---------|
| 抽取过细，浪费时间 | 严格 P0/P1/P2 优先级 |
| 抽取过粗，未来还要回读 | 关键公式 + 算法步骤必保留 |
| 概念错误 | 数学公式旁注 LaTeX 来源 |
| **公式丢失**（如用 [MATH] 替换） | **绝不丢 LaTeX**——抽取时就写好，源文件当真相 |
| 公式与代码不匹配 | 公式与 C++ 片段并列 |
| 章节间重复 | 用 cross-reference，不要重写 |
| 找不到想要的 | 写好"关键术语索引" |
| 不知道何时回读 | 严格按 9.2 决策树 |
| `_raw/` 巨大占空间 | 不存储预处理文件 |

---

## 12. 总结

- **Memory 是快速入口，源文件是真相**
- **结构化目录 + 优先级 + 日志 + 索引**
- **预处理时**用 pandoc 保留 LaTeX**，但**不存储预处理结果到 memory**
- **公式 = LaTeX + 来源 + C++**（写 memory 时直接写，不依赖预处理）
- **P0/P1/P2 优先级**——不要试图"全部抽取"
- **当 memory 不足或冲突时必须回读源文件**（用 `Read` + offset/limit 直接读原始 .xhtml/.html）
