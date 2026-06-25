# Extraction Scripts

> 辅助 reference 提取的可复用脚本。所有脚本默认 Windows + Python 3.14+。

---

## 1. `extract_epub_chapter.py` — 提取单章节（**保留 LaTeX**）

**重要**：保留所有数学公式！不替换为 `[MATH]`。

### 用法

```bash
py extract_epub_chapter.py <input.xhtml/html> <output.txt>
```

### 实现

读取 EPUB 章节文件（.xhtml 或 .html），移除 MathJax 渲染块（SVG 噪音），但**保留** `<script type="math/tex">` 中的 LaTeX 源。提取所有 `<h1-h6>` 标题和 `<p>` 段落。

### 局限

- 公式被渲染为 `<mjx-container>` 后无法恢复原始 LaTeX
- 这种情况应优先用 **pandoc 转换**（见 SOP §4.1）
- 此脚本仅作为 fallback，用于 pandoc 不可用时

### 历史

- v1（错误）：替换为 `[MATH]`，丢公式——**已废弃**
- v2（当前）：保留 `<script type="math/tex">` 中的 LaTeX

---

## 2. `batch_extract.bat` — 批量提取所有章节

```bash
cmd //c batch_extract.bat
```

遍历 `ch001.xhtml` ~ `ch026.xhtml`，逐个调用 `extract_epub_chapter.py`，输出到 `_raw/chNNN.txt`。

### 修改以适配不同 reference

修改脚本中的 `SRC` 和 `RAW` 变量。

---

## 3. **推荐**：`pandoc` 转换

```bash
# 整书转换
pandoc -f epub -t markdown --wrap=none \
       --extract-media=./_raw/media \
       -o _raw/full.md "path/to/book.epub"

# 拆分章节
csplit -f _raw/ch_ -b "%02d.md" _raw/full.md '/^# /' '{*}' 2>/dev/null
```

**优势**：
- 保留 LaTeX 公式（`$...$` / `$$...$$`）
- 保留章节结构
- 移除 HTML 噪音
- 体积减少 10-100 倍

**这是 SOP 推荐的标准做法。**

---

## 未来添加

- `extract_pdf.py`（PDF 提取，公式保留）—— 待 marker/mathpix 选定后写
- `validate_chapter.py`（检查章节 md 是否包含所有要求的段落）
