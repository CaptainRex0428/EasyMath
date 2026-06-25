"""
extract_epub_chapter.py
Extract clean text from an EPUB XHTML/HTML chapter file.

IMPORTANT: Preserves LaTeX formulas from <script type="math/tex"> tags.
This is critical — do NOT replace formulas with [MATH] placeholders.

Usage:
    python extract_epub_chapter.py <input.html> <output.txt>

The output is a clean Markdown-like text with headings (#) and paragraphs.
LaTeX formulas (from <script type="math/tex">) are preserved as $...$ inline or
$$...$$ block form so they can be re-used in summary md files.
"""
import re
import html
import sys

def clean(content: str) -> str:
    # Remove MathJax rendered SVG (visual noise)
    # But preserve <script type="math/tex"> tags with original LaTeX
    content = re.sub(r'<mjx-container[^>]*>.*?</mjx-container>', '', content, flags=re.DOTALL)
    content = re.sub(r'<svg[^>]*>.*?</svg>', '', content, flags=re.DOTALL)
    content = re.sub(r'<style[^>]*>.*?</style>', '', content, flags=re.DOTALL)
    # Remove scripts that are NOT math/tex
    content = re.sub(
        r'<script(?![^>]*math/tex)[^>]*>.*?</script>',
        '', content, flags=re.DOTALL
    )
    return content

def preserve_math_tex(content: str) -> str:
    """Convert <script type="math/tex">LaTeX</script> to $$...$$ so it's preserved."""
    # Block math: <script type="math/tex; mode=display">Latex</script>
    def replace_block(m):
        latex = m.group(2)
        return f'$${latex}$$'
    content = re.sub(
        r'<script[^>]*type=["\']math/tex;\s*mode=display["\'][^>]*>(.*?)</script>',
        replace_block, content, flags=re.DOTALL
    )
    # Inline math: <script type="math/tex">Latex</script>
    def replace_inline(m):
        latex = m.group(1)
        return f'${latex}$'
    content = re.sub(
        r'<script[^>]*type=["\']math/tex["\'][^>]*>(.*?)</script>',
        replace_inline, content, flags=re.DOTALL
    )
    return content

def extract(content: str):
    out = []
    for m in re.finditer(
        r'<h([1-6])[^>]*>(.*?)</h\1>|<p[^>]*>(.*?)</p>|<li[^>]*>(.*?)</li>',
        content, re.DOTALL
    ):
        if m.group(1):
            level = int(m.group(1))
            text = re.sub(r'<[^>]+>', '', m.group(2))
            out.append('#' * (level + 1) + ' ' + html.unescape(text).strip())
        elif m.group(3):
            text = re.sub(r'<[^>]+>', '', m.group(3))
            s = html.unescape(text).strip()
            if s:
                out.append(s)
        elif m.group(4):
            text = re.sub(r'<[^>]+>', '', m.group(4))
            s = html.unescape(text).strip()
            if s:
                out.append('- ' + s)
    return out

def main():
    if len(sys.argv) != 3:
        print("Usage: python extract_epub_chapter.py <input.html> <output.txt>")
        sys.exit(1)
    in_path, out_path = sys.argv[1], sys.argv[2]
    with open(in_path, 'r', encoding='utf-8') as f:
        content = f.read()
    content = clean(content)
    content = preserve_math_tex(content)
    blocks = extract(content)
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write('\n\n'.join(blocks))
    print(f"Wrote {len(blocks)} blocks, {sum(len(b) for b in blocks)} chars to {out_path}")

if __name__ == '__main__':
    main()
