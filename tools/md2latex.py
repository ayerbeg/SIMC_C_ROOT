#!/usr/bin/env python3
"""md2latex.py - convert the SIMC++ Markdown manuals into the bodies of the LaTeX books (task P15.D5).

usage: md2latex.py --book NAME [--book NAME ...] --source-root DIR --outdir DIR [--figdir DIR]

The Markdown files under docs/manual/user, docs/manual/developer and docs/translation stay the single source of the
text.  This script reads the book manifest docs/latex/books.json (parts, chapter files, appendices, PDF-only front
matter) and writes, for every requested book NAME, into --outdir:
    NAME_meta.tex   title, subtitle, version (git describe), date (commit date) as macros
    NAME_body.tex   front matter (the index.md of the manual + the PDF-only files), parts, chapters, appendices
The LaTeX skeletons docs/latex/<NAME>.tex (class docs/latex/simcpp.cls) \\input these two files.

Supported Markdown (GitHub flavour, the subset used by the manuals):
  ATX headings (# .. ####), paragraphs, **strong**, *emphasis*, _emphasis_, `code spans`, backslash escapes,
  fenced code blocks (``` or ~~~, optional language: bash/sh/console -> command box, cpp/c++/cmake/json -> coloured
  listing, none -> plain listing), bullet and numbered lists (nested), pipe tables (alignment row; `\\|` escapes;
  <br> inside cells), block quotes, callouts (`> **Note:**`, `> **Tip:**`, `> **Warning:**`, `> **Legacy SIMC:**`),
  links (same chapter / other chapter of the book -> \\hyperref; other files -> text + footnote with the path;
  http(s) -> \\href), images (`![caption](../figures/x.png)` alone in a paragraph -> numbered figure; a figure whose
  TikZ source docs/latex/figures/x.tex exists is included as the PDF built from it), autolinks <http://...>,
  inline <br>, HTML comments (`<!-- from: path -->` -> "source" caption of the next listing,
  `<!-- abstract: text -->` -> chapter abstract; other comments are invisible on GitHub and are skipped),
  horizontal rules.
  Chapter abstract: the `<!-- abstract: -->` comment if present, else the first paragraph after the chapter title.
  The parameter reference (tables with the header Name | Units | Default | Meaning | Used by) is typeset as reference
  entries with index entries; parameters, command-line options (`--x`) and the classes of include/simc are indexed.
Anything else (setext headings, raw HTML blocks or tags, reference-style links, footnotes, a table row with the wrong
number of cells, a link to a missing file or heading, ...) is an error: every error is printed as file:line: message
and the exit status is 1.  Nothing is dropped silently.
"""
import argparse
import json
import os
import re
import subprocess
import sys

# ----------------------------------------------------------------------------------------------------------------------
# helpers shared with tools/check_md_links.py (GitHub heading slugs)

FENCE_RE = re.compile(r'^(\s*)(`{3,}|~{3,})\s*([\w+#.-]*)\s*$')
HEADING_RE = re.compile(r'^\s{0,3}(#{1,6})\s+(.*?)\s*#*\s*$')
LIST_RE = re.compile(r'^(\s*)([*+-]|\d{1,9}[.)])(\s+|$)(.*)$')
TABLE_SEP_RE = re.compile(r'^\s*\|?\s*:?-{1,}:?\s*(\|\s*:?-{1,}:?\s*)*\|?\s*$')
HR_RE = re.compile(r'^\s{0,3}([-*_])(\s*\1){2,}\s*$')
CALLOUTS = {'Note': 'simcnote', 'Tip': 'simctip', 'Warning': 'simcwarning', 'Legacy SIMC': 'simclegacy'}
PARAM_HEADER = ['Name', 'Units', 'Default', 'Meaning', 'Used by']
SUPERSCRIPTS = {'⁰': '0', '¹': '1', '⁴': '4', '⁵': '5', '⁶': '6', '⁷': '7', '⁸': '8', '⁹': '9', '⁺': '+',
                '⁻': '−', 'ⁿ': 'n'}
SUBSCRIPTS = {'₀': '0', '₁': '1', '₂': '2', '₃': '3'}


def slugify(text):
    text = re.sub(r'<[^>]+>', '', text)
    text = re.sub(r'\[([^\]]*)\]\([^)]*\)', r'\1', text)
    text = text.replace('`', '').strip().lower()
    text = re.sub(r'[^\w\- ]', '', text, flags=re.UNICODE)
    return text.replace(' ', '-')


def label_safe(s):
    return ''.join(c if (c.isascii() and (c.isalnum() or c in '-.')) else 'u%x' % ord(c) for c in s)


class ConvError(Exception):
    pass


# ----------------------------------------------------------------------------------------------------------------------
# LaTeX escaping

TEXT_ESC = {'\\': r'\textbackslash{}', '{': r'\{', '}': r'\}', '$': r'\$', '&': r'\&', '#': r'\#', '%': r'\%',
            '_': r'\_', '~': r'\textasciitilde{}', '^': r'\textasciicircum{}'}


def esc_text(s):
    out = []
    for c in s:
        if c in TEXT_ESC:
            out.append(TEXT_ESC[c])
        elif c in SUPERSCRIPTS:
            out.append(r'\textsuperscript{%s}' % SUPERSCRIPTS[c])
        elif c in SUBSCRIPTS:
            out.append(r'\textsubscript{%s}' % SUBSCRIPTS[c])
        else:
            out.append(c)
    r = ''.join(out)
    r = r.replace('--', '-{}-').replace('--', '-{}-')      # no en/em-dash ligatures from ASCII hyphens
    r = r.replace('``', '`{}`').replace("''", "'{}'")
    return r


def esc_code(s, breaks=True):
    """Inline code: escape and allow line breaks after punctuation (identifiers and paths can be long)."""
    out = []
    for i, c in enumerate(s):
        if c in TEXT_ESC:
            out.append(TEXT_ESC[c])
        elif c == ' ':
            out.append(r'\ ')
        elif c == '-':
            out.append(r'-{}')
        else:
            out.append(c)
        if breaks and c in '/._:,(=|>' and i + 1 < len(s) and s[i + 1] not in ' :':
            out.append(r'\allowbreak{}')
    return ''.join(out)


def esc_index(s):
    """makeindex: quote the special characters ! @ | " with "."""
    return re.sub(r'([!@|"])', r'"\1', s)


def idx_key(s):
    """Sort key of an index entry: no LaTeX specials (they would be written raw into the .idx file)."""
    return esc_index(re.sub(r'[%#$&{}\\~^]', '.', s))


def plain(s):
    """Visible text of an inline Markdown string (for widths, bookmarks, captions of other books)."""
    s = re.sub(r'!?\[([^\]]*)\]\([^)]*\)', r'\1', s)
    s = s.replace('`', '').replace('**', '').replace('\\|', '|')
    s = re.sub(r'<br>', ' ', s)
    return s


# ----------------------------------------------------------------------------------------------------------------------
# blocks

class Block:
    def __init__(self, kind, line, **kw):
        self.kind, self.line = kind, line
        self.__dict__.update(kw)


def is_table_start(lines, i):
    return (lines[i].lstrip().startswith('|') and i + 1 < len(lines) and TABLE_SEP_RE.match(lines[i + 1])
            and '|' in lines[i + 1])


def split_row(line):
    s = line.strip()
    if s.startswith('|'):
        s = s[1:]
    if s.endswith('|') and not s.endswith('\\|'):
        s = s[:-1]
    cells, cur, i = [], '', 0
    while i < len(s):
        if s[i] == '\\' and i + 1 < len(s) and s[i + 1] == '|':
            cur += '\\|'
            i += 2
            continue
        if s[i] == '|':
            cells.append(cur.strip())
            cur = ''
        else:
            cur += s[i]
        i += 1
    cells.append(cur.strip())
    return cells


def parse_blocks(lines, first_line, path, errors):
    """lines: list of str (no newline). first_line: source line number of lines[0]."""
    blocks, i, n = [], 0, len(lines)
    pending_from = None
    while i < n:
        line = lines[i]
        ln = first_line + i
        if not line.strip():
            i += 1
            continue
        m = FENCE_RE.match(line)
        if m:
            indent, fence, lang = len(m.group(1)), m.group(2), m.group(3).lower()
            body, j = [], i + 1
            while j < n:
                mm = re.match(r'^\s*(`{3,}|~{3,})\s*$', lines[j])
                if mm and mm.group(1)[0] == fence[0] and len(mm.group(1)) >= len(fence):
                    break
                l2 = lines[j]
                body.append(l2[indent:] if l2[:indent].strip() == '' else l2.lstrip())
                j += 1
            if j >= n:
                errors.append(f'{path}:{ln}: unterminated fenced code block')
            blocks.append(Block('code', ln, lang=lang, text=body, source=pending_from))
            pending_from = None
            i = j + 1
            continue
        if line.lstrip().startswith('<!--'):
            text, j = line, i
            while '-->' not in text and j + 1 < n:
                j += 1
                text += '\n' + lines[j]
            if '-->' not in text:
                errors.append(f'{path}:{ln}: unterminated HTML comment')
            inner = text.strip()[4:].rsplit('-->', 1)[0].strip()
            rest = text.split('-->', 1)[1].strip() if '-->' in text else ''
            if rest:
                errors.append(f'{path}:{ln}: text after an HTML comment on the same line')
            mf = re.match(r'^from:\s*(\S+)$', inner)
            ma = re.match(r'^abstract:\s*(.*)$', inner, re.S)
            if mf:
                pending_from = mf.group(1)
            elif ma:
                blocks.append(Block('abstract', ln, text=' '.join(ma.group(1).split())))
            i = j + 1
            continue
        mh = HEADING_RE.match(line)
        if mh and not line.startswith('    '):
            blocks.append(Block('heading', ln, level=len(mh.group(1)), text=mh.group(2)))
            i += 1
            continue
        if HR_RE.match(line):
            blocks.append(Block('hr', ln))
            i += 1
            continue
        if is_table_start(lines, i):
            header = split_row(line)
            seps = split_row(lines[i + 1])
            if len(seps) != len(header):
                errors.append(f'{path}:{ln + 1}: table alignment row has {len(seps)} cells, header has {len(header)}')
            aligns = []
            for s in seps:
                s = s.strip()
                aligns.append('r' if s.endswith(':') and not s.startswith(':') else
                              'c' if s.startswith(':') and s.endswith(':') else 'l')
            rows, j = [], i + 2
            while j < n and lines[j].strip().startswith('|'):
                cells = split_row(lines[j])
                if len(cells) != len(header):
                    errors.append(f'{path}:{first_line + j}: table row has {len(cells)} cells, header has '
                                  f'{len(header)} (unescaped "|"?)')
                rows.append((first_line + j, cells))
                j += 1
            blocks.append(Block('table', ln, header=header, aligns=aligns, rows=rows))
            i = j
            continue
        if line.lstrip().startswith('>'):
            body, j = [], i
            while j < n and lines[j].strip() and lines[j].lstrip().startswith('>'):
                body.append(re.sub(r'^\s*>\s?', '', lines[j]))
                j += 1
            blocks.append(Block('quote', ln, blocks=parse_blocks(body, ln, path, errors), raw=body))
            i = j
            continue
        ml = LIST_RE.match(line)
        if ml and (ml.group(4).strip() or True):
            i = parse_list(lines, i, first_line, path, errors, blocks)
            continue
        if re.match(r'^\s*<(?!br>)[A-Za-z/]', line):
            errors.append(f'{path}:{ln}: raw HTML block is not supported: {line.strip()[:60]}')
        # paragraph
        para, j = [], i
        while j < n and lines[j].strip():
            l2 = lines[j]
            if j > i and (FENCE_RE.match(l2) or HEADING_RE.match(l2) or is_table_start(lines, j)
                          or l2.lstrip().startswith('>') or LIST_RE.match(l2) or l2.lstrip().startswith('<!--')):
                break
            if j > i and (re.match(r'^\s*=+\s*$', l2) or re.match(r'^\s*-{2,}\s*$', l2)):
                errors.append(f'{path}:{first_line + j}: setext heading is not supported')
            para.append(l2.strip())
            j += 1
        blocks.append(Block('para', ln, text='\n'.join(para)))
        i = j
    if pending_from:
        errors.append(f'{path}: "<!-- from: {pending_from} -->" is not followed by a code block')
    return blocks


def parse_list(lines, i, first_line, path, errors, blocks):
    n = len(lines)
    m = LIST_RE.match(lines[i])
    base_indent = len(m.group(1))
    ordered = m.group(2)[0].isdigit()
    start = int(m.group(2)[:-1]) if ordered else 1
    items = []
    while i < n:
        m = LIST_RE.match(lines[i])
        if not m or len(m.group(1)) != base_indent or m.group(2)[0].isdigit() != ordered:
            break
        content_indent = len(m.group(1)) + len(m.group(2)) + len(m.group(3) or ' ')
        item_lines = [m.group(4)]
        item_first = first_line + i
        j = i + 1
        while j < n:
            l2 = lines[j]
            if not l2.strip():
                # blank: the item continues only if the next non-blank line is indented into the item
                k = j
                while k < n and not lines[k].strip():
                    k += 1
                if k < n and len(lines[k]) - len(lines[k].lstrip()) >= content_indent:
                    item_lines.extend([''] * (k - j))
                    j = k
                    continue
                break
            ind = len(l2) - len(l2.lstrip())
            if ind >= content_indent:
                item_lines.append(l2[content_indent:])
            elif LIST_RE.match(l2) and ind <= base_indent + 1:
                break
            elif LIST_RE.match(l2):
                item_lines.append(l2.lstrip())     # nested list written with a smaller indent
            elif FENCE_RE.match(l2) or HEADING_RE.match(l2) or l2.lstrip().startswith('|'):
                break
            else:
                item_lines.append(l2.strip())      # lazy continuation of the paragraph
            j += 1
        items.append(parse_blocks(item_lines, item_first, path, errors))
        i = j
        while i < n and not lines[i].strip():
            k = i
            while k < n and not lines[k].strip():
                k += 1
            if k < n and LIST_RE.match(lines[k]) and len(LIST_RE.match(lines[k]).group(1)) == base_indent:
                i = k
            break
    blocks.append(Block('list', first_line + i, ordered=ordered, start=start, items=items))
    return i


# ----------------------------------------------------------------------------------------------------------------------
# the converter

class Book:
    def __init__(self, name, spec, root, outdir, figdir, all_books, errors):
        self.name, self.spec, self.root, self.outdir, self.figdir = name, spec, root, outdir, figdir
        self.all_books = all_books
        self.errors = errors
        self.dir = os.path.join(root, spec['dir'])
        self.files = []             # (abs path, role) role: index | chapter | appendix
        self.fileid = {}
        self.anchors = {}           # abs path -> {slug: label}
        self.titles = {}            # abs path -> {slug: heading text}
        self.params = {}            # parameter name -> label (only the book with the reference)
        self.cur = None
        self.fig_count = 0

    # -- set-up ------------------------------------------------------------------------------------------------------
    def collect(self):
        idx = os.path.join(self.dir, 'index.md')
        order = []
        for line in open(idx, encoding='utf-8'):
            for t in re.findall(r'\]\(([^)#\s]+\.md)(?:#[^)]*)?\)', line):
                p = os.path.normpath(os.path.join(self.dir, t))
                if os.path.dirname(p) == self.dir and p not in order and p != idx:
                    order.append(p)
        chapters = [os.path.join(self.dir, c) for part in self.spec['parts'] for c in part['chapters']]
        appendices = [os.path.join(self.dir, c) for c in self.spec.get('appendices', [])]
        if set(order) != set(chapters + appendices):
            miss = sorted(set(order) - set(chapters + appendices))
            extra = sorted(set(chapters + appendices) - set(order))
            self.errors.append(f'{idx}: the book manifest does not match the chapters linked from index.md '
                               f'(not in the manifest: {miss}; not linked from index.md: {extra})')
        if [c for c in order if c in chapters] != chapters:
            self.errors.append(f'{idx}: the chapter order of the book manifest differs from index.md')
        self.files = [(idx, 'index')] + [(c, 'chapter') for c in chapters] + [(a, 'appendix') for a in appendices]
        for k, (p, _) in enumerate(self.files):
            self.fileid[p] = '%s-%s' % (self.name, label_safe(os.path.splitext(os.path.basename(p))[0]))
            self.anchors[p], self.titles[p] = heading_anchors(p)

    def label_for(self, path, slug=None):
        if slug is None:
            return self.fileid[path]
        return self.fileid[path] + ':' + label_safe(slug)

    # -- inline ------------------------------------------------------------------------------------------------------
    def inline(self, s, ln, ctx='text'):
        """ctx: text | table | heading | caption."""
        out, buf, i, n = [], [], 0, len(s)

        def flush():
            if buf:
                out.append(esc_text(''.join(buf)))
                buf.clear()

        while i < n:
            c = s[i]
            if c == '\\' and i + 1 < n and s[i + 1] in '\\`*_{}[]()#+-.!|<>~"\'':
                buf.append(s[i + 1])
                i += 2
                continue
            if c == '`':
                k = i
                while k < n and s[k] == '`':
                    k += 1
                ticks = s[i:k]
                j = s.find(ticks, k)
                while j != -1 and j + len(ticks) < n and s[j + len(ticks)] == '`':
                    j = s.find(ticks, j + len(ticks) + 1)
                if j == -1:
                    buf.append(ticks)
                    i = k
                    continue
                code = s[k:j].replace('\n', ' ')
                if ctx == 'table':
                    code = code.replace('\\|', '|')
                if len(code) > 1 and code[0] == ' ' and code[-1] == ' ' and code.strip():
                    code = code[1:-1]
                flush()
                out.append(self.code_span(code, ctx))
                i = j + len(ticks)
                continue
            if c == '!' and i + 1 < n and s[i + 1] == '[':
                lk = parse_link(s, i + 1)
                if lk:
                    text, target, end = lk
                    flush()
                    out.append(self.image(text, target, ln, inline=True))
                    i = end
                    continue
            if c == '[':
                lk = parse_link(s, i)
                if lk:
                    text, target, end = lk
                    flush()
                    out.append(self.link(text, target, ln, ctx))
                    i = end
                    continue
                if re.match(r'\[[^\]]+\]\[[^\]]*\]', s[i:]) or re.match(r'\[\^', s[i:]):
                    self.errors.append(f'{self.cur}:{ln}: reference-style links and footnotes are not supported')
            if c == '<':
                m = re.match(r'<(https?://[^>\s]+)>', s[i:])
                if m:
                    flush()
                    out.append(r'\url{%s}' % m.group(1).replace('%', r'\%').replace('#', r'\#'))
                    i += m.end()
                    continue
                m = re.match(r'<br\s*/?>', s[i:])
                if m:
                    flush()
                    out.append(r'\newline{}' if ctx == 'table' else r'\\' if ctx == 'text' else ' ')
                    i += m.end()
                    continue
                m = re.match(r'</?[A-Za-z][A-Za-z0-9-]*(\s[^>]*)?/?>', s[i:])
                if m:
                    self.errors.append(f'{self.cur}:{ln}: inline HTML tag is not supported: {m.group(0)}')
            if c in '*_':
                k = i
                while k < n and s[k] == c:
                    k += 1
                run = k - i
                prev = s[i - 1] if i > 0 else ' '
                nxt = s[k] if k < n else ' '
                can_open = not nxt.isspace() and (c == '*' or not prev.isalnum())
                if can_open and run in (1, 2, 3):
                    close = find_closing(s, k, c, run)
                    if close != -1:
                        flush()
                        inner = self.inline(s[k:close], ln, ctx)
                        if run == 1:
                            out.append(r'\emph{%s}' % inner)
                        elif run == 2:
                            out.append(r'\textbf{%s}' % inner)
                        else:
                            out.append(r'\textbf{\emph{%s}}' % inner)
                        i = close + run
                        continue
                buf.append(s[i:k])
                i = k
                continue
            if c == '\n':
                buf.append(' ' if ctx != 'text' else '\n')
                i += 1
                continue
            buf.append(c)
            i += 1
        flush()
        return ''.join(out)

    def code_span(self, code, ctx):
        """Inline code; parameters, command-line options and classes get index entries (and parameters a link)."""
        tex = r'\code{%s}' % esc_code(code)
        if ctx in ('heading', 'caption', 'bookmark'):
            return tex
        idx = ''
        key = code.strip()
        if key in self.params:
            idx = r'\index{Parameters!%s@\code{%s}}' % (idx_key(key), esc_index(esc_code(key, False)))
            if self.params[key]:
                tex = r'\hyperref[%s]{%s}' % (self.params[key], tex)
        else:
            mo = re.match(r'^(--[a-z][a-z0-9-]*)(?:[ =].*)?$', key)
            if mo and mo.group(1) in OPTIONS:
                o = mo.group(1)
                idx = r'\index{Command-line options!%s@\code{%s}}' % (esc_index(o[2:]), esc_index(esc_code(o, False)))
            else:
                mc = re.match(r'^(?:simc::)?(?:[a-z]+::)?([A-Z][A-Za-z0-9]+)(?:::\w+(?:\(\))?|\(\))?$', key)
                if mc and mc.group(1) in CLASSES:
                    k = mc.group(1)
                    idx = r'\index{Classes and interfaces!%s@\code{%s}}' % (k, k)
                else:
                    mr = re.match(r'^SIMC_REGISTER$', key)
                    if mr:
                        idx = r'\index{SIMC_REGISTER@\code{SIMC\_REGISTER}}'
        return tex + idx

    def link(self, text, target, ln, ctx):
        inner = self.inline(text, ln, 'caption' if ctx == 'heading' else ctx)
        if ctx == 'heading':
            return inner
        if re.match(r'^(https?|ftp)://|^mailto:', target):
            if plain(text).strip() == target:
                return r'\url{%s}' % target.replace('%', r'\%').replace('#', r'\#')
            return r'\href{%s}{%s}' % (target.replace('%', r'\%').replace('#', r'\#'), inner)
        file_part, _, anchor = target.partition('#')
        base = os.path.dirname(self.cur)
        path = os.path.normpath(os.path.join(base, file_part)) if file_part else self.cur
        if not os.path.exists(path):
            self.errors.append(f'{self.cur}:{ln}: link target not found: {target}')
            return inner
        if path in self.fileid:
            if anchor:
                if anchor not in self.anchors[path]:
                    self.errors.append(f'{self.cur}:{ln}: no heading #{anchor} in {target}')
                    return inner
                return r'\hyperref[%s]{%s}' % (self.label_for(path, anchor), inner)
            return r'\hyperref[%s]{%s}' % (self.label_for(path), inner)
        # outside this book
        if anchor and path.endswith('.md') and anchor not in heading_anchors(path)[0]:
            self.errors.append(f'{self.cur}:{ln}: no heading #{anchor} in {target}')
        rel = os.path.relpath(path, self.root)
        if os.path.isdir(path):
            rel += '/'
        where = None
        for bname, b in self.all_books.items():
            if b is self:
                continue
            bdir = os.path.join(self.root, b['dir'])
            if os.path.dirname(path) == bdir and path.endswith('.md'):
                where = b['title']
                break
        note = r'\code{%s}' % esc_code(rel)
        if where:
            heading = heading_anchors(path)[1].get(anchor) if anchor else first_heading(path)
            note = r'%s, ``%s'' (%s)' % (esc_text(where), self.inline(heading or '', ln, 'caption'), note)
        if ctx in ('caption',):
            return inner
        return r'%s\footnote{%s}' % (inner, note)

    def image(self, alt, target, ln, inline=False):
        path = os.path.normpath(os.path.join(os.path.dirname(self.cur), target))
        if not os.path.exists(path):
            self.errors.append(f'{self.cur}:{ln}: image not found: {target}')
            return ''
        stem = os.path.splitext(os.path.basename(path))[0]
        tikz = os.path.join(self.root, 'docs', 'latex', 'figures', stem + '.tex')
        if os.path.exists(tikz):
            gfx = os.path.join(self.figdir, stem + '.pdf').replace('\\', '/')
        else:
            gfx = path
        if inline:
            return r'\includegraphics[height=1.2em]{%s}' % gfx
        cap = self.inline(alt, ln, 'caption')
        short = esc_text(plain(alt).split('. ')[0].split(': ')[0])
        return ('\\begin{figure}[htbp]\n\\simcfigure{%s}\n\\caption[%s]{%s}\\label{fig:%s}\n\\end{figure}\n'
                % (gfx, short, cap, label_safe(stem)))

    # -- blocks ------------------------------------------------------------------------------------------------------
    def render_blocks(self, blocks, out, depth=0):
        k = 0
        while k < len(blocks):
            b = blocks[k]
            getattr(self, 'r_' + b.kind)(b, out)
            k += 1

    def r_abstract(self, b, out):
        out.append(r'\simcabstract{%s}' % self.inline(b.text, b.line) + '\n')

    def r_heading(self, b, out):
        self.heading(b, out)

    def r_hr(self, b, out):
        out.append('\\simcrule\n')

    def r_para(self, b, out):
        t = b.text.strip()
        m = re.match(r'^!\[((?:[^\[\]]|\[[^\]]*\])*)\]\(([^)\s]+)\)$', t)
        if m:
            out.append(self.image(m.group(1), m.group(2), b.line))
            return
        out.append(self.inline(t, b.line) + '\n\n')

    def r_code(self, b, out):
        lang = b.lang
        text = '\n'.join(b.text).rstrip('\n')
        if '\\end{lstlisting}' in text:
            self.errors.append(f'{self.cur}:{b.line}: code block contains \\end{{lstlisting}}')
        width = max([len(l) for l in b.text] or [0])
        size = 'normal' if width <= 92 else 'small' if width <= 106 else 'tiny'
        env = {'bash': 'simcshell', 'sh': 'simcshell', 'console': 'simcshell', 'shell': 'simcshell',
               'cpp': 'simccode', 'c++': 'simccode', 'cmake': 'simccode', 'json': 'simccode'}.get(lang)
        if lang and env is None:
            self.errors.append(f'{self.cur}:{b.line}: unknown code block language "{lang}"')
            env = 'simcplain'
        if env is None:
            env = 'simcplain'
        style = {'cpp': 'simccpp', 'c++': 'simccpp', 'cmake': 'simccmake', 'json': 'simcjson'}.get(lang, '')
        opts = ['basicstyle=\\ttfamily\\%s' % {'normal': 'footnotesize', 'small': 'scriptsize', 'tiny': 'tiny'}[size]]
        if style:
            opts.append('style=%s' % style)
        if b.source:
            src = os.path.join(self.root, b.source)
            if not os.path.exists(src):
                self.errors.append(f'{self.cur}:{b.line}: "<!-- from: {b.source} -->": file not found')
            out.append(r'\simcsource{%s}' % esc_code(b.source) + '\n')
        out.append('\\begin{%s}[%s]\n%s\n\\end{%s}\n' % (env, ','.join(opts), text, env))

    def r_quote(self, b, out):
        inner = b.blocks
        if inner and inner[0].kind == 'para':
            m = re.match(r'^\*\*(Note|Tip|Warning|Legacy SIMC):\*\*\s*(.*)$', inner[0].text, re.S)
            if m:
                env = CALLOUTS[m.group(1)]
                first = Block('para', inner[0].line, text=m.group(2))
                sub = []
                self.render_blocks(([first] if m.group(2).strip() else []) + inner[1:], sub)
                out.append('\\begin{%s}\n%s\\end{%s}\n' % (env, ''.join(sub).rstrip() + '\n', env))
                return
        sub = []
        self.render_blocks(inner, sub)
        out.append('\\begin{quote}\n%s\\end{quote}\n' % ''.join(sub))

    def r_list(self, b, out):
        env = 'enumerate' if b.ordered else 'itemize'
        opt = '[start=%d]' % b.start if b.ordered and b.start != 1 else ''
        out.append('\\begin{%s}%s\n' % (env, opt))
        for item in b.items:
            sub = []
            self.render_blocks(item, sub)
            out.append('\\item ' + ''.join(sub).strip() + '\n')
        out.append('\\end{%s}\n' % env)

    def r_table(self, b, out):
        if [plain(h) for h in b.header] == PARAM_HEADER:
            return self.param_table(b, out)
        ncol = len(b.header)
        rows = [cells + [''] * (ncol - len(cells)) for _, cells in b.rows]
        lengths = [[len(plain(c)) for c in [b.header[j]] + [r[j] for r in rows]] for j in range(ncol)]
        natural = [max(max(len(p) for p in re.split(r'<br>', c) or ['']) if c else 0
                       for c in [plain(b.header[j])] + [r[j] for r in rows]) for j in range(ncol)]
        natural = [max(4, x + 1) for x in natural]
        words = [min(16, max([len(w) for c in [b.header[j]] + [r[j] for r in rows]
                              for w in re.split(r'[\s/_.,:]+', plain(c))] or [3])) + 1 for j in range(ncol)]
        total = [sum(l) for l in lengths]
        avg = sum(sum(l) / len(l) for l in lengths)
        landscape = ncol >= 7 or avg > 260
        size = 'footnotesize' if (landscape or ncol >= 5 or avg > 120) else 'small'
        avail = {('footnotesize', True): 205, ('footnotesize', False): 118, ('small', False): 102}[(size, landscape)]
        avail -= 2 * ncol      # column separation
        widths = allocate(natural, total, words, avail)
        scale = 1.0 if landscape else 1.0
        fracs = [w / (avail + 2 * ncol) * scale for w in widths]
        colspec = ''.join(r'>{\%s\arraybackslash}p{\dimexpr %.3f\linewidth-2\tabcolsep\relax}'
                          % ({'r': 'raggedleft', 'c': 'centering'}.get(b.aligns[j] if j < len(b.aligns) else 'l',
                                                                     'raggedright'), fracs[j])
                          for j in range(ncol))
        head = ' & '.join(r'\simcth{%s}' % self.inline(h, b.line, 'table') for h in b.header) + r' \\'
        caption = ''
        if len(rows) >= 8 and self.last_heading:
            caption = r'\caption{%s}' % self.last_heading + r'\\'
        body = []
        for (ln, _), r in zip(b.rows, rows):
            body.append(' & '.join(self.inline(c, ln, 'table') for c in r) + r' \\')
        t = []
        t.append(r'\begin{simctable}{%s}{%s}{%s}' % (size, ncol, colspec))
        t.append(caption + '%')
        t.append(r'\toprule\rowcolor{simchead}' + head + r'\midrule\endfirsthead')
        t.append(r'\multicolumn{%d}{@{}l}{\simccontinued}\\\toprule\rowcolor{simchead}' % ncol + head
                 + r'\midrule\endhead')
        t.append(r'\midrule\multicolumn{%d}{r@{}}{\simccontinues}\\\endfoot' % ncol)
        t.append(r'\bottomrule\endlastfoot')
        t.extend(body)
        t.append(r'\end{simctable}')
        s = '\n'.join(t) + '\n'
        if landscape:
            s = '\\begin{landscape}\n' + s + '\\end{landscape}\n'
        out.append(s)

    def param_table(self, b, out):
        out.append('\\begin{simcparams}\n')
        for ln, cells in b.rows:
            name = plain(cells[0]).strip()
            label = self.params.get(name) or 'param:' + label_safe(name)
            out.append(r'\simcparam{%s}{%s}{%s}{%s}{%s}{%s}{%s}' % (
                esc_code(name), label, self.inline(cells[1], ln, 'table'), self.inline(cells[2], ln, 'table'),
                self.inline(cells[4], ln, 'table'), self.inline(cells[3], ln, 'table'),
                idx_key(name) + '@' + r'\code{' + esc_index(esc_code(name, False)) + '}') + '\n')
        out.append('\\end{simcparams}\n')

    # -- headings ----------------------------------------------------------------------------------------------------
    def heading(self, b, out):
        level, text = b.level, b.text
        slug = self.slugs.pop(0) if self.slugs else slugify(text)
        label = self.label_for(self.cur, slug)
        self.last_heading = self.inline(strip_number(text)[1], b.line, 'heading')
        if level == 1:
            self.errors.append(f'{self.cur}:{b.line}: more than one level-1 heading')
            return
        num, rest = strip_number(text)
        cmd = {2: 'section', 3: 'subsection', 4: 'subsubsection'}.get(level)
        if cmd is None:
            self.errors.append(f'{self.cur}:{b.line}: heading level {level} is not supported')
            return
        self.counters[level] = self.counters.get(level, 0) + 1
        for deeper in range(level + 1, 7):
            self.counters[deeper] = 0
        expected = '.'.join(str(self.counters.get(lv, 0)) for lv in range(2, level + 1))
        expected = (str(self.chapter_no) + '.' + expected) if self.chapter_no else None
        numbered = (num is not None and expected is not None and num.rstrip('.') == expected and self.numbered)
        if not numbered:
            self.counters[level] -= 1
            rest = text
        tex = self.inline(rest, b.line, 'heading')
        bm = esc_text(plain(rest))
        if numbered:
            out.append('\\%s{\\texorpdfstring{%s}{%s}}\\label{%s}\n' % (cmd, tex, bm, label))
        else:
            toc = '\\addcontentsline{toc}{%s}{\\texorpdfstring{%s}{%s}}' % (cmd, tex, bm) if level == 2 else ''
            out.append('\\%s*{\\texorpdfstring{%s}{%s}}%s\\label{%s}\n' % (cmd, tex, bm, toc, label))

    # -- files -------------------------------------------------------------------------------------------------------
    def render_file(self, path, role, chapter_no):
        self.cur = path
        self.counters = {}
        self.chapter_no = chapter_no
        self.numbered = role == 'chapter'
        self.last_heading = ''
        with open(path, encoding='utf-8') as f:
            lines = f.read().split('\n')
        blocks = parse_blocks(lines, 1, path, self.errors)
        self.slugs = heading_slug_list(path)
        out = []
        if not blocks or blocks[0].kind != 'heading' or blocks[0].level != 1:
            self.errors.append(f'{path}:1: a chapter file must start with a level-1 heading')
            return ''
        h1 = blocks.pop(0)
        self.slugs.pop(0)
        num, title = strip_number(h1.text)
        if role == 'chapter' and (num is None or num.rstrip('.') != str(chapter_no)):
            self.errors.append(f'{path}:{h1.line}: chapter heading number "{num}" differs from its position '
                               f'{chapter_no} in the book')
        if role == 'index':
            title = self.spec.get('index_title', 'About this book')
        tex = self.inline(title, h1.line, 'heading')
        bm = esc_text(plain(title))
        self.last_heading = tex
        if role == 'index':
            out.append('\\chapter*{%s}\\addcontentsline{toc}{chapter}{%s}\\markboth{%s}{%s}\\label{%s}\n'
                       % (tex, bm, bm, bm, self.label_for(path)))
        else:
            out.append('\\chapter{\\texorpdfstring{%s}{%s}}\\label{%s}\n' % (tex, bm, self.label_for(path)))
        # abstract: explicit comment, else the first paragraph (typeset as the abstract, not repeated)
        abstract = [k for k, bl in enumerate(blocks[:3]) if bl.kind == 'abstract']
        if not abstract and role != 'index' and blocks and blocks[0].kind == 'para' and \
                not blocks[0].text.startswith('!['):
            blocks[0] = Block('abstract', blocks[0].line, text=blocks[0].text)
        self.render_blocks(blocks, out)
        return ''.join(out)

    def convert(self, front_tex):
        body = ['\\simcfrontmatter\n']
        idx = self.files[0][0]
        body.append(self.render_file(idx, 'index', None))
        for f in self.spec.get('front', []):
            if not os.path.exists(os.path.join(self.root, 'docs', 'latex', f)):
                self.errors.append(f'docs/latex/books.json: front-matter file not found: docs/latex/{f}')
            body.append('\\input{%s}\n' % f)
        if front_tex:
            body.append(front_tex)
        body.append('\\simcmainmatter\n')
        chapter_no = 0
        for part in self.spec['parts']:
            body.append('\\simcpart{%s}{%s}\n' % (esc_text(part['title']), self.inline(part.get('intro', ''), 0)))
            for c in part['chapters']:
                chapter_no += 1
                body.append(self.render_file(os.path.join(self.dir, c), 'chapter', chapter_no))
        if self.spec.get('appendices'):
            body.append('\\simcappendix\n')
            for a in self.spec['appendices']:
                body.append(self.render_file(os.path.join(self.dir, a), 'appendix', None))
        body.append('\\simcbackmatter\n')
        return ''.join(body)


def parse_link(s, i):
    """s[i] == '['. Returns (text, target, end) or None."""
    depth, j = 0, i
    while j < len(s):
        if s[j] == '\\':
            j += 2
            continue
        if s[j] == '`':
            k = s.find('`', j + 1)
            if k == -1:
                return None
            j = k + 1
            continue
        if s[j] == '[':
            depth += 1
        elif s[j] == ']':
            depth -= 1
            if depth == 0:
                break
        j += 1
    if j >= len(s) or j + 1 >= len(s) or s[j + 1] != '(':
        return None
    m = re.match(r'\(\s*<?([^)\s>]+)>?(?:\s+"[^"]*")?\s*\)', s[j + 1:])
    if not m:
        return None
    return s[i + 1:j], m.group(1), j + 1 + m.end()


def find_closing(s, k, c, run):
    j, n = k, len(s)
    while j < n:
        if s[j] == '\\':
            j += 2
            continue
        if s[j] == '`':
            m = re.match(r'`+', s[j:])
            e = s.find(m.group(0), j + len(m.group(0)))
            j = (e + len(m.group(0))) if e != -1 else j + len(m.group(0))
            continue
        if s[j] == '[':
            lk = parse_link(s, j)
            if lk:
                j = lk[2]
                continue
        if s[j] == c:
            e = j
            while e < n and s[e] == c:
                e += 1
            if e - j == run or (e - j > run and run == 1 and False):
                prev = s[j - 1]
                nxt = s[e] if e < n else ' '
                if not prev.isspace() and (c == '*' or not nxt.isalnum()):
                    return j
            j = e
            continue
        j += 1
    return -1


def strip_number(text):
    m = re.match(r'^(\d+(?:\.\d+)*\.?)\s+(.*)$', text)
    if m:
        return m.group(1), m.group(2)
    return None, text


_heading_cache = {}


def heading_slug_list(path):
    slugs, seen, in_fence = [], {}, False
    for line in open(path, encoding='utf-8'):
        if FENCE_RE.match(line.rstrip('\n')):
            in_fence = not in_fence
            continue
        if in_fence:
            continue
        m = HEADING_RE.match(line)
        if m and not line.startswith('    '):
            s = slugify(m.group(2))
            cnt = seen.get(s, 0)
            seen[s] = cnt + 1
            slugs.append(s if cnt == 0 else f'{s}-{cnt}')
    return slugs


def heading_anchors(path):
    if path in _heading_cache:
        return _heading_cache[path]
    slugs = heading_slug_list(path)
    texts, in_fence, k = {}, False, 0
    for line in open(path, encoding='utf-8'):
        if FENCE_RE.match(line.rstrip('\n')):
            in_fence = not in_fence
            continue
        if in_fence:
            continue
        m = HEADING_RE.match(line)
        if m and not line.startswith('    '):
            texts[slugs[k]] = m.group(2)
            k += 1
    _heading_cache[path] = (set(slugs), texts)
    return _heading_cache[path]


def first_heading(path):
    for line in open(path, encoding='utf-8'):
        m = HEADING_RE.match(line)
        if m:
            return m.group(2)
    return ''


def allocate(natural, total, words, avail):
    """Column widths (characters): natural widths if they fit, else water-filling weighted by the content."""
    n = len(natural)
    if sum(natural) <= avail:
        return list(natural)
    widths = [None] * n
    remaining = avail
    free = list(range(n))
    while True:
        weight = sum(total[j] for j in free) or 1
        changed = False
        for j in list(free):
            share = remaining * total[j] / weight
            if natural[j] <= share:
                widths[j] = natural[j]
                remaining -= natural[j]
                free.remove(j)
                changed = True
        if not changed:
            break
    weight = sum(total[j] for j in free) or 1
    for j in free:
        widths[j] = max(words[j], min(natural[j], remaining * total[j] / weight), 6)
    # keep the sum within the available width
    s = sum(widths)
    if s > avail:
        widths = [w * avail / s for w in widths]
    return widths


# ----------------------------------------------------------------------------------------------------------------------
# global vocabularies for the index

CLASSES, OPTIONS = set(), set()


def load_vocabulary(root):
    for d, _, names in os.walk(os.path.join(root, 'include', 'simc')):
        for nm in names:
            if nm.endswith('.h'):
                for line in open(os.path.join(d, nm), encoding='utf-8'):
                    m = re.match(r'^\s*(?:class|struct)\s+([A-Z][A-Za-z0-9]+)\b(?!;)', line)
                    if m and not line.rstrip().endswith(';'):
                        CLASSES.add(m.group(1))
    drv = os.path.join(root, 'docs', 'DRIVER.md')
    usage = os.path.join(root, 'docs', 'manual', 'user', 'VERIFIED_EXAMPLES.md')
    for p in (drv, usage):
        if os.path.exists(p):
            OPTIONS.update(re.findall(r'(--[a-z][a-z0-9-]+)', open(p, encoding='utf-8').read()))


def parameter_names(root):
    names, in_ref = [], False
    for line in open(os.path.join(root, 'docs', 'manual', 'user', '04_inputs.md'), encoding='utf-8'):
        if line.startswith('|') and [plain(c) for c in split_row(line)] == PARAM_HEADER:
            in_ref = True
            continue
        if not line.startswith('|'):
            in_ref = False
            continue
        m = re.match(r'^\|\s*`([^`]+)`\s*\|', line)
        if m and in_ref:
            names.append(m.group(1))
    return names


def git_meta(root):
    def run(*args):
        try:
            return subprocess.run(['git', '-C', root] + list(args), capture_output=True, text=True,
                                  check=True).stdout.strip()
        except (OSError, subprocess.CalledProcessError):
            return ''
    return (run('describe', '--always', '--dirty', '--tags') or 'unknown'), (run('log', '-1', '--format=%cs') or '')


# ----------------------------------------------------------------------------------------------------------------------
# translation report: the dashboard page, from the Markdown data

def dashboard(root, figdir, errors):
    cat = os.path.join(root, 'docs', 'translation', '06_catalogue.md')
    rmap = os.path.join(root, 'docs', 'translation', 'fortran_routine_map.md')
    val = os.path.join(root, 'docs', 'translation', '02_validation.md')
    classes = []
    for line in open(cat, encoding='utf-8'):
        m = re.match(r'^\|\s*(a|a-np|b|c|d|dev)\s*\|\s*(.*?)\s*\|\s*(\d+)\s*\|\s*$', line)
        if m:
            classes.append((m.group(1), m.group(2), int(m.group(3))))
    m = re.search(r'Totals:\s*(\d+) routines,\s*(\d+) cited by C\+\+ code,\s*(\d+) without',
                  open(rmap, encoding='utf-8').read())
    levels = []
    in_tab = False
    for line in open(val, encoding='utf-8'):
        if line.startswith('| level |'):
            in_tab = True
            continue
        if in_tab:
            if not line.startswith('|'):
                break
            cells = split_row(line)
            if cells and not set(cells[0]) <= set('-: '):
                levels.append(cells)
    if not classes or not m or not levels:
        errors.append(f'{cat}: cannot read the data of the dashboard page (class table, routine totals, '
                      f'validation levels)')
        return ''
    total = sum(c[2] for c in classes)
    mx = max(c[2] for c in classes) or 1
    bars = []
    for k, (c, meaning, cnt) in enumerate(classes):
        short = meaning.split(':')[0]
        bars.append(r'\simcdashbar{%d}{%s}{%s}{%d}{%.4f}' % (k, esc_text(c), esc_text(short), cnt, cnt / mx))
    lv = [r'\simcdashlevel{%s}{%s}' % (esc_text(plain(cells[0])), esc_text(plain(cells[1]))) for cells in levels]
    ladder = os.path.join(figdir, 'validation_ladder.pdf').replace('\\', '/')
    return ('\\chapter*{The translation at a glance}\\addcontentsline{toc}{chapter}{The translation at a glance}'
            '\\markboth{The translation at a glance}{}\\label{dashboard}\n'
            '\\simcabstract{The numbers of this page are read from the report itself when the book is built: the class '
            'table of the catalogue (chapter~\\ref{translation_report-06_catalogue}), the totals of the generated routine '
            'map (appendix~\\ref{translation_report-fortran_routine_map}) and the levels of the validation framework '
            '(chapter~\\ref{translation_report-02_validation}).}\n'
            '\\simcdashtiles{%d}{%s}{%s}{%s}{%d}\n'
            '\\section*{Findings by class}\n\\begin{simcdashbars}\n%s\n\\end{simcdashbars}\n'
            '\\section*{The validation ladder}\n\\begin{center}\\includegraphics[width=0.92\\linewidth]{%s}\\end{center}\n'
            '\\begin{simcdashlevels}\n%s\n\\end{simcdashlevels}\n\\clearpage\n'
            % (total, m.group(1), m.group(2), m.group(3), len(levels), '\n'.join(bars), ladder, '\n'.join(lv)))


# ----------------------------------------------------------------------------------------------------------------------

def main(argv):
    ap = argparse.ArgumentParser(description='Markdown manuals -> LaTeX book bodies (SIMC++)')
    ap.add_argument('--book', action='append', required=True)
    ap.add_argument('--source-root', required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--figdir', default='figures', help='directory of the figure PDFs as seen from LaTeX')
    a = ap.parse_args(argv[1:])
    root = os.path.abspath(a.source_root)
    with open(os.path.join(root, 'docs', 'latex', 'books.json'), encoding='utf-8') as f:
        books = json.load(f)
    load_vocabulary(root)
    params = parameter_names(root)
    version, date = git_meta(root)
    os.makedirs(a.outdir, exist_ok=True)
    errors = []
    for name in a.book:
        if name not in books:
            print(f'md2latex: unknown book {name} (known: {", ".join(books)})', file=sys.stderr)
            return 2
        spec = books[name]
        bk = Book(name, spec, root, a.outdir, a.figdir, books, errors)
        bk.collect()
        has_ref = any('04_inputs.md' in c for p in spec['parts'] for c in p['chapters'])
        bk.params = {p: ('param:' + label_safe(p)) if has_ref else None for p in params}
        front = dashboard(root, a.figdir, errors) if spec.get('dashboard') else ''
        body = bk.convert(front)
        with open(os.path.join(a.outdir, name + '_body.tex'), 'w', encoding='utf-8') as f:
            f.write('% generated by tools/md2latex.py from the Markdown sources: do not edit\n' + body)
        with open(os.path.join(a.outdir, name + '_meta.tex'), 'w', encoding='utf-8') as f:
            f.write('%% generated by tools/md2latex.py\n\\simcsetmeta{%s}{%s}{%s}{%s}{%s}\n' % (
                esc_text(spec['title']), esc_text(spec['subtitle']), esc_text(version), esc_text(date),
                esc_text(spec.get('kicker', 'SIMC++'))))
    for e in errors:
        print(e)
    print(f'md2latex: {len(a.book)} book(s), {len(errors)} error(s)')
    return 1 if errors else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
