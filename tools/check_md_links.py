#!/usr/bin/env python3
"""check_md_links.py - check the relative links of Markdown files (task P15.D2).

usage: check_md_links.py FILE_OR_DIR [FILE_OR_DIR ...]

For every Markdown file given (directories are searched recursively for *.md), every inline link or image
`[text](target)` and every reference definition `[id]: target` outside fenced code blocks and code spans is checked:
  * http(s)://, mailto: and ftp:// targets are not checked (no network access in the tests);
  * a relative path must exist (relative to the file that contains the link);
  * a `#anchor` (alone or after a .md path) must match a heading of the target file, with GitHub's slug rules
    (lower case, punctuation removed except '-' and '_', spaces -> '-', duplicates get -1, -2, ...).
Exit status 0 if every link resolves, 1 otherwise (each broken link is printed as file:line: target: reason), 2 usage.
"""
import os
import re
import sys

LINK_RE = re.compile(r'!?\[(?:[^\[\]]|\[[^\]]*\])*\]\(\s*<?([^)\s>]+)>?(?:\s+"[^"]*")?\s*\)')
REFDEF_RE = re.compile(r'^\s{0,3}\[[^\]]+\]:\s*<?(\S+?)>?(?:\s+.*)?$')
FENCE_RE = re.compile(r'^\s{0,3}(```|~~~)')
CODESPAN_RE = re.compile(r'(`+)(.+?)\1')
HEADING_RE = re.compile(r'^\s{0,3}(#{1,6})\s+(.*?)\s*#*\s*$')
EXTERNAL = ('http://', 'https://', 'mailto:', 'ftp://')


def slugify(text):
    text = re.sub(r'<[^>]+>', '', text)                 # inline HTML
    text = re.sub(r'\[([^\]]*)\]\([^)]*\)', r'\1', text)  # links -> their text
    text = text.replace('`', '').strip().lower()
    text = re.sub(r'[^\w\- ]', '', text, flags=re.UNICODE)
    return text.replace(' ', '-')


_anchor_cache = {}


def anchors_of(path):
    if path in _anchor_cache:
        return _anchor_cache[path]
    slugs, seen, in_fence = set(), {}, False
    with open(path, encoding='utf-8') as f:
        for line in f:
            if FENCE_RE.match(line):
                in_fence = not in_fence
                continue
            if in_fence:
                continue
            m = HEADING_RE.match(line)
            if not m:
                for a in re.findall(r'<a\s+(?:name|id)="([^"]+)"', line):
                    slugs.add(a)
                continue
            s = slugify(m.group(2))
            n = seen.get(s, 0)
            seen[s] = n + 1
            slugs.add(s if n == 0 else f'{s}-{n}')
    _anchor_cache[path] = slugs
    return slugs


def check_file(path):
    errors, in_fence = [], False
    base = os.path.dirname(os.path.abspath(path))
    with open(path, encoding='utf-8') as f:
        lines = f.readlines()
    for no, line in enumerate(lines, 1):
        if FENCE_RE.match(line):
            in_fence = not in_fence
            continue
        if in_fence:
            continue
        text = CODESPAN_RE.sub('', line)
        targets = [m.group(1) for m in LINK_RE.finditer(text)]
        m = REFDEF_RE.match(text)
        if m:
            targets.append(m.group(1))
        for t in targets:
            if t.startswith(EXTERNAL):
                continue
            file_part, _, anchor = t.partition('#')
            target = os.path.abspath(path) if not file_part else os.path.normpath(os.path.join(base, file_part))
            if not os.path.exists(target):
                errors.append(f'{path}:{no}: {t}: file not found ({target})')
                continue
            if anchor:
                if not target.endswith('.md') or os.path.isdir(target):
                    continue
                if anchor not in anchors_of(target):
                    errors.append(f'{path}:{no}: {t}: no heading with anchor #{anchor} in {target}')
    return errors


def main(argv):
    if len(argv) < 2:
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        return 2
    files = []
    for a in argv[1:]:
        if os.path.isdir(a):
            for root, _, names in os.walk(a):
                files += [os.path.join(root, n) for n in sorted(names) if n.endswith('.md')]
        elif os.path.isfile(a):
            files.append(a)
        else:
            print(f'check_md_links: {a}: no such file or directory', file=sys.stderr)
            return 2
    errors = []
    for f in files:
        errors += check_file(f)
    for e in errors:
        print(e)
    print(f'check_md_links: {len(files)} file(s), {len(errors)} broken link(s)')
    return 1 if errors else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
