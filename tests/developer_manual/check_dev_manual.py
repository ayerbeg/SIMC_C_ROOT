#!/usr/bin/env python3
"""check_dev_manual.py - checks of the developer manual docs/manual/developer/ (task P15.D1).

usage: check_dev_manual.py SOURCE_DIR

1. header coverage: every header under include/simc/ is named (as `simc/<path>.h`) somewhere in docs/manual/developer/*.md;
2. stage interfaces: every class declared in include/simc/stages/*.h is named in 03_architecture.md;
3. verbatim snippets: every fenced code block preceded by `<!-- from: <path> -->` must appear in <path> (relative to
   SOURCE_DIR) as a contiguous sequence of lines (leading/trailing whitespace of each line ignored);
4. the generated Fortran map (fortran_routine_map.md) has no routine with status "not cited".
Exit 0 if everything holds, 1 otherwise (each problem printed), 2 usage.
"""
import glob
import os
import re
import sys


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    src = sys.argv[1]
    man = os.path.join(src, "docs", "manual", "developer")
    texts = {p: open(p, encoding="utf-8").read() for p in sorted(glob.glob(os.path.join(man, "*.md")))}
    alltext = "\n".join(texts.values())
    problems = []

    # 1. headers
    inc = os.path.join(src, "include")
    headers = sorted(os.path.relpath(p, inc) for p in glob.glob(os.path.join(inc, "simc", "**", "*.h"), recursive=True))
    missing = [h for h in headers if h not in alltext]
    for h in missing:
        problems.append("header not mentioned in the developer manual: " + h)
    print("headers: %d, mentioned: %d" % (len(headers), len(headers) - len(missing)))

    # 2. stage interfaces
    arch = texts.get(os.path.join(man, "03_architecture.md"), "")
    classes = []
    for p in sorted(glob.glob(os.path.join(inc, "simc", "stages", "*.h"))):
        for m in re.finditer(r"^(?:class|struct)\s+(\w+)\s*(?::|\{)", open(p).read(), re.M):
            classes.append(m.group(1))
    for c in classes:
        if not re.search(r"\b" + re.escape(c) + r"\b", arch):
            problems.append("stage class not described in 03_architecture.md: " + c)
    print("stage-interface classes: %d" % len(classes))

    # 3. snippets
    nsnip = 0
    for path, text in texts.items():
        lines = text.split("\n")
        i = 0
        while i < len(lines):
            m = re.match(r"\s*<!--\s*from:\s*(\S+)\s*-->", lines[i])
            if m:
                j = i + 1
                while j < len(lines) and not lines[j].strip().startswith("```"):
                    j += 1
                k = j + 1
                while k < len(lines) and not lines[k].strip().startswith("```"):
                    k += 1
                block = [l.strip() for l in lines[j + 1:k]]
                target = os.path.join(src, m.group(1))
                nsnip += 1
                if not os.path.exists(target):
                    problems.append("%s:%d: snippet source missing: %s" % (os.path.basename(path), i + 1, m.group(1)))
                else:
                    flines = [l.strip() for l in open(target, encoding="utf-8").read().split("\n")]
                    n = len(block)
                    if not any(flines[a:a + n] == block for a in range(len(flines) - n + 1)):
                        problems.append("%s:%d: snippet is not a verbatim copy of %s" %
                                        (os.path.basename(path), i + 1, m.group(1)))
                i = k
            i += 1
    print("verbatim snippets: %d" % nsnip)

    # 4. generated map
    fmap = os.path.join(man, "fortran_routine_map.md")
    if os.path.exists(fmap):
        nc = [l for l in open(fmap, encoding="utf-8") if "| not cited |" in l]
        for l in nc:
            problems.append("routine without a C++ citation or reason: " + l.strip()[:120])
    else:
        problems.append("fortran_routine_map.md missing")

    for p in problems:
        print("PROBLEM: " + p)
    print("FAIL" if problems else "PASS")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
