#!/usr/bin/env python3
"""check_dev_manual.py - checks of the developer manual docs/manual/developer/ (task P15.D1).

usage: check_dev_manual.py SOURCE_DIR [developer|translation]

developer (default), docs/manual/developer/:
1. header coverage: every header under include/simc/ is named (as `simc/<path>.h`) somewhere in docs/manual/developer/*.md;
2. stage interfaces: every class declared in include/simc/stages/*.h is named in 03_architecture.md;
3. verbatim snippets (below).
translation, docs/translation/ (the translation report):
3. verbatim snippets: every fenced code block preceded by `<!-- from: <path> -->` must appear in <path> (relative to
   SOURCE_DIR) as a contiguous sequence of lines (leading/trailing whitespace of each line ignored);
4. the generated Fortran map (fortran_routine_map.md) has no routine with status "not cited";
5. the generated catalogue (06_catalogue.md) is up to date with docs/FORTRAN_DIFFERENCES.md (sections 6.1-6.3 are
   regenerated with tools/gen_catalogue.py and compared; 6.4 depends on project_control/DECISIONS.md and is not compared).
Exit 0 if everything holds, 1 otherwise (each problem printed), 2 usage.
"""
import glob
import os
import re
import sys


def main():
    if len(sys.argv) not in (2, 3) or (len(sys.argv) == 3 and sys.argv[2] not in ("developer", "translation")):
        print(__doc__)
        return 2
    src = sys.argv[1]
    mode = sys.argv[2] if len(sys.argv) == 3 else "developer"
    man = os.path.join(src, "docs", "manual", "developer") if mode == "developer" else os.path.join(src, "docs", "translation")
    texts = {p: open(p, encoding="utf-8").read() for p in sorted(glob.glob(os.path.join(man, "*.md")))}
    alltext = "\n".join(texts.values())
    problems = []

    inc = os.path.join(src, "include")
    if mode == "developer":
        check_headers(src, man, texts, alltext, problems)
    nsnip = check_snippets(src, texts, problems)
    print("verbatim snippets: %d" % nsnip)
    if mode == "translation":
        check_translation(src, man, problems)
    for p in problems:
        print("PROBLEM: " + p)
    print("FAIL" if problems else "PASS")
    return 1 if problems else 0


def check_headers(src, man, texts, alltext, problems):
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


def check_snippets(src, texts, problems):
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
    return nsnip


def check_translation(src, man, problems):
    # 4. generated map
    fmap = os.path.join(man, "fortran_routine_map.md")
    if os.path.exists(fmap):
        nc = [l for l in open(fmap, encoding="utf-8") if "| not cited |" in l]
        for l in nc:
            problems.append("routine without a C++ citation or reason: " + l.strip()[:120])
    else:
        problems.append("fortran_routine_map.md missing")
    # 5. catalogue up to date
    import subprocess
    import tempfile
    cat = os.path.join(man, "06_catalogue.md")
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "cat.md")
        r = subprocess.run([sys.executable, "-I", os.path.join(man, "tools", "gen_catalogue.py"), "--out", out],
                           capture_output=True, text=True)
        if r.returncode != 0 or not os.path.exists(out):
            problems.append("gen_catalogue.py failed: " + r.stderr.strip()[-300:])
        elif not os.path.exists(cat):
            problems.append("06_catalogue.md missing")
        else:
            cut = lambda t: t.split("## 6.4")[0]
            if cut(open(out, encoding="utf-8").read()) != cut(open(cat, encoding="utf-8").read()):
                problems.append("06_catalogue.md is out of date with docs/FORTRAN_DIFFERENCES.md: run "
                                "python3 -I docs/translation/tools/gen_catalogue.py")
            else:
                print("catalogue: up to date")


if __name__ == "__main__":
    sys.exit(main())
