#!/usr/bin/env python3
"""run_doxygen.py - CTest docs_doxygen (task P15.D3): build the API reference and require 0 Doxygen warnings for include/.

usage: run_doxygen.py DOXYGEN DOXYFILE WARNING_LOG INCLUDE_DIR

Runs DOXYGEN on DOXYFILE (configured by docs/doxygen/CMakeLists.txt; its WARN_LOGFILE is WARNING_LOG), then reads
WARNING_LOG. Every warning located in a file under INCLUDE_DIR is printed and makes the test fail; the warnings of
other inputs (the Markdown manuals rendered as pages) are counted and reported but do not fail the test.
Exit status: 0 = no warning for include/, 1 = warnings (or Doxygen failed), 77 = Doxygen not installed (CTest: skipped),
2 = usage.
"""
import os
import shutil
import subprocess
import sys


def main():
    if len(sys.argv) != 5:
        print(__doc__)
        return 2
    exe, doxyfile, log, inc = sys.argv[1:]
    path = exe if os.path.isabs(exe) and os.access(exe, os.X_OK) else shutil.which(exe)
    if not path:
        print("doxygen not found (%s): docs_doxygen skipped; the coverage is still checked by docs_doc_coverage" % exe)
        return 77
    if os.path.exists(log):
        os.remove(log)
    r = subprocess.run([path, doxyfile], cwd=os.path.dirname(doxyfile), capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-3000:])
        print(r.stderr[-3000:])
        print("doxygen failed with exit status %d" % r.returncode)
        return 1
    lines = open(log, encoding="utf-8", errors="replace").read().splitlines() if os.path.exists(log) else []
    inc = os.path.realpath(inc)
    mine, other = [], 0
    for l in lines:
        if not l.strip() or l.startswith(" "):
            continue  # continuation lines of a multi-line warning
        f = l.split(":", 1)[0]
        if f and os.path.realpath(f).startswith(inc + os.sep):
            mine.append(l)
        else:
            other += 1
    for l in mine:
        print(l)
    print("doxygen %s: %d warnings for include/, %d for the other inputs (Markdown pages; not counted) -> %s"
          % (subprocess.run([path, "--version"], capture_output=True, text=True).stdout.strip(), len(mine), other,
             "PASS" if not mine else "FAIL"))
    print("HTML: %s" % os.path.join(os.path.dirname(doxyfile), "html", "index.html"))
    return 1 if mine else 0


if __name__ == "__main__":
    sys.exit(main())
