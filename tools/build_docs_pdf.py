#!/usr/bin/env python3
"""build_docs_pdf.py - build the PDF editions of the SIMC++ books (task P15.D5).

usage: build_docs_pdf.py --source-root DIR --build-dir DIR [--pdf-dir DIR] [--book NAME ...] [--png] [--figures-only]

Steps (all products under --build-dir, normally <cmake build>/docs):
  1. every TikZ figure docs/latex/figures/<name>.tex is compiled alone (xelatex) into figures/<name>.pdf;
     with --png the figures are also exported to docs/figures/<name>.png (pdftoppm, 150 dpi) for the Markdown view
     (this writes into the source tree: run it by hand after changing a figure and commit the PNGs);
  2. tools/md2latex.py converts the Markdown chapters into <book>_body.tex / <book>_meta.tex (books.json);
  3. xelatex, makeindex (docs/latex/simcpp.ist), xelatex, xelatex (until no "Rerun" request, at most 5 passes) for
     docs/latex/<book>.tex; the PDF is copied to <pdf-dir>/<output>.pdf (books.json "output"; default --build-dir).
Fails (exit 1) on a converter error, a LaTeX error, an undefined reference or citation, or a missing figure;
exit 77 if xelatex is not installed (CTest SKIP). Overfull boxes are counted and reported, not fatal.
"""
import argparse
import glob
import json
import os
import re
import shutil
import subprocess
import sys

SKIP = 77


def run(cmd, cwd, env, log=None):
    r = subprocess.run(cmd, cwd=cwd, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                       errors='replace')
    if log:
        with open(log, 'w', encoding='utf-8') as f:
            f.write(r.stdout)
    return r


def latex_problems(logfile):
    text = open(logfile, encoding='utf-8', errors='replace').read()
    errors = []
    for m in re.finditer(r'^! (.*)$', text, re.M):
        errors.append('LaTeX error: ' + m.group(1))
    for m in re.finditer(r"(LaTeX|Package hyperref) Warning: (Reference|Hyper reference|Citation) `([^']*)' on page "
                         r"(\S+) undefined", text):
        errors.append(f'undefined {m.group(2).lower()} {m.group(3)} (page {m.group(4)})')
    for m in re.finditer(r'LaTeX Warning: File `([^\']*)\' not found', text):
        errors.append('file not found: ' + m.group(1))
    if re.search(r'Missing character: There is no', text):
        for m in sorted(set(re.findall(r'Missing character: There is no (.) \(U\+([0-9A-F]+)\) in font ([^!\n]*)', text))):
            errors.append(f'missing glyph U+{m[1]} in font {m[2]}')
    overfull = len(re.findall(r'^Overfull \\hbox \((\d+\.\d+)pt', text, re.M))
    big = [float(x) for x in re.findall(r'^Overfull \\hbox \((\d+\.\d+)pt', text, re.M) if float(x) > 3.0]
    rerun = bool(re.search(r'Rerun to get|Label\(s\) may have changed|Package rerunfilecheck Warning', text))
    return errors, overfull, len(big), rerun


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument('--source-root', required=True)
    ap.add_argument('--build-dir', required=True)
    ap.add_argument('--pdf-dir', help='where the final PDFs are copied (default: --build-dir)')
    ap.add_argument('--book', action='append')
    ap.add_argument('--png', action='store_true', help='also export the figures to docs/figures/*.png')
    ap.add_argument('--figures-only', action='store_true')
    a = ap.parse_args(argv[1:])
    xelatex = shutil.which('xelatex')
    if not xelatex:
        print('build_docs_pdf: xelatex not found: the PDF manuals are not built (SKIP)')
        return SKIP
    makeindex = shutil.which('makeindex')
    root = os.path.abspath(a.source_root)
    bdir = os.path.abspath(a.build_dir)
    latexdir = os.path.join(root, 'docs', 'latex')
    figsrc = os.path.join(latexdir, 'figures')
    figout = os.path.join(bdir, 'figures')
    os.makedirs(figout, exist_ok=True)
    env = dict(os.environ)
    env['TEXINPUTS'] = os.pathsep.join([latexdir + '//', os.path.join(root, 'docs', 'figures'), bdir, ''])
    env['max_print_line'] = '1000'
    failures = []

    # 1. figures
    for tex in sorted(glob.glob(os.path.join(figsrc, '*.tex'))):
        name = os.path.splitext(os.path.basename(tex))[0]
        r = run([xelatex, '-interaction=nonstopmode', '-halt-on-error', '-output-directory', figout, tex], figout, env)
        log = os.path.join(figout, name + '.log')
        errs = latex_problems(log)[0] if os.path.exists(log) else ['no log']
        if r.returncode != 0 or errs:
            failures.append(f'figure {name}: ' + '; '.join(errs or ['xelatex exit %d' % r.returncode]))
            continue
        if a.png:
            pdftoppm = shutil.which('pdftoppm')
            if not pdftoppm:
                failures.append('pdftoppm not found: cannot export the PNG figures')
            else:
                dst = os.path.join(root, 'docs', 'figures', name)
                run([pdftoppm, '-png', '-r', '150', '-singlefile', os.path.join(figout, name + '.pdf'), dst], figout, env)
                print(f'build_docs_pdf: exported docs/figures/{name}.png')
    if a.figures_only:
        for f in failures:
            print(f)
        return 1 if failures else 0

    # 2. conversion
    with open(os.path.join(latexdir, 'books.json'), encoding='utf-8') as f:
        books = json.load(f)
    names = a.book or list(books)
    conv = [sys.executable, '-I', os.path.join(root, 'tools', 'md2latex.py'), '--source-root', root,
            '--outdir', bdir, '--figdir', figout]
    for n in names:
        conv += ['--book', n]
    r = run(conv, bdir, env)
    print(r.stdout.rstrip())
    if r.returncode != 0:
        print('build_docs_pdf: converter errors (above)')
        return 1

    # 3. LaTeX
    for n in names:
        job = n
        src = os.path.join(latexdir, n + '.tex')
        ok, overfull, big = True, 0, 0
        for p in range(1, 6):
            r = run([xelatex, '-interaction=nonstopmode', '-halt-on-error', '-jobname', job, src], bdir, env)
            log = os.path.join(bdir, job + '.log')
            errs, overfull, big, rerun = latex_problems(log)
            if r.returncode != 0 and not errs:
                errs = ['xelatex exit status %d' % r.returncode]
            if r.returncode != 0:
                errs = [e for e in errs if e.startswith('LaTeX error')] or errs
                failures.append(f'{n}: pass {p}: ' + '; '.join(errs[:10]))
                ok = False
                break
            if p == 1 and makeindex:
                ri = run([makeindex, '-q', '-s', os.path.join(latexdir, 'simcpp.ist'), job + '.idx'], bdir, env)
                if ri.returncode != 0:
                    failures.append(f'{n}: makeindex failed: {ri.stdout[-500:]}')
            if p >= 3 and not rerun:
                break
        if not ok:
            continue
        if errs:
            failures.append(f'{n}: ' + '; '.join(sorted(set(errs))[:20]) + (' ...' if len(set(errs)) > 20 else ''))
        if rerun:
            failures.append(f'{n}: references still changing after 5 passes')
        os.makedirs(a.pdf_dir or bdir, exist_ok=True)
        out = os.path.join(os.path.abspath(a.pdf_dir or bdir), books[n]['output'] + '.pdf')
        shutil.copyfile(os.path.join(bdir, job + '.pdf'), out)
        pages = '?'
        info = shutil.which('pdfinfo')
        if info:
            m = re.search(r'^Pages:\s+(\d+)', run([info, out], bdir, env).stdout, re.M)
            pages = m.group(1) if m else '?'
        print(f'build_docs_pdf: {out}: {pages} pages, {p} xelatex passes, {overfull} overfull boxes ({big} > 3pt)')
    for f in failures:
        print('build_docs_pdf: FAIL: ' + f)
    return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
