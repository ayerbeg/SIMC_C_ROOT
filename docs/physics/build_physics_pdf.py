#!/usr/bin/env python3
"""build_physics_pdf.py - build the SIMC++ physics notes and the bound "SIMC++ Physics Reference".

usage: build_physics_pdf.py --source-root DIR --build-dir DIR [--pdf-dir DIR] [--publish] [--note SOURCE ...] [--no-book]
       build_physics_pdf.py --source-root DIR --lint-only [--note SOURCE ...]
       build_physics_pdf.py --source-root DIR --lint-file FILE [--lint-file FILE ...]
       build_physics_pdf.py --source-root DIR --check-copies

Single source: docs/physics/latex/<NN>_<name>.tex (listed in docs/physics/latex/notes.json, built in that order;
00_common is always built first, the reaction notes read its labels through xr-hyper). The book
docs/physics/latex/physics_reference.tex inputs the same files (note 00 as Part I, the reactions as the chapters of
Part II, in the order of notes.json) through physics_reference_notes.tex, which this script writes.

Steps (all intermediate files under --build-dir):
  1. lint every note: every display is a feq/feqa with TWO citations, {C++ file: function}{Fortran file:lines} (the
     Fortran one may be "none" for a SIMC++-only formula); no plain equation/align/gather/multline/\\[ \\]/$$, no
     \\label (use \\lbl), no thebibliography/\\bibitem (use simcrefs), \\simcnotemeta after \\begin{document}, and no
     internal finding id (D-0xx, U-Xn, Pn.Fn, "Pn Qn") outside the simcfindings table and comments;
  2. simcnote_version.tex and simcpp_version.tex (git describe --always --tags, commit date);
  3. xelatex until the references are stable (at most 5 passes) for every note, then for the book; fails on a LaTeX
     error, an undefined reference or citation, a missing file or glyph;
  4. the PDFs are copied to --pdf-dir (default: the build directory). With --publish they are also copied to
     docs/physics/pdf/ (the masters) and every note to each source folder listed under "copies" in notes.json (the
     plug-in folders of its reaction): the copies are the same bytes as the master and are written only here.
--check-copies compares the committed masters with their copies (CTest docs_physics_pdf_copies), checks that each copy
folder has a README.md naming the PDF, and that the book is present.
Exit status: 0 ok, 1 failure, 77 xelatex not installed (CTest SKIP).
"""
import argparse
import filecmp
import json
import os
import re
import shutil
import subprocess
import sys

SKIP = 77
BOOK_SOURCE = 'physics_reference.tex'
DISPLAY_FORBIDDEN = [r'\\begin\{(equation|align|gather|multline|eqnarray|displaymath)\*?\}', r'\\\[', r'\$\$']
CXX_CITATION = re.compile(r'\.(cpp|h)\b')
F_CITATION = re.compile(r'\.(f|inc)\b|\bnone\b')
FINDING_ID = re.compile(r'\bD-\d{3}\b|\bU-(?:A|B|D|DEV)\d+\b|\bP\d+(?:\.\w+)?\s+Q\d+\b|\bP\d+\.F\d+\b')


def run(cmd, cwd, env=None):
    return subprocess.run(cmd, cwd=cwd, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                          errors='replace')


def strip_comments(text):
    out = []
    for line in text.split('\n'):
        m = re.search(r'(?<!\\)%', line)
        out.append(line[:m.start()] if m else line)
    return '\n'.join(out)


def balanced(text, i):
    """text[i] is just after an opening brace: return (content, index after the closing brace)."""
    depth, j = 1, i
    while j < len(text) and depth:
        depth += {'{': 1, '}': -1}.get(text[j], 0)
        j += 1
    return text[i:j - 1], j


def lint(path):
    """Return (number of cited displays, list of problems) for one note source."""
    name = os.path.basename(path)
    text = strip_comments(open(path, encoding='utf-8').read())
    problems = []

    def where(pos):
        return f'{name}:{text.count(chr(10), 0, pos) + 1}'
    for pat in DISPLAY_FORBIDDEN:
        for m in re.finditer(pat, text):
            problems.append(f'{where(m.start())}: display without code citations ({m.group(0)}); use '
                            '\\begin{feq}{C++ file: function}{Fortran file:lines}')
    n = 0
    for m in re.finditer(r'\\begin\{(feqa?)\}', text):
        n += 1
        i = m.end()
        args = []
        for _ in range(2):
            k = i
            while k < len(text) and text[k] in ' \t\n':
                k += 1
            if k < len(text) and text[k] == '{':
                arg, i = balanced(text, k + 1)
                args.append(arg)
            else:
                break
        if len(args) < 2:
            problems.append(f'{where(m.start())}: {m.group(1)} needs two citations {{C++}}{{Fortran}}')
            continue
        if not CXX_CITATION.search(args[0]):
            problems.append(f'{where(m.start())}: first citation "{args[0]}" names no C++ file (.cpp/.h)')
        if not F_CITATION.search(args[1]):
            problems.append(f'{where(m.start())}: second citation "{args[1]}" names no Fortran file (or "none")')
    for m in re.finditer(r'\\label\{', text):
        problems.append(f'{where(m.start())}: use \\lbl{{...}} instead of \\label (labels are prefixed per note)')
    for m in re.finditer(r'\\begin\{thebibliography\}|\\bibitem', text):
        problems.append(f'{where(m.start())}: use the simcrefs environment and \\simcref/\\citer')
    b, meta = text.find('\\begin{document}'), text.find('\\simcnotemeta')
    if b < 0 or meta < b:
        problems.append(f'{name}: \\simcnotemeta must come after \\begin{{document}} (the book skips the preamble)')
    body = re.sub(r'\\begin\{simcfindings\}.*?\\end\{simcfindings\}', '', text, flags=re.S)
    for m in FINDING_ID.finditer(body):
        problems.append(f'{where(m.start())}: finding id "{m.group(0)}" outside the findings table (describe the '
                        'fact in words; the id belongs to the simcfindings table of the last section)')
    return n, problems


def latex_problems(logfile):
    text = open(logfile, encoding='utf-8', errors='replace').read()
    errors = []
    for m in re.finditer(r'^! (.*)$', text, re.M):
        errors.append('LaTeX error: ' + m.group(1))
    for m in re.finditer(r"(LaTeX|Package hyperref) Warning: (Reference|Hyper reference|Citation) `([^']*)' on page "
                         r"(\S+) undefined", text):
        errors.append(f'undefined {m.group(2).lower()} {m.group(3)} (page {m.group(4)})')
    for m in re.finditer(r"LaTeX Warning: File `([^']*)' not found", text):
        errors.append('file not found: ' + m.group(1))
    for m in re.finditer(r'Package xr(-hyper)? Warning: No file (\S+)', text):
        errors.append('xr: no file ' + m.group(2))
    for m in sorted(set(re.findall(r'Missing character: There is no (.) \(U\+([0-9A-F]+)\) in font ([^!\n]*)', text))):
        errors.append(f'missing glyph U+{m[1]} in font {m[2]}')
    overfull = [float(x) for x in re.findall(r'^Overfull \\hbox \((\d+\.\d+)pt', text, re.M)]
    rerun = bool(re.search(r'Rerun to get|Label\(s\) may have changed|Package rerunfilecheck Warning', text))
    return errors, len(overfull), len([x for x in overfull if x > 3.0]), rerun


def git_meta(root):
    def g(*args):
        try:
            return subprocess.run(['git', '-C', root] + list(args), capture_output=True, text=True,
                                  check=True).stdout.strip()
        except (OSError, subprocess.CalledProcessError):
            return ''
    return (g('describe', '--always', '--tags') or 'unknown'), (g('log', '-1', '--format=%cs') or '')


def check_copies(root, manifest):
    notes, book = manifest['notes'], manifest.get('book', {})
    bad, ncopies = [], 0
    pdfdir = os.path.join(root, 'docs', 'physics', 'pdf')
    for n in notes:
        master = os.path.join(pdfdir, n['output'] + '.pdf')
        if not os.path.exists(master):
            bad.append('missing master ' + os.path.relpath(master, root))
            continue
        for d in n.get('copies', []):
            ncopies += 1
            copy = os.path.join(root, d, n['output'] + '.pdf')
            if not os.path.exists(copy):
                bad.append('missing copy ' + os.path.relpath(copy, root))
            elif not filecmp.cmp(master, copy, shallow=False):
                bad.append('copy differs from the master: ' + os.path.relpath(copy, root))
            readme = os.path.join(root, d, 'README.md')
            if not os.path.exists(readme) or n['output'] + '.pdf' not in open(readme, encoding='utf-8').read():
                bad.append(f'{os.path.relpath(readme, root)} does not name {n["output"]}.pdf')
    if book and not os.path.exists(os.path.join(pdfdir, book['output'] + '.pdf')):
        bad.append(f'missing book docs/physics/pdf/{book["output"]}.pdf')
    for b in bad:
        print('build_physics_pdf: FAIL: ' + b)
    if not bad:
        print(f'build_physics_pdf: {ncopies} copies identical to their masters; book present')
    return 1 if bad else 0


def compile_tex(xelatex, src, job, bdir, env):
    """xelatex until stable. Return (ok, errors, overfull, big, passes)."""
    for p in range(1, 6):
        r = run([xelatex, '-interaction=nonstopmode', '-halt-on-error', '-jobname', job, src], bdir, env)
        errs, overfull, big, rerun = latex_problems(os.path.join(bdir, job + '.log'))
        if r.returncode != 0:
            errs = [e for e in errs if e.startswith('LaTeX error')] or errs or ['xelatex exit %d' % r.returncode]
            return False, [f'pass {p}: ' + '; '.join(errs[:10])], overfull, big, p
        if p >= 2 and not rerun:
            break
    if rerun:
        errs = errs + ['references still changing after 5 passes']
    return not errs, errs, overfull, big, p


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument('--source-root', required=True)
    ap.add_argument('--build-dir')
    ap.add_argument('--pdf-dir')
    ap.add_argument('--publish', action='store_true',
                    help='copy the PDFs to docs/physics/pdf and to the "copies" folders of notes.json')
    ap.add_argument('--note', action='append', help='build only this note source (00_common.tex is always built); '
                    'implies --no-book')
    ap.add_argument('--no-book', action='store_true', help='do not build the bound book')
    ap.add_argument('--lint-only', action='store_true')
    ap.add_argument('--lint-file', action='append', help='lint this file (e.g. the template) and exit')
    ap.add_argument('--check-copies', action='store_true')
    a = ap.parse_args(argv[1:])
    root = os.path.abspath(a.source_root)
    ldir = os.path.join(root, 'docs', 'physics', 'latex')
    manifest = json.load(open(os.path.join(ldir, 'notes.json'), encoding='utf-8'))
    notes = manifest['notes']
    book = manifest.get('book')
    if a.check_copies:
        return check_copies(root, manifest)
    if a.lint_file:
        bad = []
        for f in a.lint_file:
            k, probs = lint(f if os.path.isabs(f) else os.path.join(ldir, f))
            bad += probs
            print(f'build_physics_pdf: {os.path.basename(f)}: {k} displays')
        for b in bad:
            print('build_physics_pdf: FAIL: ' + b)
        return 1 if bad else 0

    failures = []
    wanted = set(a.note or [n['source'] for n in notes]) | {'00_common.tex'}
    do_book = bool(book) and not a.note and not a.no_book
    # 1. lint
    neq = 0
    for n in notes:
        if n['source'] not in wanted:
            continue
        k, probs = lint(os.path.join(ldir, n['source']))
        neq += k
        failures += probs
    if failures:
        for f in failures:
            print('build_physics_pdf: FAIL: ' + f)
        return 1
    print(f'build_physics_pdf: lint ok, {neq} displayed equations, each with a C++ and a Fortran citation')
    if a.lint_only:
        return 0
    if not a.build_dir:
        ap.error('--build-dir is required')

    xelatex = shutil.which('xelatex')
    if not xelatex:
        print('build_physics_pdf: xelatex not found: the physics notes are not built (SKIP)')
        return SKIP
    bdir = os.path.abspath(a.build_dir)
    os.makedirs(bdir, exist_ok=True)
    version, date = git_meta(root)
    v = version.replace('_', '\\_')
    with open(os.path.join(bdir, 'simcnote_version.tex'), 'w', encoding='utf-8') as f:
        f.write('%% written by docs/physics/build_physics_pdf.py\n\\renewcommand\\simcnoteversion{%s}'
                '\\renewcommand\\simcnotedate{%s}\n' % (v, date))
    with open(os.path.join(bdir, 'simcpp_version.tex'), 'w', encoding='utf-8') as f:
        f.write('%% written by docs/physics/build_physics_pdf.py\n\\renewcommand\\simcbuildversion{%s}'
                '\\renewcommand\\simcbuilddate{%s}\n' % (v, date))
    env = dict(os.environ)
    env['TEXINPUTS'] = os.pathsep.join([ldir, os.path.join(root, 'docs', 'common', 'latex') + '//',
                                        os.path.join(root, 'docs', 'common', 'figures'), bdir, ''])
    env['max_print_line'] = '1000'
    info = shutil.which('pdfinfo')

    def pages_of(pdf):
        if not info:
            return '?'
        m = re.search(r'^Pages:\s+(\d+)', run([info, pdf], bdir).stdout, re.M)
        return m.group(1) if m else '?'

    def deliver(built, output, copies):
        outs = [os.path.join(os.path.abspath(a.pdf_dir or bdir), output + '.pdf')]
        if a.publish:
            outs.append(os.path.join(root, 'docs', 'physics', 'pdf', output + '.pdf'))
            outs += [os.path.join(root, d, output + '.pdf') for d in copies]
        for out in outs:
            os.makedirs(os.path.dirname(out), exist_ok=True)
            if os.path.abspath(out) != os.path.abspath(built):
                shutil.copyfile(built, out)
        return ', '.join(os.path.relpath(o, root) if o.startswith(root) else o for o in outs)

    # 3. the notes
    for n in notes:
        if n['source'] not in wanted:
            continue
        job = os.path.splitext(n['source'])[0]
        ok, errs, overfull, big, p = compile_tex(xelatex, os.path.join(ldir, n['source']), job, bdir, env)
        if not ok:
            failures.append(f'{job}: ' + '; '.join(sorted(set(errs))[:20]))
            continue
        built = os.path.join(bdir, job + '.pdf')
        where = deliver(built, n['output'], n.get('copies', []))
        print(f'build_physics_pdf: {job}: {pages_of(built)} pages, {p} xelatex passes, {overfull} overfull boxes '
              f'({big} > 3pt) -> {where}')

    # 3b. the book
    if do_book:
        with open(os.path.join(bdir, 'physics_reference_notes.tex'), 'w', encoding='utf-8') as f:
            f.write('% written by docs/physics/build_physics_pdf.py from notes.json: the reaction chapters\n')
            for n in notes:
                if n['source'] != '00_common.tex':
                    f.write('\\simcinputnote{%s}\n' % os.path.splitext(n['source'])[0])
        job = os.path.splitext(BOOK_SOURCE)[0]
        ok, errs, overfull, big, p = compile_tex(xelatex, os.path.join(ldir, BOOK_SOURCE), job, bdir, env)
        if not ok:
            failures.append(f'{job}: ' + '; '.join(sorted(set(errs))[:20]))
        else:
            built = os.path.join(bdir, job + '.pdf')
            where = deliver(built, book['output'], [])
            print(f'build_physics_pdf: {job}: {pages_of(built)} pages, {p} xelatex passes, {overfull} overfull boxes '
                  f'({big} > 3pt) -> {where}')

    for f in failures:
        print('build_physics_pdf: FAIL: ' + f)
    return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
