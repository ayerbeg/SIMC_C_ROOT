#!/usr/bin/env python3
"""build_physics_pdf.py - build the SIMC++ physics notes (docs/physics/latex/*.tex) into PDFs.

usage: build_physics_pdf.py --source-root DIR --build-dir DIR [--pdf-dir DIR] [--publish] [--note SOURCE ...]
       build_physics_pdf.py --source-root DIR --check-copies

The notes and their order are listed in docs/physics/latex/notes.json. 00_common is always built first: the reaction
notes refer to its labels (xr-hyper reads 00_common.aux from the build directory).

Steps (all intermediate files under --build-dir):
  1. lint every note source: a displayed equation must be a  feq  or  feqa  environment whose mandatory argument is a
     code citation (a file name such as event.f:518, radc.f:..., a .cpp/.h file or a FORTRAN_DIFFERENCES id);
     equation/align/gather/multline/\\[...\\]/$$ displays are refused (rule: every formula cites the code);
  2. simcnote_version.tex (git describe --always --tags and the commit date of the source tree);
  3. xelatex until the references are stable (at most 5 passes); fails on a LaTeX error, an undefined reference or
     citation, a missing file or glyph;
  4. the PDF is copied to --pdf-dir (default: the build directory). With --publish it is also copied to
     docs/physics/pdf/<output>.pdf and to every directory listed under "copies" in notes.json (for example
     src/physics/proton/): the copies are the same bytes as the master, written by this script only.
--check-copies only compares the committed master PDFs of docs/physics/pdf with their copies (CTest
docs_physics_pdf_copies) and fails if a copy is missing or differs.
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
DISPLAY_FORBIDDEN = [r'\\begin\{(equation|align|gather|multline|eqnarray|displaymath)\*?\}', r'\\\[', r'\$\$']
CITATION = re.compile(r'(\.(f|inc|cpp|h|theory|dat|inp)\b)|FORTRAN_DIFFERENCES|D-\d{3}|U-[A-Z]\d')


def run(cmd, cwd, env=None):
    return subprocess.run(cmd, cwd=cwd, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                          errors='replace')


def strip_comments(text):
    out = []
    for line in text.split('\n'):
        m = re.search(r'(?<!\\)%', line)
        out.append(line[:m.start()] if m else line)
    return '\n'.join(out)


def lint(path):
    """Return (number of cited displays, list of problems) for one note source."""
    text = strip_comments(open(path, encoding='utf-8').read())
    problems = []
    for pat in DISPLAY_FORBIDDEN:
        for m in re.finditer(pat, text):
            line = text.count('\n', 0, m.start()) + 1
            problems.append(f'{os.path.basename(path)}:{line}: display without code citation ({m.group(0)}); '
                            'use \\begin{feq}{file:line} or \\begin{feqa}{file:line}')
    n = 0
    for m in re.finditer(r'\\begin\{(feqa?)\}\{', text):
        # read the balanced argument
        i, depth = m.end(), 1
        while i < len(text) and depth:
            depth += {'{': 1, '}': -1}.get(text[i], 0)
            i += 1
        arg = text[m.end():i - 1]
        n += 1
        line = text.count('\n', 0, m.start()) + 1
        if not CITATION.search(arg):
            problems.append(f'{os.path.basename(path)}:{line}: the citation "{arg}" of a {m.group(1)} names no code file')
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


def check_copies(root, notes):
    bad = []
    for n in notes:
        master = os.path.join(root, 'docs', 'physics', 'pdf', n['output'] + '.pdf')
        if not os.path.exists(master):
            bad.append('missing master ' + os.path.relpath(master, root))
            continue
        for d in n.get('copies', []):
            copy = os.path.join(root, d, n['output'] + '.pdf')
            if not os.path.exists(copy):
                bad.append('missing copy ' + os.path.relpath(copy, root))
            elif not filecmp.cmp(master, copy, shallow=False):
                bad.append('copy differs from the master: ' + os.path.relpath(copy, root))
            readme = os.path.join(root, d, 'README.md')
            if not os.path.exists(readme) or n['output'] + '.pdf' not in open(readme, encoding='utf-8').read():
                bad.append(f'{os.path.relpath(readme, root)} does not name {n["output"]}.pdf')
    for b in bad:
        print('build_physics_pdf: FAIL: ' + b)
    if not bad:
        print(f'build_physics_pdf: {sum(len(n.get("copies", [])) for n in notes)} copies identical to their masters')
    return 1 if bad else 0


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument('--source-root', required=True)
    ap.add_argument('--build-dir')
    ap.add_argument('--pdf-dir')
    ap.add_argument('--publish', action='store_true',
                    help='copy the PDFs to docs/physics/pdf and to the "copies" directories of notes.json')
    ap.add_argument('--note', action='append', help='build only this note source (00_common.tex is always built)')
    ap.add_argument('--check-copies', action='store_true')
    a = ap.parse_args(argv[1:])
    root = os.path.abspath(a.source_root)
    ldir = os.path.join(root, 'docs', 'physics', 'latex')
    notes = json.load(open(os.path.join(ldir, 'notes.json'), encoding='utf-8'))['notes']
    if a.check_copies:
        return check_copies(root, notes)
    if not a.build_dir:
        ap.error('--build-dir is required')

    failures = []
    wanted = set(a.note or [n['source'] for n in notes]) | {'00_common.tex'}
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
    print(f'build_physics_pdf: lint ok, {neq} displayed equations, each with a code citation')

    xelatex = shutil.which('xelatex')
    if not xelatex:
        print('build_physics_pdf: xelatex not found: the physics notes are not built (SKIP)')
        return SKIP
    bdir = os.path.abspath(a.build_dir)
    os.makedirs(bdir, exist_ok=True)
    version, date = git_meta(root)
    with open(os.path.join(bdir, 'simcnote_version.tex'), 'w', encoding='utf-8') as f:
        f.write('%% written by docs/physics/build_physics_pdf.py\n\\renewcommand\\simcnoteversion{%s}'
                '\\renewcommand\\simcnotedate{%s}\n' % (version.replace('_', '\\_'), date))
    env = dict(os.environ)
    env['TEXINPUTS'] = os.pathsep.join([ldir, os.path.join(root, 'docs', 'common', 'latex') + '//',
                                        os.path.join(root, 'docs', 'common', 'figures'), bdir, ''])
    env['max_print_line'] = '1000'

    info = shutil.which('pdfinfo')
    for n in notes:
        if n['source'] not in wanted:
            continue
        job = os.path.splitext(n['source'])[0]
        src = os.path.join(ldir, n['source'])
        ok = True
        for p in range(1, 6):
            r = run([xelatex, '-interaction=nonstopmode', '-halt-on-error', '-jobname', job, src], bdir, env)
            errs, overfull, big, rerun = latex_problems(os.path.join(bdir, job + '.log'))
            if r.returncode != 0:
                errs = [e for e in errs if e.startswith('LaTeX error')] or errs or ['xelatex exit %d' % r.returncode]
                failures.append(f'{job}: pass {p}: ' + '; '.join(errs[:10]))
                ok = False
                break
            if p >= 2 and not rerun:
                break
        if not ok:
            continue
        if errs:
            failures.append(f'{job}: ' + '; '.join(sorted(set(errs))[:20]))
            continue
        if rerun:
            failures.append(f'{job}: references still changing after 5 passes')
        built = os.path.join(bdir, job + '.pdf')
        outs = [os.path.join(os.path.abspath(a.pdf_dir or bdir), n['output'] + '.pdf')]
        if a.publish:
            outs.append(os.path.join(root, 'docs', 'physics', 'pdf', n['output'] + '.pdf'))
            outs += [os.path.join(root, d, n['output'] + '.pdf') for d in n.get('copies', [])]
        for out in outs:
            os.makedirs(os.path.dirname(out), exist_ok=True)
            if os.path.abspath(out) != os.path.abspath(built):
                shutil.copyfile(built, out)
        pages = '?'
        if info:
            m = re.search(r'^Pages:\s+(\d+)', run([info, built], bdir).stdout, re.M)
            pages = m.group(1) if m else '?'
        where = ', '.join(os.path.relpath(o, root) if o.startswith(root) else o for o in outs)
        print(f'build_physics_pdf: {job}: {pages} pages, {p} xelatex passes, {overfull} overfull boxes ({big} > 3pt)'
              f' -> {where}')
    for f in failures:
        print('build_physics_pdf: FAIL: ' + f)
    return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
