#!/usr/bin/env python3
"""check_parameters.py - the user manual documents every input parameter and every command-line option (task P15.D2).

usage: check_parameters.py <inp2json> <simc++> <parameter reference .md> <manual dir>

1. Every input variable of SimcParameters (the 107 rows of ParameterTable(), printed by `inp2json --mapping`) must have
   a row in the parameter reference: a Markdown table line whose first cell is the Fortran name in backticks.
   Rows naming a '%' variable that does not exist, or a parameter with the wrong case, are errors too.
2. Every option printed by `simc++ --help` (--xxx and -x) must appear in backticks somewhere in the manual directory.
Exit status 0 if complete, 1 otherwise (the missing items are printed).
"""
import os
import re
import subprocess
import sys


def main(argv):
    if len(argv) != 5:
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        return 2
    inp2json, simcpp, ref, mandir = argv[1:]
    mapping = subprocess.run([inp2json, '--mapping'], capture_output=True, text=True, check=True).stdout
    names = []
    for line in mapping.splitlines():
        cells = [c.strip() for c in line.strip().strip('|').split('|')]
        if len(cells) > 2 and cells[0].isdigit():
            names.append(cells[1].strip('`'))
    if len(names) < 100:
        print(f'inp2json --mapping gave only {len(names)} parameters: cannot check')
        return 1
    documented = set()
    with open(ref, encoding='utf-8') as f:
        for line in f:
            m = re.match(r'^\|\s*`([^`]+)`\s*\|', line)
            if m:
                documented.add(m.group(1))
    lower = {n.lower(): n for n in names}
    missing = [n for n in names if n not in documented]
    # a backticked first cell with '%' (a Fortran structure name) or a case variant of a parameter must be exact
    unknown = sorted(d for d in documented if d not in names and ('%' in d or d.lower() in lower))
    bad = 0
    for n in missing:
        print(f'parameter not documented in {ref}: {n}')
        bad += 1
    for d in unknown:
        print(f'{ref} documents a parameter that SimcParameters does not have: {d}')
        bad += 1
    print(f'parameters: {len(names)} in SimcParameters, {len(names) - len(missing)} documented')

    helptext = subprocess.run([simcpp, '--help'], capture_output=True, text=True).stderr
    options = sorted(set(re.findall(r'(?<![\w-])(--[a-z][a-z-]*|-[a-z])\b', helptext)))
    text = ''
    for root, _, files in os.walk(mandir):
        for n in files:
            if n.endswith('.md'):
                with open(os.path.join(root, n), encoding='utf-8') as f:
                    text += f.read()
    for o in options:
        if not re.search(r'`[^`\n]*' + re.escape(o) + r'(?![\w-])[^`\n]*`', text):
            print(f'command-line option not documented in {mandir}: {o}')
            bad += 1
    print(f'command-line options: {len(options)} in simc++ --help')
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
