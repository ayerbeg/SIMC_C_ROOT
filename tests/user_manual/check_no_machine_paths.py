#!/usr/bin/env python3
"""check_no_machine_paths.py - the documentation contains no machine-specific paths (decision D-073).

usage: check_no_machine_paths.py FILE_OR_DIR [FILE_OR_DIR ...]

Every text file given (directories are searched recursively) must not contain "/home/" or "NewDisk". Documentation uses
paths relative to the source tree or the layout $SIMCPP_SRC, $SIMCPP_BUILD, $SIMCPP_PREFIX, $SIMC_FORTRAN, $SIMC_REFERENCE
(user manual, chapter 2.1). Binary files are skipped. Exit status 0 if clean, 1 otherwise (file:line: text), 2 usage.
"""
import os
import sys

PATTERNS = ('/home/', 'NewDisk')


def files_of(arg):
    if os.path.isdir(arg):
        for root, _, names in os.walk(arg):
            for n in sorted(names):
                yield os.path.join(root, n)
    else:
        yield arg


def main(argv):
    if len(argv) < 2:
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        return 2
    bad, nfiles = 0, 0
    for a in argv[1:]:
        if not os.path.exists(a):
            print(f'check_no_machine_paths: {a}: no such file or directory', file=sys.stderr)
            return 2
        for f in files_of(a):
            try:
                with open(f, encoding='utf-8') as fh:
                    lines = fh.readlines()
            except (UnicodeDecodeError, OSError):
                continue
            nfiles += 1
            for no, line in enumerate(lines, 1):
                if any(p in line for p in PATTERNS):
                    print(f'{f}:{no}: {line.rstrip()[:160]}')
                    bad += 1
    print(f'check_no_machine_paths: {nfiles} text file(s), {bad} line(s) with a machine-specific path')
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
