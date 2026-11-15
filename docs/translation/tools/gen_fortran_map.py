#!/usr/bin/env python3
"""Generate docs/manual/developer/fortran_routine_map.md: every routine of the Fortran SIMC program -> the C++ files that
cite it (`// Fortran: file.f:routine ...` comments), plus the L1 reference tables that exist for it.

Usage: gen_fortran_map.py [--fortran DIR] [--reference DIR] [--out FILE]
Defaults: /home/ayerbe/NewDisk/simc_gfortran, /home/ayerbe/NewDisk/reference, docs/manual/developer/fortran_routine_map.md
(relative to the SIMC++ source tree, found from this script's location). The Fortran tree is only read.

The routine list is the set of SUBROUTINE/FUNCTION/PROGRAM statements of the files the Fortran Makefile links into
`simc` (my_objs + simc.f + the module files). A C++ file "cites" routine R of file F when one of its lines contains
`F:R` (F = file name without directory, case-insensitive), or contains both the file name F and the word R. The status
column is "ported" when at least one C++ file cites the routine, else the curated reason of NOT_PORTED below (or
"not cited": a routine without a citation and without a reason, which a reviewer should look at).
"""
import argparse
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.normpath(os.path.join(HERE, "..", "..", "..", ".."))

# Routines that have no C++ counterpart on purpose (checked by hand against the code, task P15.D1).
NOT_PORTED = {
    "ntupleinit": "HBOOK/.bin ntuple not written; the tree is written directly by RootTreeSink + NtupleLayout (D-002)",
    "ntupleclose": "idem (RootTreeSink::endRun closes the TFile)",
    "ntuple_close": "idem",
    "ntuple_init": "idem",
    "sgrnd": "RanluxRandom::seed (RLUXGO); default mt64 seeded by Mt64Random::seed",
    "grnd": "IRandom::uniform (RanluxRandom: grnd = ranlux single value; Mt64Random)",
    "save_random_state": "not ported: random_state_file is ignored with a NOTE (RANLUX state not portable, D-002)",
    "start_random_state": "idem",
    "abend": "CERNLIB error stop: C++ exceptions (core/Errors.h)",
    "kerset": "CERNLIB error-message table (used by fint/lfit error paths only): not needed",
    "kermtr": "idem",
    "results_ntu_write1": "never called in simc.f (dead code)",
    "decay": "event.f:DECAY is never called (dead code); the decays are in transp.f/project.f (optics::CosyEngine)",
    "start_file_random_state": "not ported as such: random_state_file is ignored with a NOTE (D-002); RanluxRandom "
                               "can start from a saved state given in JSON (core/random/RanluxRandom.h)",
    "trgxtrack": "not called by SIMC (trg_track.f library routine; SIMC uses trgTrackToPlane/trgRK4)",
    "trgrk4bdl": "only called by trgXTrack (not reached)",
    "print_coord1": "debug print, every call commented out in trg_track.f",
    "print_coord2": "idem",
    "print_coord3": "idem",
    "nextun": "free Fortran unit finder of the CTEQ5 table reader; C++ uses std::ifstream",
}
# F1F2IN21_v1.0.f: semi_physics.f:peepiX calls F1F2IN21 with A = 1 only (include/simc/physics/sidis/F1F2In21.h).
for _r in ("f1f2qe21 formfacts rescsd smearsub smearsublowq offshellqe getfy getfyoff qenuc21off mec2021 gsmearing "
           "sfcross csfitcomp vcoul dsimps cdbonn av18 wjc1 wjc2 pinterp hunt polintnum fnp_nmc").split():
    NOT_PORTED[_r] = "F1F2IN21 deuteron/nucleus/QE branch, not reached by SIMC (peepiX calls F1F2IN21 with A = 1; " \
                     "physics/sidis/F1F2In21.h)"
NOT_PORTED["sqesub"] = "called once by peepiX but only fills deuteron tables never used with A = 1: not ported " \
                       "(physics/sidis/F1F2In21.h, FORTRAN_DIFFERENCES.md)"
# L1 tables whose name starts with a routine name but belong to another routine
TABLE_EXCLUDE = {"decay": "*", "inter": {"inter_prime"}}

LINKED_EXTRA = ["simc.f", "modules.f", "histograms_module.f"]
DIRS = dict(A="shared/", H="hms/", S="sos/", R="hrsr/", L="hrsl/", SH="shms/", T="cteq5/", C="calo/", CH="cern/",
            D="fdss/")
DECL = re.compile(r"^[ \t]+(?:(?:real|integer|logical|double\s+precision|complex|character)(?:\*\d+)?\s+)?"
                  r"(subroutine|function|program|block\s*data)\s+(\w+)", re.I)


def fortran_files(fdir):
    mk = open(os.path.join(fdir, "Makefile")).read()
    files = list(LINKED_EXTRA)
    for m in re.finditer(r"^OBJ\w*\s*=\s*(.*)$", mk, re.M):
        for o in m.group(1).split():
            o = re.sub(r"\$\((\w+)\)", lambda x: DIRS[x.group(1)], o)
            files.append(o[:-2] + ".f")
    return files


def routines(path):
    lines = open(path, errors="replace").read().split("\n")
    out = []
    for i, l in enumerate(lines, 1):
        if not l or l[0] in "cC*!":
            continue
        m = DECL.match(l.split("!")[0])
        if m:
            out.append([m.group(2), i, None])
    for k, r in enumerate(out):
        r[2] = (out[k + 1][1] - 1) if k + 1 < len(out) else len(lines)
        # trim trailing blank / comment lines
        while r[2] > r[1] and (not lines[r[2] - 1].strip() or lines[r[2] - 1][0] in "cC*!"):
            r[2] -= 1
    return out


def cpp_lines():
    res = []
    for top in ("src", "include"):
        for d, _, fs in os.walk(os.path.join(SRC, top)):
            for f in sorted(fs):
                if f.endswith((".cpp", ".h")):
                    p = os.path.join(d, f)
                    rel = os.path.relpath(p, SRC)
                    res.append((rel, open(p, errors="replace").read().lower().split("\n")))
    return sorted(res)


def tables(refdir):
    t = []
    fdir = os.path.join(refdir, "functions")
    if not os.path.isdir(fdir):
        return t
    for d, ds, fs in os.walk(fdir):
        if "input.txt" in fs or "output.txt" in fs or "README.txt" in fs:
            t.append(os.path.relpath(d, fdir))
    return sorted(t)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fortran", default="/home/ayerbe/NewDisk/simc_gfortran")
    ap.add_argument("--reference", default="/home/ayerbe/NewDisk/reference")
    ap.add_argument("--out", default=os.path.join(SRC, "docs", "manual", "developer", "fortran_routine_map.md"))
    a = ap.parse_args()

    cpp = cpp_lines()
    tabs = tables(a.reference)
    rows, nport, nnot = [], 0, 0
    for f in fortran_files(a.fortran):
        p = os.path.join(a.fortran, f)
        if not os.path.exists(p):
            continue
        base = os.path.basename(f).lower()
        for name, l0, l1 in routines(p):
            r = name.lower()
            word = re.compile(r"\b" + re.escape(r) + r"\b")
            cites = []
            for rel, lines in cpp:
                if any((base + ":" + r) in l or (base in l and word.search(l)) for l in lines):
                    cites.append(rel)
            ex = TABLE_EXCLUDE.get(r, set())
            tl = [t for t in tabs if ex != "*" and os.path.basename(t) not in ex and
                  (os.path.basename(t) == r or os.path.basename(t).startswith(r + "_"))]
            if r in NOT_PORTED:   # curated: wins over a citation (a comment may cite a routine to say it is not ported)
                status = NOT_PORTED[r]
                nnot += 1
            elif cites:
                status = "ported"
                nport += 1
            else:
                status = "not cited"
                nnot += 1
            shown = cites[:4] + (["(+%d more)" % (len(cites) - 4)] if len(cites) > 4 else [])
            tshown = tl[:4] + (["(+%d more)" % (len(tl) - 4)] if len(tl) > 4 else [])
            rows.append("| %s | %s | %d-%d | %s | %s | %s |" % (
                f, name, l0, l1, "<br>".join("`%s`" % c if not c.startswith("(") else c for c in shown) or "-",
                status, ", ".join("`%s`" % t if not t.startswith("(") else t for t in tshown) or "-"))
    with open(a.out, "w") as o:
        o.write("# Fortran routine → C++ map (generated)\n\n")
        o.write("Generated by `docs/manual/developer/tools/gen_fortran_map.py` from the Fortran sources linked into "
                "`simc` (Makefile `my_objs` + `simc.f` + modules) and the `// Fortran: file.f:routine` citations in "
                "`src/` and `include/`. Do not edit by hand: re-run the script after a port "
                "(see [chapter 4](04_fortran_map.md)).\n\n")
        o.write("Columns: Fortran file, routine, line range of the routine in the file, the C++ files that cite it "
                "(first four), status, and the L1 reference tables under `reference/functions/` whose directory is "
                "named after the routine (first four; replays and captures are listed in [chapter 7](07_validation.md)).\n\n")
        o.write("Totals: %d routines, %d cited by C++ code, %d without a citation (reason in the status column).\n\n"
                % (nport + nnot, nport, nnot))
        o.write("| Fortran file | routine | lines | C++ (citing files) | status | L1 tables |\n|---|---|---|---|---|---|\n")
        o.write("\n".join(rows) + "\n")
    print("wrote %s: %d routines, %d ported, %d not" % (a.out, nport + nnot, nport, nnot))


if __name__ == "__main__":
    sys.exit(main())
