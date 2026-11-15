#!/usr/bin/env python3
"""check_doc_coverage.py - Doxygen documentation coverage of the public headers, without Doxygen (task P15.D3).

usage: check_doc_coverage.py INCLUDE_DIR [--list] [--quiet]

Scans every *.h under INCLUDE_DIR with a small C++ declaration scanner and reports each entity that Doxygen (with
WARN_IF_UNDOCUMENTED=YES, WARN_NO_PARAMDOC=YES, EXTRACT_PRIVATE=NO, INHERIT_DOCS=YES) would report:

  * a file without a `@file` block;
  * a class / struct / union / enum / enumerator / type alias / function / variable / macro at namespace, class
    (public or protected section) or enum scope without a doc comment: a preceding `///`, `//!`, `/** */`, `/*! */`
    block or a trailing `///<` / `//!<` / `/**< */` comment;
  * a documented function whose named parameters lack an `@param` (or `\\param`), or whose non-void return value lacks
    an `@return` / `@returns` / `@retval` (constructors, destructors and conversion operators excepted);
    `@copydoc` satisfies both.

Exceptions that follow Doxygen: private members (EXTRACT_PRIVATE=NO); undocumented `override`/`final` functions
(they inherit the documentation of the base class, INHERIT_DOCS=YES); include guards; forward declarations; inheriting
constructors (`using Base::Base;`); anything inside an anonymous namespace or a function body. A comment block that only
holds grouping commands (`@addtogroup`, `@{`, `@}`, `@name`) documents nothing.

The scanner is heuristic but conservative: it understands the subset of C++ used in include/simc. When Doxygen is
available, the CTest docs_doxygen runs it as well (docs/doxygen/run_doxygen.cmake).

Exit status: 0 if every entity is documented, 1 otherwise (each problem printed as file:line: message), 2 usage.
"""
import os
import re
import sys

DOC_START = ("///", "//!", "/**", "/*!")
GROUP_ONLY = re.compile(r"^\s*(?:[@\\](?:addtogroup\s+\w+(?:\s+.*)?|\{|\}|name\b.*)\s*)*$", re.S)
SPECIFIERS = {"static", "inline", "virtual", "explicit", "constexpr", "friend", "extern", "mutable", "typename",
              "[[nodiscard]]", "consteval", "constinit"}


class Tok:
    """A lexical piece: kind in {'code', 'doc', 'tdoc' (trailing doc), 'cmt' (plain comment), 'pp' (preprocessor)}."""

    def __init__(self, kind, text, line):
        self.kind, self.text, self.line = kind, text, line


def lex(src):
    """Split the source into code characters, comments and preprocessor lines (string literals blanked)."""
    toks, i, n, line = [], 0, len(src), 1
    code, code_line = [], 1
    at_line_start = True

    def flush():
        if code:
            toks.append(Tok("code", "".join(code), code_line))
            code.clear()

    while i < n:
        c = src[i]
        if at_line_start and c == "#":
            j = i
            while True:
                k = src.find("\n", j)
                if k < 0:
                    k = n
                if k > 0 and src[k - 1] == "\\":
                    j = k + 1
                    continue
                break
            flush()
            text = src[i:k]
            toks.append(Tok("pp", text, line))
            line += text.count("\n")
            i = k
            continue
        if c == "\n":
            line += 1
            at_line_start = True
            if not code:
                code_line = line
            code.append(c)
            i += 1
            continue
        if c in " \t":
            if not code:
                code_line = line
            code.append(c)
            i += 1
            continue
        at_line_start = False
        if src.startswith("//", i):
            k = src.find("\n", i)
            k = n if k < 0 else k
            text = src[i:k]
            flush()
            kind = "cmt"
            if text.startswith(("///<", "//!<")):
                kind = "tdoc"
            elif text.startswith(("///", "//!")) and not text.startswith("////"):
                kind = "doc"
            toks.append(Tok(kind, text, line))
            i = k
            continue
        if src.startswith("/*", i):
            k = src.find("*/", i + 2)
            k = n if k < 0 else k + 2
            text = src[i:k]
            flush()
            kind = "cmt"
            if text.startswith(("/**<", "/*!<")):
                kind = "tdoc"
            elif (text.startswith("/**") and not text.startswith("/**/")) or text.startswith("/*!"):
                kind = "doc"
            toks.append(Tok(kind, text, line))
            line += text.count("\n")
            i = k
            continue
        if c in "\"'":
            # string / char literal (raw strings are not used in the headers)
            q, j = c, i + 1
            while j < n and src[j] != q:
                j += 2 if src[j] == "\\" else 1
            if not code:
                code_line = line
            code.append(q + q)
            i = j + 1
            continue
        if not code:
            code_line = line
        code.append(c)
        i += 1
    flush()
    return toks


def chars(toks):
    """Yield (char, line, attached) where attached is a list of comment tokens preceding the char."""
    for t in toks:
        if t.kind == "code":
            line = t.line
            for ch in t.text:
                yield ("c", ch, line)
                if ch == "\n":
                    line += 1
        else:
            yield (t.kind, t, t.line)


def strip_template(s):
    s = s.strip()
    while s.startswith("template"):
        m = re.match(r"template\s*<", s)
        if not m:
            break
        depth, j = 1, m.end()
        while j < len(s) and depth:
            depth += {"<": 1, ">": -1}.get(s[j], 0)
            j += 1
        s = s[j:].strip()
    return s


def split_top(s, sep=","):
    out, depth, cur = [], 0, []
    for ch in s:
        if ch in "([{<":
            depth += 1
        elif ch in ")]}>":
            depth -= 1
        if ch == sep and depth == 0:
            out.append("".join(cur))
            cur = []
        else:
            cur.append(ch)
    out.append("".join(cur))
    return out


def param_names(plist):
    names = []
    plist = plist.strip()
    if plist in ("", "void"):
        return names
    for p in split_top(plist):
        p = split_top(p, "=")[0].strip()
        if not p or p == "...":
            continue
        if re.search(r"\(\s*\*\s*(\w+)\s*\)", p):  # function pointer parameter
            names.append(re.search(r"\(\s*\*\s*(\w+)\s*\)", p).group(1))
            continue
        p = re.sub(r"\[[^\]]*\]$", "", p).strip()
        toks = re.findall(r"[A-Za-z_]\w*(?:::\w+)*|<[^<>]*(?:<[^<>]*>[^<>]*)*>|[&*]+", p)
        toks = [t for t in toks if t not in ("const", "volatile", "struct", "class", "enum", "typename")]
        if len(toks) >= 2 and re.match(r"^[A-Za-z_]\w*$", toks[-1]):
            names.append(toks[-1])
    return names


class Scope:
    def __init__(self, kind, name="", access="public"):
        self.kind, self.name, self.access = kind, name, access


class Checker:
    def __init__(self, path, rel):
        self.path, self.rel = path, rel
        self.problems = []
        self.entities = 0

    def problem(self, line, msg):
        self.problems.append("%s:%d: %s" % (self.rel, line, msg))

    def run(self):
        src = open(self.path, encoding="utf-8").read()
        toks = lex(src)
        if not any(t.kind == "doc" and re.search(r"[@\\]file\b", t.text) for t in toks):
            self.problem(1, "no @file block")
        self.scan(toks)

    # ------------------------------------------------------------------------------------------------------------
    def scan(self, toks):
        stack = [Scope("file")]
        stmt, stmt_line = [], 0
        doc = []            # doc blocks collected before the current statement
        last = None         # last completed entity (receives trailing docs)
        cur = None          # entity dict of the statement in progress (receives in-statement trailing docs)
        pending_class = None
        prev_pp = None
        guards = set()
        it = list(chars(toks))
        i = 0
        anon_ns_depth = 0

        def new_entity(kind, name, line, extra=None):
            e = dict(kind=kind, name=name, line=line, docs=list(doc), extra=extra or {})
            return e

        def finish(e):
            nonlocal last
            if e is None:
                return
            self.check(e, stack)
            last = e

        while i < len(it):
            kind, val, line = it[i]
            i += 1
            top = stack[-1]
            if kind == "pp":
                text = val.text
                m = re.match(r"#\s*ifndef\s+(\w+)", text)
                if m:
                    guards.add(m.group(1))
                m = re.match(r"#\s*define\s+(\w+)(\(([^)]*)\))?", text)
                if m and m.group(1) not in guards:
                    e = dict(kind="macro", name=m.group(1), line=val.line, docs=list(doc),
                             extra=dict(params=[p.strip() for p in (m.group(3) or "").split(",") if p.strip()]
                                        if m.group(2) else [], ret=False))
                    if stack[-1].kind != "body":
                        finish(e)
                doc = []
                continue
            if kind == "doc":
                if not "".join(stmt).strip():
                    doc.append(val.text)
                    continue
                # a doc comment inside a statement: e.g. between two enumerators; treat as leading doc of the rest
                if top.kind == "enum":
                    doc.append(val.text)
                continue
            if kind == "tdoc":
                if top.kind == "enum" and cur is not None:
                    cur["docs"].append(val.text)
                elif "".join(stmt).strip() and cur is not None:
                    cur["docs"].append(val.text)
                elif last is not None:
                    last["docs"].append(val.text)
                    # re-check the last entity now that it has a trailing doc (enumerators are checked at the end)
                    if not last.get("deferred"):
                        self.recheck(last)
                continue
            if kind == "cmt":
                continue
            ch = val
            if top.kind == "body":
                if ch == "{":
                    stack.append(Scope("body"))
                elif ch == "}":
                    stack.pop()
                    if stack[-1].kind != "body" and top.name == "fn":
                        # end of an inline function definition
                        stmt, doc, cur = [], [], None
                    elif stack[-1].kind != "body" and top.name == "init":
                        stmt.append("{}")
                continue
            if top.kind == "enum":
                if ch in ",}":
                    text = "".join(stmt).strip()
                    if text:
                        name = re.match(r"(\w+)", text)
                        if cur is None:
                            cur = new_entity("enumerator", "", stmt_line)
                        if not cur["name"]:
                            cur["name"] = name.group(1) if name else text
                        cur["deferred"] = True
                        if self.pending_enum_visible:
                            self.pending.append(cur)
                        last = cur
                    stmt, doc, cur = [], [], None
                    if ch == "}":
                        stack.pop()
                        stmt = ["@enum_done"]
                    continue
                if not "".join(stmt).strip() and not ch.isspace():
                    stmt_line = line
                    cur = new_entity("enumerator", "", line)
                if cur is not None and not cur["name"] and re.match(r"\w", ch) and not "".join(stmt).strip():
                    pass
                stmt.append(ch)
                if cur is not None and not cur["name"]:
                    m = re.match(r"\s*(\w+)", "".join(stmt))
                    if m and not re.match(r"\w", ch):
                        cur["name"] = m.group(1)
                continue
            # namespace / class scope
            if ch.isspace() and not "".join(stmt).strip():
                continue
            if not "".join(stmt).strip():
                stmt_line = line
            if ch == ":" and top.kind == "class":
                text = "".join(stmt).strip()
                if text in ("public", "private", "protected"):
                    top.access = text
                    stmt, doc = [], []
                    continue
            if ch == "{":
                text = strip_template("".join(stmt))
                if re.match(r"namespace\b", text) or re.match(r'extern\s+""', text):
                    anon = re.fullmatch(r"namespace\s*", text) is not None
                    stack.append(Scope("ns" if not anon else "anon"))
                    stmt, doc, cur = [], [], None
                    continue
                m = re.match(r"(?:(?:struct|class|union)\s+(?:\[\[[^\]]*\]\]\s*)?(\w+)?)\s*(?:final\s*)?(?::[^{]*)?$",
                             text)
                if m and re.match(r"(struct|class|union)\b", text) and "(" not in text:
                    ck = text.split()[0]
                    e = new_entity(ck, m.group(1) or "", stmt_line)
                    if m.group(1):
                        finish(e) if self.visible(stack) else None
                    stack.append(Scope("class", m.group(1) or "", "private" if ck == "class" else "public"))
                    stack[-1].entity = e
                    stmt, doc, cur = [], [], None
                    continue
                m = re.match(r"enum\s+(?:class\s+|struct\s+)?(\w+)", text)
                if m:
                    e = new_entity("enum", m.group(1), stmt_line)
                    if self.visible(stack):
                        finish(e)
                    stack.append(Scope("enum", m.group(1)))
                    self.pending_enum_visible = self.visible(stack[:-1])
                    stmt, doc, cur = [], [], None
                    continue
                # function body or brace initializer
                if self.is_function(text) and text.count("(") == text.count(")"):
                    e = self.function_entity(text, stmt_line, doc, stack)
                    if e is not None and self.visible(stack):
                        finish(e)
                    stack.append(Scope("body", "fn"))
                    continue
                stack.append(Scope("body", "init"))
                continue
            if ch == "}":
                # end of class / namespace
                if top.kind == "class":
                    stack.pop()
                    ent = getattr(top, "entity", None)
                    stmt = ["@class_done"]
                    doc = []  # a variable declared with the class (`struct {...} x;`) needs its own trailing ///<
                    cur = None
                    continue
                stack.pop()
                stmt, doc, cur = [], [], None
                continue
            if ch == ";":
                text = strip_template("".join(stmt)).strip()
                text = text.replace("@class_done", "").replace("@enum_done", "").strip()
                stmt_done_line = stmt_line
                if text and self.visible(stack) and stack[-1].kind != "anon":
                    e = self.declaration(text, stmt_done_line, doc, stack)
                    if e is not None:
                        if cur is not None:
                            e["docs"] += cur["docs"]
                        finish(e)
                stmt, doc, cur = [], [], None
                continue
            stmt.append(ch)
            if cur is None:
                cur = dict(kind="?", name="", line=line, docs=[], extra={})

    # ------------------------------------------------------------------------------------------------------------
    pending = None

    def visible(self, stack):
        for s in stack:
            if s.kind == "anon" or s.kind == "body":
                return False
            if s.kind == "class" and s.access == "private":
                return False
        return True

    @staticmethod
    def is_function(text):
        t = re.sub(r"<[^<>]*(?:<[^<>]*>[^<>]*)*>", "<>", text)
        if "(" not in t:
            return False
        pre = t.split("(")[0]
        if "=" in pre and "operator" not in pre:
            return False
        return re.search(r"(?:\w|operator\s*\S+|~\w+)\s*$", pre.strip()) is not None

    def function_entity(self, text, line, doc, stack):
        t = re.sub(r"\[\[[^\]]*\]\]", "", text).strip()
        # name and parameter list
        depth, start = 0, None
        for j, ch in enumerate(t):
            if ch == "(":
                if depth == 0 and start is None:
                    start = j
                depth += 1
            elif ch == ")":
                depth -= 1
                if depth == 0 and start is not None:
                    end = j
                    break
        else:
            return None
        pre, plist, post = t[:start].strip(), t[start + 1:end], t[end + 1:]
        m = re.search(r"(operator\s*(?:\(\)|[^\s(]+)|~?\w+)$", pre)
        name = m.group(1) if m else pre
        rettype = pre[:m.start()].strip() if m else ""
        words = [w for w in rettype.split() if w not in SPECIFIERS]
        rettype = " ".join(words)
        cls = stack[-1].name if stack[-1].kind == "class" else ""
        base = name.split("::")[-1]
        ctor = (base == cls) or base.startswith("~") or name.startswith("operator ") and not rettype or not rettype
        override = re.search(r"\b(override|final)\b", post) is not None
        deleted = re.search(r"=\s*delete", post) is not None
        ret = not ctor and rettype not in ("void",) and not deleted
        return dict(kind="function", name=name, line=line, docs=list(doc),
                    extra=dict(params=param_names(plist), ret=ret, override=override))

    def declaration(self, text, line, doc, stack):
        if re.match(r"(class|struct|union|enum)\s+(class\s+)?\w+\s*$", text):
            return None  # forward declaration
        if re.match(r"(friend|static_assert|using\s+namespace)\b", text):
            return None
        m = re.match(r"using\s+(\w+)\s*=", text)
        if m:
            return dict(kind="alias", name=m.group(1), line=line, docs=list(doc), extra={})
        if re.match(r"using\s+[\w:]+::\w+\s*$", text):
            return None  # using-declaration (inheriting constructors)
        m = re.match(r"typedef\b.*?(\w+)\s*$", text)
        if m:
            return dict(kind="alias", name=m.group(1), line=line, docs=list(doc), extra={})
        if self.is_function(text) and not re.search(r"\(\s*\*\s*\w+\s*\)\s*\(", text):
            return self.function_entity(text, line, doc, stack)
        # variable(s)
        t = re.sub(r"\{\}", "", text)
        t = split_top(t, "=")[0]
        m = re.search(r"\(\s*\*\s*(\w+)\s*\)", t)
        name = m.group(1) if m else (re.findall(r"(\w+)\s*(?:\[[^\]]*\])?\s*$", split_top(t)[0]) or [t])[0]
        return dict(kind="variable", name=name, line=line, docs=list(doc), extra={})

    # ------------------------------------------------------------------------------------------------------------
    def check(self, e, stack=None):
        self.entities += 1
        e["checked"] = True
        self.evaluate(e, record=True)

    def recheck(self, e):
        # a trailing doc arrived after the entity was checked: drop its earlier problems and evaluate again
        key = "%s:%d: " % (self.rel, e["line"])
        self.problems = [p for p in self.problems if not (p.startswith(key) and ("'" + e["name"] + "'") in p)]
        self.evaluate(e, record=True)

    def evaluate(self, e, record):
        docs = [d for d in e["docs"] if not GROUP_ONLY.match(re.sub(r"^(/\*[*!]<?|//[/!]<?)|\*/$", "", d,
                                                                       flags=re.M).replace("*", " ").strip())]
        if e["kind"] == "function" and e["extra"].get("override") and not docs:  # INHERIT_DOCS
            return
        if not docs:
            self.problem(e["line"], "%s '%s' is not documented" % (e["kind"], e["name"]))
            return
        if e["kind"] in ("function", "macro"):
            text = "\n".join(docs)
            if re.search(r"[@\\]copydoc\b", text):
                return
            for p in e["extra"].get("params", []):
                if not re.search(r"[@\\]param(?:\[[^\]]*\])?\s+(?:\w+\s*,\s*)*" + re.escape(p) + r"\b", text):
                    self.problem(e["line"], "%s '%s': parameter '%s' not documented (@param)" %
                                 (e["kind"], e["name"], p))
            if e["extra"].get("ret") and not re.search(r"[@\\](return|returns|retval)\b", text):
                self.problem(e["line"], "%s '%s': return value not documented (@return)" % (e["kind"], e["name"]))


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if len(args) != 1 or not os.path.isdir(args[0]):
        print(__doc__)
        return 2
    inc = args[0]
    total, problems, nfiles = 0, [], 0
    for d, _, fs in sorted(os.walk(inc)):
        for f in sorted(fs):
            if not f.endswith(".h"):
                continue
            p = os.path.join(d, f)
            c = Checker(p, os.path.relpath(p, os.path.dirname(os.path.abspath(inc))))
            c.pending = []
            c.run()
            # enumerators collected during the scan
            for e in c.pending:
                c.check(e)
            nfiles += 1
            total += c.entities
            problems += c.problems
            if "--list" in sys.argv:
                print("%-60s %4d entities, %3d problems" % (c.rel, c.entities, len(c.problems)))
    if "--quiet" not in sys.argv:
        for p in problems:
            print(p)
    print("doc coverage: %d headers, %d entities, %d problems -> %s" %
          (nfiles, total, len(problems), "PASS" if not problems else "FAIL"))
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
