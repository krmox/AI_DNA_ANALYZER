#!/usr/bin/env python3
"""Convert research/PAPER_MANUSCRIPT_V2.md to research/latex/paper.tex.

No content is edited: the Markdown is the single source of truth. Re-run after any
change to the manuscript:  python3 md2tex.py
Engine: pdflatex (all Unicode is mapped to LaTeX commands; see UNI below).
"""
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE.parent / "PAPER_MANUSCRIPT_V2.md"
DST = HERE / "paper.tex"

# (text-mode, math-mode) replacements for non-ASCII characters
UNI = {
    "–": ("--", "-"), "—": ("---", "-"), "−": ("$-$", "-"), "×": ("$\\times$", "\\times "),
    "≈": ("$\\approx$", "\\approx "), "≤": ("$\\le$", "\\le "), "≥": ("$\\ge$", "\\ge "),
    "Δ": ("$\\Delta$", "\\Delta "), "ε": ("$\\varepsilon$", "\\varepsilon "), "δ": ("$\\delta$", "\\delta "),
    "τ": ("$\\tau$", "\\tau "), "θ": ("$\\theta$", "\\theta "), "±": ("$\\pm$", "\\pm "),
    "→": ("$\\to$", "\\to "), "‡": ("\\ddag{}", "\\ddag "), "†": ("\\dag{}", "\\dag "),
    "…": ("\\ldots{}", "\\ldots "), "’": ("'", "'"), "‘": ("`", "`"), "“": ("``", "``"), "”": ("''", "''"),
    "·": ("$\\cdot$", "\\cdot "), "≠": ("$\\ne$", "\\ne "), "≫": ("$\\gg$", "\\gg "), "≪": ("$\\ll$", "\\ll "),
    "‑": ("-", "-"), "\u00a0": ("~", "~"), "∞": ("$\\infty$", "\\infty "), "°": ("$^{\\circ}$", "^{\\circ}"),
    "μ": ("$\\mu$", "\\mu "), "σ": ("$\\sigma$", "\\sigma "), "é": ("\\'e", "e"), "ö": ('\\"o', "o"),
    "ü": ('\\"u', "u"), "ä": ('\\"a', "a"), "½": ("$\\frac12$", "\\frac12 "), "✓": ("", ""),
    "§": ("\\S{}", "\\S "), "⁻": ("", ""), "₀": ("$_0$", "_0"), "₁": ("$_1$", "_1"),
}
SUP = {"⁻": "-", "⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4", "⁵": "5", "⁶": "6", "⁷": "7", "⁸": "8", "⁹": "9"}
ESC = {"\\": "\\textbackslash{}", "&": "\\&", "%": "\\%", "$": "\\$", "#": "\\#", "_": "\\_", "{": "\\{", "}": "\\}",
       "~": "\\textasciitilde{}", "^": "\\textasciicircum{}"}

PH = re.compile(r"\[(CITATION REQUIRED|REFERENCE NEEDED|AUTHOR TO|DATA VERIFICATION|UNRESOLVED)")


def sup_runs(t):
    return re.sub("[" + "".join(SUP) + "]+", lambda m: "$^{" + "".join(SUP[c] for c in m.group()) + "}$", t)


def uni(t, math=False):
    if math:
        t = re.sub("[" + "".join(SUP) + "]+", lambda m: "^{" + "".join(SUP[c] for c in m.group()) + "}", t)
    else:
        t = sup_runs(t)
    out = []
    for c in t:
        if ord(c) < 128:
            out.append(c)
        elif c in UNI:
            out.append(UNI[c][1 if math else 0])
        else:
            out.append("<<UNMAPPED U+%04X>>" % ord(c))
    return "".join(out)


def esc_text(t):
    t = uni(t)
    return "".join(ESC.get(c, c) for c in t) if "<<" not in t else "".join(ESC.get(c, c) for c in t)


def esc_text_keep_math(t):
    """Escape plain text but keep the $...$ wrappers that uni() inserted (they are added after escaping)."""
    parts = []
    for c in t:
        parts.append(ESC.get(c, c) if ord(c) < 128 else c)
    return uni("".join(parts))


def allow_breaks(s):
    s = s.replace("/", "/\\allowbreak{}").replace("\\_", "\\_\\allowbreak{}").replace(".", ".\\allowbreak{}")
    # long hex runs (SHA-256 digests) have no break points and overrun the margin
    return re.sub(r"[0-9a-fA-F]{16,}", lambda m: "\\allowbreak{}".join(m.group(0)[i:i + 8] for i in range(0, len(m.group(0)), 8)), s)


def inline(t):
    """Markdown inline -> LaTeX. Order: protect code/math/urls, escape, then emphasis, then restore."""
    store = []

    def keep(s):
        store.append(s)
        return "\x00%d\x00" % (len(store) - 1)

    def code(m):
        body = m.group(1)
        if PH.match(body):
            # one box per word so a long placeholder can wrap instead of overrunning the margin
            return keep(" ".join("\\placeholder{" + esc_text_keep_math(w) + "}" for w in body.split(" ")))
        return keep("\\texttt{" + allow_breaks(esc_text_keep_math(body)) + "}")

    t = re.sub(r"`([^`]+)`", code, t)
    t = re.sub(r"\$([^$\n]+)\$", lambda m: keep("$" + uni(m.group(1), True) + "$"), t)

    def url(m):
        u = m.group(0).rstrip(".,;)")
        tail = m.group(0)[len(u):]
        return keep("\\url{" + u + "}") + tail
    t = re.sub(r"https?://[^\s`]+", url, t)
    t = re.sub(r"!\[[^\]]*\]\([^)]*\)", "", t)
    # straight ASCII quotes render as vertical ticks in print; pair them into typographic quotes
    # (code/math/urls are already protected above, so this only touches prose)
    if t.count('"') % 2 == 0:
        parts = t.split('"')
        t = parts[0] + "".join(("``" if i % 2 else "''") + p for i, p in enumerate(parts[1:], 1))
    t = esc_text_keep_math(t)
    t = re.sub(r"\*\*(.+?)\*\*", r"\\textbf{\1}", t)
    t = re.sub(r"(?<![\w*])\*(?!\s)(.+?)(?<!\s)\*(?![\w*])", r"\\emph{\1}", t)
    # doi/identifier strings: allow line breaks
    t = re.sub(r"(doi:\S+)", lambda m: allow_breaks(m.group(1)), t)
    t = re.sub(r"(\d) \\%", r"\1\\,\\%", t)
    return re.sub("\x00(\\d+)\x00", lambda m: store[int(m.group(1))], t)


def split_row(line):
    line = line.strip()
    if line.startswith("|"):
        line = line[1:]
    if line.endswith("|"):
        line = line[:-1]
    cells, cur, in_m, in_c = [], [], False, False
    for ch in line:
        if ch == "$" and not in_c:
            in_m = not in_m
        elif ch == "`":
            in_c = not in_c
        if ch == "|" and not in_m and not in_c:
            cells.append("".join(cur).strip()); cur = []
        else:
            cur.append(ch)
    cells.append("".join(cur).strip())
    return cells


def table(lines, caption):
    rows = [split_row(l) for l in lines]
    head, align, body = rows[0], rows[1], rows[2:]
    n = len(head)
    body = [r + [""] * (n - len(r)) if len(r) < n else r[:n] for r in body]
    # Columns are either fixed-width (short content) or flexible X columns. X columns get widths
    # proportional to their content length (tabularx \hsize weights that sum to the number of X columns),
    # so a wide cell such as "+0.000103 [+0.000025, +0.000180]" is not squeezed by equal-share columns.
    lens = [max([len(head[j])] + [len(r[j]) for r in body]) for j in range(n)]
    flex = [j for j in range(n) if lens[j] > 12]
    if not flex:
        flex = [max(range(n), key=lambda j: lens[j])]
    def longest_word(j):
        return max(len(w) for c in [head[j]] + [r[j] for r in body] for w in (c.split() or [""]))
    wt = {j: min(max(lens[j], 12, 2.2 * longest_word(j)), 60) for j in flex}
    scale = len(flex) / float(sum(wt.values()))
    fw = {j: wt[j] * scale for j in flex}
    # A column scaled well below the group average can end up narrower than its own header
    # word (e.g. "Evaluator"), which then overflows/clips instead of wrapping. Floor every
    # flex column at 65% of the average share and take the difference from the columns that
    # have room to spare, proportional to their surplus over that floor.
    min_share = 0.65
    deficit = sum(max(0.0, min_share - fw[j]) for j in flex)
    if deficit > 0:
        donors = [j for j in flex if fw[j] > min_share]
        surplus_total = sum(fw[j] - min_share for j in donors)
        for j in flex:
            if fw[j] < min_share:
                fw[j] = min_share
            elif surplus_total > 0:
                fw[j] -= deficit * (fw[j] - min_share) / surplus_total
    spec = []
    for j in range(n):
        a = align[j].strip()
        if j in flex:
            spec.append(">{\\hsize=%.3f\\hsize\\raggedright\\arraybackslash}X" % fw[j])
        elif a.endswith(":") and not a.startswith(":"):
            spec.append("r")
        else:
            spec.append("l")
    size = "\\scriptsize" if n >= 8 else "\\footnotesize"
    out = ["{" + size, "\\setlength{\\tabcolsep}{3.2pt}", "\\rowcolors{2}{rowgrey}{white}",
           "\\begin{xltabular}{\\linewidth}{" + "".join(spec) + "}"]
    hdr = " & ".join("\\textbf{" + inline(h) + "}" if h else "" for h in head) + " \\\\"
    if caption:
        cap = inline(caption)
        out += ["\\caption{" + cap + "}\\\\", "\\toprule", hdr, "\\midrule", "\\endfirsthead",
                "\\caption[]{\\textit{(continued)}}\\\\", "\\toprule", hdr, "\\midrule", "\\endhead"]
    else:
        out += ["\\toprule", hdr, "\\midrule", "\\endhead"]
    out += ["\\bottomrule", "\\endlastfoot"] if False else []
    for r in body:
        out.append(" & ".join(inline(c) for c in r) + " \\\\")
    out += ["\\bottomrule", "\\end{xltabular}", "}", ""]
    return "\n".join(out)


def convert():
    md = SRC.read_text(encoding="utf-8").split("\n")
    out, i = [], 0
    title = ""
    abstract_open = False
    fig_counter_set = False
    pending_caption = None   # table caption
    pending_fig = None       # figure path
    seen_abstract = False
    in_refs = False
    last_table_no, last_fig_no = 0, 0
    front = []  # paragraphs between title and Abstract
    keywords = ""
    while i < len(md):
        ln = md[i]
        s = ln.strip()
        if not s or s == "---":
            i += 1; continue
        if s.startswith("# ") and not title:
            title = s[2:].strip(); i += 1; continue
        m = re.match(r"^(#{2,4})\s+(.*)$", s)
        if m:
            lvl, txt = len(m.group(1)), m.group(2).strip()
            num = re.match(r"^\d+(\.\d+)*\.?\s+", txt)
            body = txt[num.end():] if num else txt
            in_refs = body.startswith("References")
            if body == "Abstract":
                out.append("\\begin{abstract}"); seen_abstract = True; abstract_open = True
                i += 1; continue
            cmd = {2: "section", 3: "subsection", 4: "subsubsection"}[lvl]
            star = "" if num else "*"
            if abstract_open:
                out.append("\\end{abstract}")
                if keywords:
                    out.append("\\noindent\\textbf{Keywords:} " + keywords + "\\par\\bigskip")
                abstract_open = False
            out.append("\\%s%s{%s}" % (cmd, star, inline(body)))
            if star and cmd == "section":
                out.append("\\addcontentsline{toc}{section}{%s}" % inline(body))
            i += 1; continue
        if s.startswith("$$") and s.endswith("$$") and len(s) > 4:
            out.append("\\[\n" + uni(s[2:-2].strip(), True) + "\n\\]")
            i += 1; continue
        if s.startswith("$$"):
            buf = []
            if s != "$$":
                buf.append(s[2:])
            i += 1
            while i < len(md) and md[i].strip() != "$$":
                buf.append(md[i]); i += 1
            i += 1
            out.append("\\[\n" + uni("\n".join(buf), True) + "\n\\]")
            continue
        if s.startswith("!["):
            mm = re.match(r"!\[[^\]]*\]\(([^)]+)\)", s)
            pending_fig = "../" + mm.group(1)
            i += 1; continue
        if s.startswith("|"):
            blk = []
            while i < len(md) and md[i].strip().startswith("|"):
                blk.append(md[i]); i += 1
            out.append(table(blk, pending_caption)); pending_caption = None
            continue
        mt = re.match(r"^\*\*Table (\d+)\.\s*(.*)$", s)
        if mt and not s.endswith("|"):
            no = int(mt.group(1))
            assert no == last_table_no + 1, "table numbering gap: %d after %d" % (no, last_table_no)
            last_table_no = no
            rest = mt.group(2)
            # caption = bold title + any trailing text
            rest = re.sub(r"\*\*", "", rest, count=1) if "**" in rest else rest
            rest = rest.replace("**", "")
            pending_caption = rest
            i += 1; continue
        mf = re.match(r"^\*\*Figure (\d+)\.\*\*\s*(.*)$", s)
        if mf and pending_fig:
            no = int(mf.group(1)); assert no > last_fig_no
            last_fig_no = no
            out.append("\\setcounter{figure}{%d}" % (no - 1))
            out.append("\\begin{figure}[tbp]\n\\centering\n\\includegraphics[width=\\linewidth]{%s}\n\\caption{%s}\n\\end{figure}\n\\FloatBarrier"
                       % (pending_fig, inline(mf.group(2))))
            pending_fig = None; i += 1; continue
        if s.startswith("**Keywords:**"):
            keywords = inline(s[len("**Keywords:**"):].strip()); i += 1; continue
        if re.match(r"^- ", s):
            items = []
            while i < len(md) and md[i].strip().startswith("- "):
                items.append(md[i].strip()[2:]); i += 1
            out.append("\\begin{itemize}[leftmargin=1.6em,itemsep=2pt,topsep=3pt]\n" +
                       "\n".join("\\item " + inline(x) for x in items) + "\n\\end{itemize}")
            continue
        mr = re.match(r"^(\d+)\.\s+(.*)$", s)
        if mr and in_refs:
            out.append("\\refitem{%s}{%s}" % (mr.group(1), inline(mr.group(2)))); i += 1; continue
        if not seen_abstract:
            front.append(s); i += 1; continue
        out.append(inline(s) + "\n"); i += 1
    return title, front, out


def main():
    title, front, body = convert()
    pre = (HERE / "preamble.tex").read_text(encoding="utf-8")
    note = ""
    if front:
        note = "\\begin{quote}\\small " + " ".join(inline(f) for f in front) + "\\end{quote}"
    doc = pre.replace("%%TITLE%%", inline(title)).replace("%%NOTE%%", note)
    doc = doc.replace("%%BODY%%", "\n".join(body))
    DST.write_text(doc, encoding="utf-8")
    bad = re.findall(r"<<UNMAPPED U\+[0-9A-F]+>>", doc)
    print("wrote", DST, "| unmapped:", sorted(set(bad)))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
