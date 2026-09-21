#!/usr/bin/env python3
"""Build a single-file, print-ready HTML reading copy of PAPER_MANUSCRIPT_V2.md.

Output: research/PAPER_MANUSCRIPT_V2.html (figures embedded as base64; equations via KaTeX from a CDN,
so equations need an internet connection when the file is opened). Open in a browser; Ctrl+P -> Save as PDF.
The Markdown remains the single source of truth: re-run after edits.
"""
import base64
import html
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RES = HERE.parent
SRC = RES / "PAPER_MANUSCRIPT_V2.md"
DST = RES / "PAPER_MANUSCRIPT_V2.html"

PH = re.compile(r"\[(CITATION REQUIRED|REFERENCE NEEDED|AUTHOR TO|DATA VERIFICATION|UNRESOLVED)")


def inline(t):
    store = []

    def keep(s):
        store.append(s)
        return "\x00%d\x00" % (len(store) - 1)

    def code(m):
        b = m.group(1)
        if PH.match(b):
            return keep('<span class="ph">' + html.escape(b) + "</span>")
        return keep("<code>" + html.escape(b) + "</code>")
    t = re.sub(r"`([^`]+)`", code, t)
    t = re.sub(r"\$([^$\n]+)\$", lambda m: keep("\\(" + html.escape(m.group(1), quote=False) + "\\)"), t)
    t = re.sub(r"(https?://[^\s`<>]+[^\s`<>.,;)])", lambda m: keep('<a href="%s">%s</a>' % (m.group(1), html.escape(m.group(1)))), t)
    t = html.escape(t, quote=False)
    t = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", t)
    t = re.sub(r"(?<![\w*])\*(?!\s)(.+?)(?<!\s)\*(?![\w*])", r"<em>\1</em>", t)
    return re.sub("\x00(\\d+)\x00", lambda m: store[int(m.group(1))], t)


def split_row(line):
    line = line.strip().strip("|") if False else line.strip()
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
    right = [a.strip().endswith(":") and not a.strip().startswith(":") for a in align]
    out = ['<figure class="tbl">']
    if caption:
        out.append("<figcaption>%s</figcaption>" % inline(caption))
    out.append("<div class='scroll'><table><thead><tr>" + "".join("<th>%s</th>" % inline(h) for h in head) + "</tr></thead><tbody>")
    for r in body:
        r = (r + [""] * n)[:n]
        out.append("<tr>" + "".join('<td class="%s">%s</td>' % ("r" if right[j] else "", inline(c)) for j, c in enumerate(r)) + "</tr>")
    out.append("</tbody></table></div></figure>")
    return "".join(out)


def img_tag(path):
    f = (RES / path).resolve()
    data = base64.b64encode(f.read_bytes()).decode()
    return '<img alt="%s" src="data:image/png;base64,%s">' % (html.escape(f.name), data)


CSS = """
:root{--ink:#1f2937;--muted:#6b7280;--accent:#0072B2;--line:#e5e7eb;--bg:#ffffff;--soft:#f3f4f6}
@media (prefers-color-scheme:dark){:root:not(.light){--ink:#e5e7eb;--muted:#9ca3af;--accent:#56B4E9;--line:#374151;--bg:#111827;--soft:#1f2937}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font:17px/1.62 Georgia,'Iowan Old Style','Times New Roman',serif}
main{max-width:860px;margin:0 auto;padding:36px 20px 80px}
h1{font:700 27px/1.25 system-ui,sans-serif;margin:.2em 0 .3em}
h2{font:700 22px/1.3 system-ui,sans-serif;color:var(--accent);margin:2em 0 .4em;border-bottom:1px solid var(--line);padding-bottom:.2em}
h3{font:700 18px/1.3 system-ui,sans-serif;margin:1.6em 0 .3em}
h4{font:700 16px/1.3 system-ui,sans-serif;margin:1.3em 0 .3em}
p{margin:.65em 0}
.byline{font:15px system-ui,sans-serif;color:var(--muted);margin:0 0 1.2em}
.note{font:14px/1.5 system-ui,sans-serif;color:var(--muted);border-left:3px solid var(--line);padding:.2em 0 .2em 12px;margin:1em 0}
.abstract{background:var(--soft);border-radius:8px;padding:6px 18px;margin:1.2em 0;font-size:16px}
code{font:13.5px ui-monospace,Menlo,Consolas,monospace;background:var(--soft);padding:1px 4px;border-radius:3px;word-break:break-word}
.ph{background:#FFF3B0;color:#1f2937;font:13px system-ui,sans-serif;padding:1px 5px;border-radius:3px}
ul{padding-left:1.3em} li{margin:.25em 0}
figure{margin:1.4em 0}
figure.tbl figcaption,figure.fig figcaption{font:14px/1.5 system-ui,sans-serif;color:var(--muted);margin:.4em 0}
figure.tbl figcaption{margin-bottom:.5em} figure.tbl figcaption strong,figure.fig figcaption strong{color:var(--accent)}
.scroll{overflow-x:auto}
table{border-collapse:collapse;width:100%;font:13px/1.4 system-ui,sans-serif}
th{text-align:left;border-top:2px solid var(--ink);border-bottom:1px solid var(--ink);padding:5px 8px;vertical-align:bottom}
td{padding:5px 8px;border-bottom:1px solid var(--line);vertical-align:top}
td.r{text-align:right;font-variant-numeric:tabular-nums}
tbody tr:nth-child(even){background:var(--soft)}
table{border-bottom:2px solid var(--ink)}
img{max-width:100%;height:auto;display:block;margin:0 auto}
.ref{padding-left:2.4em;text-indent:-2.4em;margin:.35em 0;font-size:15px}
.kw{font:15px system-ui,sans-serif;color:var(--muted)}
hr{border:0;border-top:1px solid var(--line);margin:1.6em 0}
.katex-display{overflow-x:auto;overflow-y:hidden}
@media print{
 @page{size:A4;margin:18mm 16mm}
 body{font-size:10.5pt;background:#fff;color:#000} main{max-width:none;padding:0}
 h2{break-after:avoid} h3,h4{break-after:avoid}
 figure,table,tr{break-inside:avoid} figure.tbl{break-inside:auto}
 table{font-size:8pt} .scroll{overflow:visible} .abstract{background:#f6f6f6}
 a{color:inherit;text-decoration:none}
}
"""

KATEX = """
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.css">
<script defer src="https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.js"></script>
<script defer src="https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/contrib/auto-render.min.js"
 onload="renderMathInElement(document.body,{delimiters:[{left:'$$',right:'$$',display:true},{left:'\\\\[',right:'\\\\]',display:true},{left:'\\\\(',right:'\\\\)',display:false}]})"></script>
"""


def convert():
    md = SRC.read_text(encoding="utf-8").split("\n")
    out, i = [], 0
    title, front, pending_cap, pending_fig = "", [], None, None
    in_abs = False
    in_refs = False
    keywords = ""
    seen_abs = False
    while i < len(md):
        s = md[i].strip()
        if not s or s == "---":
            i += 1; continue
        if s.startswith("# ") and not title:
            title = s[2:].strip(); i += 1; continue
        m = re.match(r"^(#{2,4})\s+(.*)$", s)
        if m:
            lvl, txt = len(m.group(1)), m.group(2).strip()
            if in_abs:
                out.append("</div>")
                if keywords:
                    out.append('<p class="kw"><strong>Keywords:</strong> %s</p>' % inline(keywords))
                in_abs = False
            in_refs = txt.startswith("References")
            if txt == "Abstract":
                out.append('<h2>Abstract</h2><div class="abstract">'); in_abs = True; seen_abs = True
            else:
                out.append("<h%d>%s</h%d>" % (lvl, inline(txt), lvl))
            i += 1; continue
        if s.startswith("$$"):
            buf = [s[2:]] if s != "$$" else []
            if s.endswith("$$") and len(s) > 4:
                out.append("<p>$$%s$$</p>" % html.escape(s[2:-2].strip(), quote=False)); i += 1; continue
            i += 1
            while i < len(md) and md[i].strip() != "$$":
                buf.append(md[i]); i += 1
            i += 1
            out.append("<p>$$%s$$</p>" % html.escape("\n".join(buf), quote=False))
            continue
        if s.startswith("!["):
            pending_fig = re.match(r"!\[[^\]]*\]\(([^)]+)\)", s).group(1); i += 1; continue
        if s.startswith("|"):
            blk = []
            while i < len(md) and md[i].strip().startswith("|"):
                blk.append(md[i]); i += 1
            out.append(table(blk, pending_cap)); pending_cap = None; continue
        if re.match(r"^\*\*Table \d+\.", s) and not s.endswith("|"):
            pending_cap = s; i += 1; continue
        mf = re.match(r"^\*\*(Figure \d+\.)\*\*\s*(.*)$", s)
        if mf and pending_fig:
            out.append('<figure class="fig">%s<figcaption><strong>%s</strong> %s</figcaption></figure>'
                       % (img_tag(pending_fig), mf.group(1), inline(mf.group(2))))
            pending_fig = None; i += 1; continue
        if s.startswith("**Keywords:**"):
            keywords = s[len("**Keywords:**"):].strip(); i += 1; continue
        if s.startswith("- "):
            items = []
            while i < len(md) and md[i].strip().startswith("- "):
                items.append(md[i].strip()[2:]); i += 1
            out.append("<ul>" + "".join("<li>%s</li>" % inline(x) for x in items) + "</ul>"); continue
        mr = re.match(r"^(\d+)\.\s+(.*)$", s)
        if mr and in_refs:
            out.append('<p class="ref">[%s]&nbsp;%s</p>' % (mr.group(1), inline(mr.group(2)))); i += 1; continue
        if not seen_abs:
            front.append(s); i += 1; continue
        out.append("<p>%s</p>" % inline(s)); i += 1
    return title, front, out


def main():
    title, front, body = convert()
    note = "".join('<div class="note">%s</div>' % inline(f) for f in front)
    doc = ("<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>"
           "<title>%s</title><style>%s</style>%s</head><body><main><h1>%s</h1>"
           "<p class='byline'>Mark Kostogryzov &middot; Independent researcher &middot; m.2researcher141941241@proton.me</p>%s%s</main></body></html>"
           % ("Binomial screen with confidence routing to a Poisson-binomial SNP caller", CSS, KATEX, inline(title), note, "\n".join(body)))
    DST.write_text(doc, encoding="utf-8")
    print("wrote", DST, round(len(doc) / 1e6, 2), "MB")


if __name__ == "__main__":
    sys.exit(main())
