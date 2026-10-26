const fs = require("fs");
const D = require("docx");
const C = require("./content.js");
const { Document, Packer, Paragraph, TextRun, ImageRun, Table, TableRow, TableCell, WidthType, ShadingType,
  BorderStyle, AlignmentType, Footer, Header, PageNumber, SectionType, VerticalAlign, CommentRangeStart,
  CommentRangeEnd, CommentReference, HeadingLevel, LevelFormat } = D;

const FONT = "Times New Roman";
const PAGE_W = 11906, PAGE_H = 16838, MARGIN = 1020, GAP = 340;
const FULL = PAGE_W - 2 * MARGIN;
const COL = Math.floor((FULL - GAP) / 2);
const BODY = 20, SMALL = 16, CAP = 16, REF = 16;

// ---------- inline markup ----------
function runs(text, base = {}) {
  text = text.replace(/`/g, "");
  const out = [];
  const re = /(\*\*.+?\*\*|\*.+?\*|\^.+?\^|~.+?~|\[\[c\d+\]\]|\[\[\/c\d+\]\])/g;
  for (const part of text.split(re)) {
    if (!part) continue;
    let m;
    if ((m = part.match(/^\[\[c(\d+)\]\]$/))) out.push(new CommentRangeStart(+m[1]));
    else if ((m = part.match(/^\[\[\/c(\d+)\]\]$/))) {
      out.push(new CommentRangeEnd(+m[1]));
      out.push(new TextRun({ children: [new CommentReference(+m[1])] }));
    }
    else if (part.startsWith("**")) out.push(new TextRun({ ...base, text: part.slice(2, -2), bold: true }));
    else if (part.startsWith("*")) out.push(new TextRun({ ...base, text: part.slice(1, -1), italics: true }));
    else if (part.startsWith("^")) out.push(new TextRun({ ...base, text: part.slice(1, -1), superScript: true }));
    else if (part.startsWith("~")) out.push(new TextRun({ ...base, text: part.slice(1, -1), subScript: true }));
    else out.push(new TextRun({ ...base, text: part }));
  }
  return out;
}
const base = (size = BODY) => ({ font: FONT, size });
const P = (t, o = {}) => new Paragraph({
  children: runs(t, base(o.size || BODY)),
  alignment: o.align || AlignmentType.JUSTIFIED,
  spacing: { before: o.before ?? 0, after: o.after ?? 70, line: o.line ?? 240 },
  keepNext: o.keepNext, indent: o.indent,
});
const H1 = (t) => new Paragraph({ heading: HeadingLevel.HEADING_1, keepNext: true, spacing: { before: 130, after: 50 },
  children: [new TextRun({ text: t, font: FONT, size: 21, bold: true })] });
const H2 = (t) => new Paragraph({ heading: HeadingLevel.HEADING_2, keepNext: true, spacing: { before: 70, after: 30 },
  children: [new TextRun({ text: t, font: FONT, size: BODY, bold: true, italics: true })] });

// ---------- figures ----------
function pngSize(path) {
  const b = fs.readFileSync(path);
  return { w: b.readUInt32BE(16), h: b.readUInt32BE(20) };
}
function figure(file, widthDxa, captionText) {
  const { w, h } = pngSize(file);
  const px = Math.round(widthDxa / 15); // 1 px = 15 dxa
  const py = Math.round(px * h / w);
  return [
    new Paragraph({ alignment: AlignmentType.CENTER, spacing: { before: 60, after: 20 }, keepNext: true,
      children: [new ImageRun({ type: "png", data: fs.readFileSync(file), transformation: { width: px, height: py },
        altText: { title: "Figure", description: captionText.slice(0, 120), name: file } })] }),
    new Paragraph({ children: runs(captionText, base(CAP)), alignment: AlignmentType.JUSTIFIED, spacing: { after: 100, line: 208 } }),
  ];
}

// ---------- tables ----------
const none = { style: BorderStyle.NONE, size: 0, color: "FFFFFF" };
const rule = (sz) => ({ style: BorderStyle.SINGLE, size: sz, color: "000000" });
function table(widths, rows, opts = {}) {
  const total = widths.reduce((a, b) => a + b, 0);
  const fs_ = opts.size || SMALL - 1;
  return new Table({
    width: { size: total, type: WidthType.DXA }, columnWidths: widths,
    borders: { top: rule(8), bottom: rule(8), left: none, right: none, insideHorizontal: none, insideVertical: none },
    rows: rows.map((r, ri) => new TableRow({
      tableHeader: ri === 0, cantSplit: true,
      children: r.map((c, ci) => new TableCell({
        width: { size: widths[ci], type: WidthType.DXA }, verticalAlign: VerticalAlign.CENTER,
        margins: { top: 15, bottom: 15, left: 40, right: 40 },
        shading: ri === 0 ? { type: ShadingType.CLEAR, fill: "EEEEEE", color: "auto" } : (opts.hl && opts.hl.includes(ri) ? { type: ShadingType.CLEAR, fill: "F4F4F4", color: "auto" } : undefined),
        borders: { top: none, left: none, right: none, bottom: ri === 0 ? rule(4) : none },
        children: [new Paragraph({ spacing: { after: 0, line: 200 }, alignment: ci === 0 || opts.leftAll ? AlignmentType.LEFT : (opts.center ? AlignmentType.CENTER : AlignmentType.LEFT),
          children: runs(String(c), { font: FONT, size: fs_, bold: ri === 0 }) })],
      })),
    })),
  });
}
const caption = (t) => new Paragraph({ children: runs(t, base(CAP)), alignment: AlignmentType.JUSTIFIED, spacing: { before: 60, after: 30, line: 208 }, keepNext: true });
const afterTable = (t) => new Paragraph({ children: runs(t, base(CAP - 1)), spacing: { before: 20, after: 100, line: 190 } });

module.exports = { D, C, runs, base, P, H1, H2, figure, table, caption, afterTable, FONT, PAGE_W, PAGE_H, MARGIN, GAP, FULL, COL, BODY, SMALL, CAP, REF };
