# Editorial report: Bioinformatics Advances working draft

Deliverable: `AI_DNA_ANALYZER_BioinformaticsAdvances_working_draft.docx` (PDF preview alongside; figure/build scripts in `build/`).
Source (not modified): `research/PAPER_MANUSCRIPT_V2.md`, built as `research/latex/paper_biorxiv_final.pdf` (36 pp.) plus `Supplementary_Information.pdf`. Result CSV/JSON files were read for figure data. `AI_DNA_ANALYZER_Scientific_Report.pdf` (15 pp., older) was not used.
Nothing was re-run or re-computed from raw data. No experiment, test or pipeline was executed for this draft.

## Narrative chosen
Three-stage development with unequal weight: V1 as the frozen, honestly scored baseline (~25 % of results text); V2 as an exact re-implementation with a limited runtime claim (~15 %); V2.x as the main quantitative result (~40 %, holdout + one pre-specified run); external callers as a bounded context section (~10 %).

## 1. Condensed substantially
- V1 methods/provenance (constant selection, extraction details, backends), the 8-experiment table (kept, shortened), failure taxonomy (3 classes → 1 paragraph), router-capture analysis, safety layer (1 sentence), genotype methods A/B/C (1 sentence), bottleneck/Amdahl analysis (2 sentences).
- V2 correctness gates and profiling (1 paragraph + table). V2.x ablation of 16 A–D combinations (reduced to components in Fig. 2b).
- Limitations (27-row table → 1 paragraph). Discussion (4 paragraphs).

## 2. Retained vs moved out of main text
Retained: v14 NEGATIVE verdict; v20 DEGRADED (ΔF1 −5.64e-5, 10 discordant loci, McNemar p 0.021); 10 lost SNPs; 1.65× as the only measured routing gain; unmatched V1/V2 runtime; 0.97195→0.98902 holdout; base-quality rule ≈92 % of gain; no HG003 control; DeepVariant/GATK higher; scope limits.
Moved out (not in draft; still in the 36-page manuscript/Supplement): chr21:32–44 Mb external comparison (DeepVariant 1.6.1, Clair3, GATK 4.5); HG005 results (record-only, raw artefacts lost); genotype-method table; failure-taxonomy and rescue tables; native-backend equivalence-gate history and test counts; Table 5 truth-SNP denominators; chr21/chr20 300× V1 hap.py detail; per-cell v14/v19 results; neural-phase history.

## 3. Claims needing author verification (also as Word comments)
- Figure 2a values are transcribed from the manuscript's Table 4, not recomputed from `results/bench_v13–v20`.
- "Base-quality rule A ≈ 92 % of the A+C gain" is computed here (0.01574/0.01706); the manuscript does not state it.
- V1 18,578 s baseline is from a frozen record under other conditions; ratios (~17.5×, ~102×) indicative only.
- Statement that no V1 constant was changed for V2.x, and the disclosure about B-variant holdout numbers printed in an exploratory script.
- Test-suite counts (490 passed) are from the 2026-09-20 freeze report; not re-run.
- Two Zenodo DOIs exist (v1.0.0 …22895746, v2.0.0 …22978775); confirm which to cite.
- AI-assistance statement vs. journal policy.

## 4. Figures
Real data: Fig. 2b (`V2X_BOOTSTRAP_DELTAS_VS_V2.csv`, holdout rows); Fig. 3a (precision/recall of Table 5 callers; chr20 V2 control/A+C from `V2X_ACCURACY_PERFORMANCE_BENCHMARK.csv`, DeepVariant/GATK values from the manuscript's Table 16); Fig. 3b (FP/FN counts from the same CSV and `final_independent_eval/results_chr8_analysis.json`; the 79 % outside-window figure from the manuscript); Fig. 2a (transcribed from manuscript Table 4, see above).
Schematic: Fig. 1 (no data).
No placeholders were needed. Not plotted (no per-cell values in final tables): v15–v18 cells.

## 5. Reference verification
All 22 DOIs in the source bibliography were resolved against Crossref on 2026-10-25 (title/first author/journal/volume/pages matched). 18 references are in the draft. Not machine-verified: arXiv items (Mamba, FreeBayes; both dropped from the draft). Flags: Poplin 2017 is a bioRxiv preprint; Zheng 2022 is the long-read Clair3 paper (Clair3 itself is not discussed in the draft except through this citation of tiered design); miniSNV is a nanopore caller, cited only as an example of tiered processing. Amdahl, Efron/Efron–Tibshirani, Strelka2, Mamba, FreeBayes, Dwarshuis, HTSlib-adjacent software refs were dropped with the text that used them.

## 6. Author decisions still needed
1. Affiliation, ORCID, e-mail, funding (placeholders).
2. Confirm journal article type, page/word limit and template; replace layout with the OUP template.
3. Abstract format (structured, 184 words now) against current guidelines.
4. Whether V1 detail (decision-rule defect, margin origin) belongs in main text or Supplement.
5. Which Zenodo DOI(s) to cite; whether to add a V2 control on HG003 chr8 (would turn the final run into a measured gain).
6. Whether to keep the "no claim of novelty" sentence (source says only a brief literature search was done).
7. Whether to include the chr21 external comparison (incl. Clair3) in a Supplement.

## 7. Page count and formatting
Rendered with LibreOffice: **7 pages** (A4, two columns, 10-pt Times body, 8-pt tables/captions, 8-pt references, 3 figures, 5 tables). Word may differ by roughly ±0.5 page. This is an assumed layout, not the OUP template; the OUP guidelines page returned HTTP 403, so journal limits, structured-abstract format and author-year referencing are assumptions. 17+2 Word comments mark verification points. The DOCX passed the skill's schema validation; it was not opened in Microsoft Word.
