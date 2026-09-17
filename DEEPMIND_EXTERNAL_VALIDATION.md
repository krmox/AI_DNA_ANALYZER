# DeepMind External Dataset Validation — Audit Report

Status: **NOT SUITABLE — no experiment run**

This is an additional, standalone audit requested alongside the main FINAL_RESEARCH_REPORT.md.
It is **not** one of the six consolidated benchmark rounds and must not be merged with them.
No cascade run occurred; no CSV rows were added; no data was downloaded.

---

## 1. Dataset

Two categories of candidate were evaluated. Neither qualifies as a usable, genuinely-DeepMind,
BAM-based variant-calling benchmark dataset.

### Candidate A — "Google Brain Genomics Sequencing Dataset for Benchmarking and Development"
- Source: AWS Registry of Open Data — `registry.opendata.aws/google-brain-genomics-public/`
- Content: nine human samples as three parent-child trios; WGS at 40x/30x/20x and WES at
  100x/75x/50x; Illumina short-read for all samples, PacBio HiFi for four samples; FASTQ on S3,
  BAM/VCF at `gs://google-brain-genomics-public`.
- License: CC0 1.0 (public domain).
- Truth: references Genome in a Bottle (GIAB) truth sets.

### Candidate B — AlphaGenome / AlphaMissense (genuine Google DeepMind products)
- Source: `github.com/google-deepmind/alphagenome`, `deepmind.google/science/alphagenome/`,
  DeepMind blog (AlphaGenome Atlas, published ~2026-09-08).
- Content: precomputed functional/pathogenicity impact scores (AlphaGenome Variant Impact score)
  for all ~9 billion possible single-nucleotide substitutions across the reference genome; a
  1-petabyte scored atlas. AlphaMissense (2023) covers protein-coding missense variants
  specifically.

---

## 2. Provenance

**Candidate A is not a DeepMind dataset.** The AWS registry entry attributes it explicitly to the
**Google Brain group**, with the associated publication authored by Baid, Nattestad, Kolesnikov,
Goel, Yang, Chang, and Carroll (2020) — the DeepVariant team, part of Google Research/Google
Health, not DeepMind. Google Brain and DeepMind merged organizationally into "Google DeepMind" in
2023, but that merger does not retroactively make a 2020 Google Brain dataset a DeepMind-authored
one, and the dataset predates the merger by three years. Labeling this a "DeepMind dataset" would
misrepresent its provenance.

**Candidate B is genuinely Google DeepMind-authored** (AlphaGenome, AlphaMissense), but see §3 —
it is not structurally usable for this task.

No other public dataset attributable to DeepMind/Google DeepMind with raw sequencing reads (BAM/
CRAM/FASTQ) and a GIAB-compatible small-variant truth set was found in this search pass.

---

## 3. Scientific compatibility

| Check | Candidate A (Google Brain trios) | Candidate B (AlphaGenome/AlphaMissense) |
|---|---|---|
| GRCh38 or compatible reference | Likely yes (not confirmed without downloading metadata) | Yes (reference-coordinate scored atlas) |
| High-confidence truth set | Yes — cites GIAB truth sets | No — not a calling truth set; it's a predicted functional-impact score, not TP/FP/FN ground truth for read-based calling |
| Can derive TP/FP/FN/F1 | Plausible in principle | **No** — there are no raw reads and no calling-style truth to score against |
| Truth format compatible with our scoring | Would need to confirm VCF/BED format match | N/A — no calling truth exists |
| BAM/CRAM or raw reads available | Yes (BAM/VCF on GCS, FASTQ on S3) | **No** — no raw reads, no alignment data at all |
| Can run frozen cascade | Technically yes, if reads obtained | **No** — cascade requires BAM pileups; AlphaGenome/AlphaMissense provide no read-level input |
| Meaningful PB comparison possible | Plausible in principle | **No** |

**Candidate B verdict: `DEEPMIND DATASET = NOT SUITABLE`.** AlphaGenome and AlphaMissense are
variant-effect prediction models over the reference genome, not read-level variant callers, and
supply no BAM/CRAM/FASTQ and no calling-style truth set. They cannot be scored for TP/FP/FN
against sample reads under this protocol. Adapting the protocol to use predicted pathogenicity
scores as a substitute for calling truth would change the scientific question being asked and was
explicitly avoided.

**Candidate A verdict: technically closer to compatible, but fails independence (see §7) and its
provenance is Google Brain, not DeepMind** — so even a positive technical result here would not
constitute a "DeepMind validation" as requested. Given that failure, no download or cascade run
was attempted for Candidate A either; pursuing it would not answer the question this audit was
commissioned to answer.

---

## 4. Protocol

Not applicable — no experiment was run for either candidate. See §2–3 for why each was excluded
before reaching a disk-budget or execution step.

---

## 5. Results

Not applicable — no cascade run, no PB-only run, no comparison performed.

---

## 6. Statistical analysis

Not applicable.

---

## 7. Failure analysis

Not applicable — no per-locus results were generated.

---

## 8. Independence

This is the decisive finding for Candidate A, so it is documented in detail even though no run
occurred:

- Candidate A explicitly cites **GIAB truth sets** and is built around **parent-child trios**,
  which is exactly the structure of the GIAB Ashkenazi trio (HG002/HG003/HG004) already used as
  the primary benchmark population in this project (see FINAL_RESEARCH_REPORT.md).
- The AWS registry entry does not list explicit HG sample IDs, so it cannot be confirmed from the
  metadata alone whether Candidate A uses the *same* HG002/HG003/HG004 individuals, a different
  GIAB trio (e.g. HG005/HG006/HG007 — note HG005 is the Chinese trio son, already separately
  identified as the target for the blocked independent-sample experiment), or genuinely different
  individuals. Resolving this would require downloading and inspecting file headers/sample
  metadata, which was not done given Candidate A's disqualifying provenance issue (§2) and the
  disk-budget concerns already documented in FINAL_RESEARCH_REPORT.md (13 GB free / 94 GB).
- Sequencing platform overlaps with the project's existing protocol (Illumina short-read; some
  PacBio HiFi), so even in the best case, platform alone would not establish independence — batch,
  aligner, and individual identity would all need separate confirmation.
- **Conclusion: Candidate A cannot be called "independent" of HG002/HG003/HG004 without further
  verification, and given the provenance disqualification in §2, that verification was not
  pursued.**

Candidate B has no independence question to evaluate — it was excluded on structural grounds
(§3) before independence became relevant.

---

## 9. Limitations

- This audit is a literature/registry search, not an exhaustive survey of every dataset DeepMind
  or Google DeepMind has ever released; a dataset meeting the criteria could exist and not have
  surfaced in this search.
- Candidate A was not downloaded or inspected at the file level; the "likely overlaps with
  existing GIAB trio" conclusion in §8 is inferred from metadata (trio structure + GIAB truth
  citation), not confirmed by reading sample headers.
- No disk-budget calculation was performed for Candidate A beyond the general 13 GB/94 GB baseline
  already recorded in FINAL_RESEARCH_REPORT.md, since the dataset was disqualified on provenance
  and independence grounds before a download decision was needed.

---

## 10. Interpretation

No dataset was found that is both (a) genuinely attributable to DeepMind or Google DeepMind and
(b) structurally usable for a BAM-based small-variant calling benchmark against the frozen
cascade. The one dataset with the right *shape* (reads + GIAB truth) is Google Brain-authored, not
DeepMind-authored, and its trio structure makes independence from the existing HG002/HG003/HG004
benchmark population doubtful without further verification that was not pursued. The two genuine
DeepMind products found (AlphaGenome, AlphaMissense) are variant-effect predictors without raw
reads or calling-style truth, and are not applicable to this protocol at all.

**Formulation for external use, per instructions:** no additional evaluation was performed on a
DeepMind/Google DeepMind dataset in this pass, because no suitable dataset meeting both the
provenance and technical-compatibility requirements was identified. This should not be described
as "DeepMind validated the model" or "DeepMind data confirmed the results" — no such claim is
supported, positively or negatively, by this audit.

**Final classification: `DEEPMIND DATASET = NOT SUITABLE`** (Candidate B, structural) and
**`BLOCKED — DATA ACCESS / INDEPENDENCE UNVERIFIED`** (Candidate A, provenance + independence,
not pursued further).
