# Final factual audit of PAPER_MANUSCRIPT_V2.md (2026-09-26)

Scope: the rewritten manuscript (V1 → V2 → V2.x, final HG003 chr8 evaluation, chr20 300× comparison).

* **Numbers.** `audit_manuscript_numbers.py` checks 188 values of the new sections (Tables 11–16, chr8 result, bootstrap intervals, frame audit, depth, timing) against the CSV/JSON/hap.py artefacts: 188/188 pass. V1-era numbers (Tables 3–9, Figures 2–5) were carried over unchanged from the previously audited text.
* **Pronouns.** No "I/we/our/my" in the manuscript (regex scan).
* **Strong words.** "accurate/fast/efficient/generalizes/competitive/scalable/robust" do not occur as claims. "Faster" appears only in negations or in a question; runtime comparison with DeepVariant/GATK is stated as conditional and unmatched. "AI" is explained as a project name; no neural network exists in any version; V2.x is described as one logistic regression.
* **Weakened/added statements.** chr8 result is labelled "consistent with" the holdout, not a demonstrated gain (no V2 control was run); it is not called a biological replication. DeepVariant and GATK scoring higher than V2.x is stated explicitly. Platform/library similarity of HG003 was reduced to what the BAM header and directory name verify. Whole-chromosome V2.x numbers on chr20 are flagged as mixed development+holdout.
* **Known soft spots (kept, disclosed).** V1→V2 speedup ratios use a V1 record from another day (±15 %); chr8 runtime is a single cold-ish run; V2.x code is uncommitted; HG005 numbers remain record-only.
