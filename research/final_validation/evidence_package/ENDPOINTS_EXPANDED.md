# Expanded benchmark — endpoint evaluation (rules: expanded_benchmark_preregistration.md)

HG001 and HG007 are **previously evaluated** samples; this is NOT an independent-genome validation. Δ = absolute F1 difference. CIs = 2-Mb block bootstrap within 2 contigs/sample; pooled = 4 cells of 2 samples (not a population interval).

## P1 — gain over R1.1 at 30× (PASS needs pooled Δ ≥ +0.003 and CI lower > 0)

* HG001: Δ=+0.00109, 95% CI [+0.00043, +0.00210] (2 cells, 6 cell-seeds) → **P1 FAIL** — any improvement (P1b): PASS
* HG007: Δ=+0.00160, 95% CI [+0.00072, +0.00290] (2 cells, 6 cell-seeds) → **P1 FAIL** — any improvement (P1b): PASS
* POOLED(4 cells): Δ=+0.00134, 95% CI [+0.00073, +0.00211] (4 cells, 12 cell-seeds) → **P1 FAIL** — any improvement (P1b): PASS

## P2 — non-inferiority to bcftools (PASS if CI lower ≥ −0.0005); superiority if lower > 0

### 10×
* HG001: Δ=+0.03496, 95% CI [+0.03306, +0.03671] (2 cells, 6 cell-seeds) → **P2 PASS**, superior (lower>0)
* HG007: Δ=+0.03531, 95% CI [+0.03315, +0.03744] (2 cells, 6 cell-seeds) → **P2 PASS**, superior (lower>0)
* POOLED(4 cells): Δ=+0.03513, 95% CI [+0.03374, +0.03655] (4 cells, 12 cell-seeds) → **P2 PASS**, superior (lower>0)
### 15×
* HG001: Δ=+0.01308, 95% CI [+0.01210, +0.01418] (2 cells, 6 cell-seeds) → **P2 PASS**, superior (lower>0)
* HG007: Δ=+0.01372, 95% CI [+0.01217, +0.01563] (2 cells, 6 cell-seeds) → **P2 PASS**, superior (lower>0)
* POOLED(4 cells): Δ=+0.01340, 95% CI [+0.01246, +0.01450] (4 cells, 12 cell-seeds) → **P2 PASS**, superior (lower>0)
### 20×
* HG001: Δ=+0.00584, 95% CI [+0.00491, +0.00725] (2 cells, 6 cell-seeds) → **P2 PASS**, superior (lower>0)
* HG007: Δ=+0.00649, 95% CI [+0.00502, +0.00860] (2 cells, 6 cell-seeds) → **P2 PASS**, superior (lower>0)
* POOLED(4 cells): Δ=+0.00616, 95% CI [+0.00525, +0.00738] (4 cells, 12 cell-seeds) → **P2 PASS**, superior (lower>0)
### 30×
* HG001: Δ=+0.00229, 95% CI [+0.00114, +0.00423] (2 cells, 6 cell-seeds) → **P2 PASS**, superior (lower>0)
* HG007: Δ=+0.00268, 95% CI [+0.00109, +0.00514] (2 cells, 6 cell-seeds) → **P2 PASS**, superior (lower>0)
* POOLED(4 cells): Δ=+0.00249, 95% CI [+0.00141, +0.00398] (4 cells, 12 cell-seeds) → **P2 PASS**, superior (lower>0)
### 50×
* HG001: Δ=+0.00147, 95% CI [+0.00018, +0.00352] (2 cells, 2 cell-seeds) → **P2 PASS**, superior (lower>0)
* HG007: Δ=+0.00202, 95% CI [+0.00034, +0.00464] (2 cells, 2 cell-seeds) → **P2 PASS**, superior (lower>0)
* POOLED(4 cells): Δ=+0.00175, 95% CI [+0.00059, +0.00341] (4 cells, 4 cell-seeds) → **P2 PASS**, superior (lower>0)

## P4 — no low-depth regression vs R1.1 (point Δ ≥ −0.002)

* 10× HG001: Δ=+0.00071, 95% CI [+0.00034, +0.00130] (2 cells, 6 cell-seeds) → **P4 PASS**
* 10× HG007: Δ=+0.00116, 95% CI [+0.00062, +0.00196] (2 cells, 6 cell-seeds) → **P4 PASS**
* 10× POOLED(4 cells): Δ=+0.00094, 95% CI [+0.00057, +0.00141] (4 cells, 12 cell-seeds) → **P4 PASS**
* 15× HG001: Δ=+0.00125, 95% CI [+0.00071, +0.00209] (2 cells, 6 cell-seeds) → **P4 PASS**
* 15× HG007: Δ=+0.00195, 95% CI [+0.00118, +0.00308] (2 cells, 6 cell-seeds) → **P4 PASS**
* 15× POOLED(4 cells): Δ=+0.00160, 95% CI [+0.00109, +0.00229] (4 cells, 12 cell-seeds) → **P4 PASS**

## P3 — absolute level (mean over seeds per cell; R2 ≥ 0.990; bcftools ceiling check)

| depth | sample | contig | R2 F1 | R1.1 F1 | bcftools F1 | R2 ≥ 0.990 | bcftools ≥ 0.990 | claim |
|---|---|---|---|---|---|---|---|---|
| 30× | HG001 | chr17 | 0.99217 | 0.99054 | 0.98886 | True | False | PASS (reached by R2 only; data-limited claim not allowed) |
| 30× | HG001 | chr19 | 0.99325 | 0.99271 | 0.99198 | True | True | PASS (ceiling check ok) |
| 30× | HG007 | chr17 | 0.99073 | 0.98847 | 0.98654 | True | False | PASS (reached by R2 only; data-limited claim not allowed) |
| 30× | HG007 | chr19 | 0.99381 | 0.99287 | 0.99263 | True | True | PASS (ceiling check ok) |
| 50× | HG001 | chr17 | 0.99261 | 0.99140 | 0.99004 | True | True | PASS (ceiling check ok) |
| 50× | HG001 | chr19 | 0.99421 | 0.99400 | 0.99384 | True | True | PASS (ceiling check ok) |
| 50× | HG007 | chr17 | 0.99201 | 0.99008 | 0.98842 | True | False | PASS (reached by R2 only; data-limited claim not allowed) |
| 50× | HG007 | chr19 | 0.99553 | 0.99511 | 0.99508 | True | True | PASS (ceiling check ok) |

## Mean SNP F1 per cell/depth (seeds averaged; n seeds)

| sample | contig | depth | n seeds | R2 | R1.1 | bcftools | provenance |
|---|---|---|---|---|---|---|---|
| HG001 | chr17 | 10× | 3 | 0.94634 | 0.94539 | 0.91324 | EXISTING |
| HG001 | chr17 | 15× | 3 | 0.97821 | 0.97660 | 0.96527 | NEW |
| HG001 | chr17 | 20× | 3 | 0.98800 | 0.98624 | 0.98156 | EXISTING |
| HG001 | chr17 | 30× | 3 | 0.99217 | 0.99054 | 0.98886 | EXISTING |
| HG001 | chr17 | 50× | 1 | 0.99261 | 0.99140 | 0.99004 | NEW |
| HG001 | chr19 | 10× | 3 | 0.94150 | 0.94103 | 0.90468 | EXISTING |
| HG001 | chr19 | 15× | 3 | 0.97729 | 0.97640 | 0.96408 | NEW |
| HG001 | chr19 | 20× | 3 | 0.98787 | 0.98697 | 0.98262 | EXISTING |
| HG001 | chr19 | 30× | 3 | 0.99325 | 0.99271 | 0.99198 | EXISTING |
| HG001 | chr19 | 50× | 1 | 0.99421 | 0.99400 | 0.99384 | NEW |
| HG007 | chr17 | 10× | 3 | 0.94242 | 0.94083 | 0.90828 | EXISTING |
| HG007 | chr17 | 15× | 3 | 0.97534 | 0.97276 | 0.96116 | NEW |
| HG007 | chr17 | 20× | 3 | 0.98578 | 0.98305 | 0.97814 | EXISTING |
| HG007 | chr17 | 30× | 3 | 0.99073 | 0.98847 | 0.98654 | EXISTING |
| HG007 | chr17 | 50× | 1 | 0.99201 | 0.99008 | 0.98842 | NEW |
| HG007 | chr19 | 10× | 3 | 0.94051 | 0.93978 | 0.90403 | EXISTING |
| HG007 | chr19 | 15× | 3 | 0.97607 | 0.97475 | 0.96281 | NEW |
| HG007 | chr19 | 20× | 3 | 0.98782 | 0.98628 | 0.98248 | EXISTING |
| HG007 | chr19 | 30× | 3 | 0.99381 | 0.99287 | 0.99263 | EXISTING |
| HG007 | chr19 | 50× | 1 | 0.99553 | 0.99511 | 0.99508 | NEW |

## Resource guardrails (R2 vs R1.1; median over seeds). Verdicts use NEW rows only (sequential, load recorded); EXISTING rows shown for context

| sample | contig | depth | prov | wall Δ% (Δs) | CPU Δ% | RSS Δ% | verdict |
|---|---|---|---|---|---|---|---|
| HG001 | chr17 | 10× | EXISTING | +0.2% (+0.0s) | -0.9% | -0.8% | within (context only) |
| HG001 | chr17 | 15× | NEW | +0.6% (+0.1s) | +0.7% | +1.1% | PASS |
| HG001 | chr17 | 20× | EXISTING | +6.7% (+0.9s) | +0.6% | +0.6% | within (context only) |
| HG001 | chr17 | 30× | EXISTING | +1.7% (+0.3s) | +0.8% | -0.7% | within (context only) |
| HG001 | chr17 | 50× | NEW | +2.2% (+0.5s) | +1.0% | +2.1% | PASS |
| HG001 | chr19 | 10× | EXISTING | +1.7% (+0.1s) | +0.3% | +0.1% | within (context only) |
| HG001 | chr19 | 15× | NEW | +0.9% (+0.1s) | +0.2% | +0.3% | PASS |
| HG001 | chr19 | 20× | EXISTING | +3.2% (+0.3s) | +0.2% | +1.2% | within (context only) |
| HG001 | chr19 | 30× | EXISTING | +0.1% (+0.0s) | +0.5% | +0.9% | within (context only) |
| HG001 | chr19 | 50× | NEW | +18.2% (+3.4s) | -0.7% | -1.4% | FAIL |
| HG007 | chr17 | 10× | EXISTING | -3.7% (-0.4s) | -0.1% | +0.8% | within (context only) |
| HG007 | chr17 | 15× | NEW | +1.4% (+0.2s) | +0.6% | +0.2% | PASS |
| HG007 | chr17 | 20× | EXISTING | +1.7% (+0.2s) | +0.9% | +0.0% | within (context only) |
| HG007 | chr17 | 30× | EXISTING | +0.2% (+0.0s) | +0.7% | +0.2% | within (context only) |
| HG007 | chr17 | 50× | NEW | +0.9% (+0.2s) | +0.7% | +0.6% | PASS |
| HG007 | chr19 | 10× | EXISTING | +4.4% (+0.3s) | +0.3% | -0.2% | within (context only) |
| HG007 | chr19 | 15× | NEW | +0.7% (+0.1s) | +0.2% | -0.9% | PASS |
| HG007 | chr19 | 20× | EXISTING | +0.8% (+0.1s) | +0.1% | -0.1% | within (context only) |
| HG007 | chr19 | 30× | EXISTING | -2.4% (-0.4s) | +0.4% | -0.9% | within (context only) |
| HG007 | chr19 | 50× | NEW | -0.1% (-0.0s) | -0.4% | -0.5% | PASS |

Metric consistency check (F1 from TP/FP/FN vs hap.py F1, tol 2e-4): 0 inconsistencies
