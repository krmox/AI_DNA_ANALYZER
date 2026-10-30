# DeepVariant comparison (separate table; descriptive only; same hap.py 0.3.15 / truth / BED / reference / BAM / seed 101 for all callers)

No non-inferiority or superiority claim vs DeepVariant is preregistered. CIs: 2-Mb block bootstrap, single seed. HG001/HG007 are previously evaluated samples.

| sample | contig | depth | caller | TP | FP | FN | precision | recall | F1 | F1 95% CI | wall s | CPU s | peak RSS MB |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| HG001 | chr19 | 30× | R2 | 66790 | 369 | 552 | 0.99451 | 0.99180 | 0.99315 | [0.99069, 0.99525] | 12 | 86 | 199 |
| HG001 | chr19 | 30× | R1.1 | 66680 | 337 | 662 | 0.99497 | 0.99017 | 0.99256 | [0.99005, 0.99470] | 12 | 86 | 198 |
| HG001 | chr19 | 30× | bcftools 1.21 | 66574 | 309 | 768 | 0.99538 | 0.98860 | 0.99198 | [0.98953, 0.99414] | 68 | 494 | 102 |
| HG001 | chr19 | 30× | DeepVariant 1.10.0 | 66976 | 228 | 366 | 0.99661 | 0.99456 | 0.99559 | [0.99317, 0.99779] | 1688 | 10659 | 9214 |
| HG007 | chr19 | 30× | R2 | 67460 | 312 | 531 | 0.99540 | 0.99219 | 0.99379 | [0.99239, 0.99507] | 14 | 86 | 195 |
| HG007 | chr19 | 30× | R1.1 | 67288 | 267 | 703 | 0.99605 | 0.98966 | 0.99284 | [0.99123, 0.99432] | 15 | 86 | 197 |
| HG007 | chr19 | 30× | bcftools 1.21 | 67258 | 241 | 733 | 0.99643 | 0.98922 | 0.99281 | [0.99096, 0.99446] | 79 | 499 | 103 |
| HG007 | chr19 | 30× | DeepVariant 1.10.0 | 67690 | 172 | 301 | 0.99747 | 0.99557 | 0.99652 | [0.99498, 0.99786] | 1879 | 10885 | 8935 |
