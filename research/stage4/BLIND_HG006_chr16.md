# R2 blind-new: HG006 chr16 (hap.py 0.3.15, SNP PASS; point estimates, bootstrap CIs NOT yet computed)
Seeds 102-103 for 10/15/20/50x not run (reduced schedule, PREREG addendum 2).

| cond | n | R1.1 F1 | R2 F1 | bcftools F1 | R2-R1.1 | R2-bcftools |
|---|---|---|---|---|---|---|
| 10x | 1 | 0.94467 | 0.94600 | 0.91525 | +0.0013 | +0.0308 |
| 15x | 1 | 0.97448 | 0.97615 | 0.96419 | +0.0017 | +0.0120 |
| 20x | 1 | 0.98388 | 0.98596 | 0.97977 | +0.0021 | +0.0062 |
| 30x | 3 | 0.98878 | 0.99029 | 0.98681 | +0.0015 | +0.0035 |
| 50x | 1 | 0.99027 | 0.99131 | 0.98842 | +0.0010 | +0.0029 |
| native ~106x | 1 | 0.99140 | 0.99132 | 0.98907 | -0.0001 | +0.0023 |

30x per seed R2: 0.99056 / 0.99029 / 0.99003.
Point-estimate verdicts (CI pending): P1 (>= +0.003 at 30x) NOT met (+0.0015). P2 met at every depth (point). P3 F1>=0.99 reached at 30x/50x/native where bcftools does not at 30x/50x (bcftools reaches only native 0.9891 <0.99). P4 met.
Cost: CPU R2 ~ R1.1 (+-2%), RSS +~9 MB; wall noisy under load (machine busy).
OOD guard did not trigger (read_len 148, eps_hat clamped at 0.002).
