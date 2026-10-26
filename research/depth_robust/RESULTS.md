# Stage 1 results (per PREREGISTRATION.md; all numbers from hap.py 0.3.15 SNP PASS, raw outputs in /mnt/archive/AI_DNA_ANALYZER_depth_robust/happy/out)
| Set | Method | TP | FN | FP | P | R | F1 |
|---|---|---|---|---|---|---|---|
| chr8 HG003 30x | V2 control | 182592 | 3930 | 2686 | 0.9855 | 0.9789 | 0.98220 |
| chr8 HG003 30x | old A+C | 182644 | 3878 | 5065 | 0.9730 | 0.9792 | 0.97610 |
| chr8 HG003 30x | R1 (finalist) | 182310 | 4212 | 441 | 0.9976 | 0.9774 | 0.98740 |
| chr8 HG003 30x | R2 | 182395 | 4127 | 493 | 0.9973 | 0.9779 | 0.98749 |
| chr8 HG003 30x | bcftools QUAL>=30 | 184706 | 1816 | 517 | 0.9972 | 0.9903 | 0.99372 |
| chr20 HG002 odd blocks 30x (seed 7) | V2 control | 33339 | 691 | 876 | 0.9744 | 0.9797 | 0.97704 |
| chr20 odd 30x | old A+C | 33342 | 688 | 1673 | 0.9522 | 0.9798 | 0.96580 |
| chr20 odd 30x | R1 | 33291 | 739 | 115 | 0.9966 | 0.9783 | 0.98734 |
| chr20 odd 30x | R2 | 33308 | 722 | 131 | 0.9961 | 0.9788 | 0.98736 |
| chr20 odd 30x | bcftools QUAL>=30 | 33674 | 356 | 176 | 0.9948 | 0.9895 | 0.99216 |
| chr8 HG003 300x | V2 control | 182959 | 3563 | 3902 | 0.9791 | 0.9809 | 0.98001 |
| chr8 HG003 300x | old A+C (earlier single run) | 182904 | 3618 | 446 | 0.9976 | 0.9806 | 0.98901 |
| chr8 HG003 300x | R1 | 182613 | 3909 | 302 | 0.9983 | 0.9790 | 0.98860 |
| chr8 HG003 300x | R2 | 182525 | 3997 | 345 | 0.9981 | 0.9786 | 0.98825 |

Paired block bootstrap (bootstrap_30x.txt): chr8 30x R1-control +0.00519 [+0.00432,+0.00625]; chr20 odd 30x R1-control +0.01030 [+0.00723,+0.01393].
Endpoints: primary (chr8 30x, >= +0.010) NOT MET (+0.0052). Secondary 1 (chr8 300x R1 vs old A+C within -0.002): MET (-0.0004). Secondary 2 (chr20 30x, >= +0.010, CI>0): MET (+0.0103).
bcftools at 300x on chr8 was not run (about 3.5 h single thread).
