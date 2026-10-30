## Benchmark: HG004 chr20, GIAB v4.2.1 confident regions, hap.py (xcmp)

| Caller | Wall (s) | CPU (s) | Peak RSS (MB) | SNP F1 | INDEL F1 | INDEL Precision | INDEL Recall | INDEL TP | INDEL FP | INDEL FN | GT concordance |
|---|---|---|---|---|---|---|---|---|---|---|---|
| val_baseline_30x | 20.4 (n=3) | 109.3 | 484 | N/A (no SNP calls) | 0.8918 | 0.8560 | 0.9307 | 9848 | 1657 | 733 | 0.9918 (82 GT err) |
| val_A_30x | 19.2 (n=3) | 105.4 | 500 | N/A (no SNP calls) | 0.9609 | 0.9594 | 0.9623 | 10182 | 445 | 399 | 0.9933 (69 GT err) |
| val_D_30x | 21.5 (n=3) | 121.4 | 491 | N/A (no SNP calls) | 0.9637 | 0.9689 | 0.9586 | 10143 | 337 | 438 | 0.9967 (34 GT err) |
| val_bcftools_30x | 306.7 (n=1) | 325.7 | 102 | 0.9905 | 0.9552 | 0.9619 | 0.9486 | 10037 | 417 | 544 | 0.9772 (234 GT err) |

## INDEL (PASS) by length class, genotype and repeat context

| Caller | Stratum | Truth | Recall | Precision | F1 |
|---|---|---|---|---|---|
| val_baseline_30x | homopolymer_ge7 | 4948 | 0.9194 | 0.8116 | 0.8621 |
| val_baseline_30x | not_in_TR_or_HP | 3722 | 0.9860 | 0.9745 | 0.9802 |
| val_baseline_30x | tandem_repeat_51to200 | 549 | 0.8288 | 0.7423 | 0.7831 |
| val_baseline_30x | tandem_repeat_le50 | 1172 | 0.8498 | 0.7568 | 0.8006 |
| val_baseline_30x | D16_PLUS | 96 | 0.8958 | 0.9149 | 0.9053 |
| val_baseline_30x | D1_5 | 4941 | 0.9405 | 0.8508 | 0.8934 |
| val_baseline_30x | D6_15 | 477 | 0.8386 | 0.8130 | 0.8256 |
| val_baseline_30x | I16_PLUS | 94 | 0.7340 | 0.7841 | 0.7582 |
| val_baseline_30x | I1_5 | 4808 | 0.9226 | 0.8789 | 0.9003 |
| val_baseline_30x | I6_15 | 495 | 0.8162 | 0.7580 | 0.7860 |
| val_A_30x | homopolymer_ge7 | 4948 | 0.9598 | 0.9500 | 0.9549 |
| val_A_30x | not_in_TR_or_HP | 3722 | 0.9844 | 0.9860 | 0.9852 |
| val_A_30x | tandem_repeat_51to200 | 549 | 0.9071 | 0.9250 | 0.9160 |
| val_A_30x | tandem_repeat_le50 | 1172 | 0.9283 | 0.9368 | 0.9325 |
| val_A_30x | D16_PLUS | 96 | 0.9583 | 0.9894 | 0.9736 |
| val_A_30x | D1_5 | 4941 | 0.9658 | 0.9557 | 0.9607 |
| val_A_30x | D6_15 | 477 | 0.9581 | 0.9525 | 0.9553 |
| val_A_30x | I16_PLUS | 94 | 0.8936 | 0.9767 | 0.9333 |
| val_A_30x | I1_5 | 4808 | 0.9563 | 0.9649 | 0.9606 |
| val_A_30x | I6_15 | 495 | 0.9596 | 0.9425 | 0.9510 |
| val_D_30x | homopolymer_ge7 | 4948 | 0.9553 | 0.9657 | 0.9605 |
| val_D_30x | not_in_TR_or_HP | 3722 | 0.9841 | 0.9863 | 0.9852 |
| val_D_30x | tandem_repeat_51to200 | 549 | 0.8925 | 0.9319 | 0.9118 |
| val_D_30x | tandem_repeat_le50 | 1172 | 0.9215 | 0.9492 | 0.9351 |
| val_D_30x | D16_PLUS | 96 | 0.9583 | 0.9894 | 0.9736 |
| val_D_30x | D1_5 | 4941 | 0.9595 | 0.9667 | 0.9631 |
| val_D_30x | D6_15 | 477 | 0.9602 | 0.9585 | 0.9593 |
| val_D_30x | I16_PLUS | 94 | 0.8936 | 0.9767 | 0.9333 |
| val_D_30x | I1_5 | 4808 | 0.9536 | 0.9739 | 0.9636 |
| val_D_30x | I6_15 | 495 | 0.9636 | 0.9483 | 0.9559 |
| val_bcftools_30x | homopolymer_ge7 | 4948 | 0.9402 | 0.9565 | 0.9483 |
| val_bcftools_30x | not_in_TR_or_HP | 3722 | 0.9839 | 0.9921 | 0.9880 |
| val_bcftools_30x | tandem_repeat_51to200 | 549 | 0.8561 | 0.8940 | 0.8746 |
| val_bcftools_30x | tandem_repeat_le50 | 1172 | 0.9155 | 0.9299 | 0.9227 |
| val_bcftools_30x | D16_PLUS | 96 | 0.8229 | 0.8404 | 0.8316 |
| val_bcftools_30x | D1_5 | 4941 | 0.9468 | 0.9657 | 0.9561 |
| val_bcftools_30x | D6_15 | 477 | 0.9099 | 0.9165 | 0.9132 |
| val_bcftools_30x | I16_PLUS | 94 | 0.7660 | 0.8554 | 0.8082 |
| val_bcftools_30x | I1_5 | 4808 | 0.9555 | 0.9694 | 0.9624 |
| val_bcftools_30x | I6_15 | 495 | 0.9253 | 0.9332 | 0.9292 |
