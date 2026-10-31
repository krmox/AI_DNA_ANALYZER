## Benchmark: HG004 chr20, GIAB v4.2.1 confident regions, hap.py (xcmp)

| Caller | Wall (s) | CPU (s) | Peak RSS (MB) | SNP F1 | INDEL F1 | INDEL Precision | INDEL Recall | INDEL TP | INDEL FP | INDEL FN |
|---|---|---|---|---|---|---|---|---|---|---|
| bcftools_15x | 371.1 (n=1) | 339.4 | 110 | 0.9692 | 0.8857 | 0.9242 | 0.8502 | 9352 | 799 | 1648 |
| nindels_15x | 12.7 (n=3) | 71.0 | 329 | N/A (no SNP calls) | 0.8491 | 0.8274 | 0.8719 | 9591 | 2001 | 1409 |
| bcftools_30x | 660.5 (n=1) | 581.0 | 110 | 0.9911 | 0.9330 | 0.9483 | 0.9181 | 10099 | 578 | 901 |
| nindels_30x | 29.2 (n=3) | 160.7 | 502 | N/A (no SNP calls) | 0.8743 | 0.8334 | 0.9194 | 10113 | 2023 | 887 |
| bcftools_67x | 1019.8 (n=1) | 1019.2 | 111 | 0.9926 | 0.9548 | 0.9627 | 0.9470 | 10417 | 426 | 583 |
| nindels_67x | 64.9 (n=3) | 367.5 | 960 | N/A (no SNP calls) | 0.8390 | 0.7678 | 0.9248 | 10173 | 3078 | 827 |

## INDEL (PASS) by length class, genotype and repeat context

| Caller | Stratum | Truth | Recall | Precision | F1 |
|---|---|---|---|---|---|
| bcftools_15x | homopolymer_ge7 | 4908 | 0.7553 | 0.8918 | 0.8179 |
| bcftools_15x | not_in_TR_or_HP | 3876 | 0.9585 | 0.9792 | 0.9687 |
| bcftools_15x | tandem_repeat_51to200 | 653 | 0.8637 | 0.9035 | 0.8832 |
| bcftools_15x | tandem_repeat_le50 | 1449 | 0.8806 | 0.8925 | 0.8865 |
| bcftools_15x | D16_PLUS | 136 | 0.8676 | 0.9365 | 0.9008 |
| bcftools_15x | D1_5 | 5066 | 0.8279 | 0.9196 | 0.8713 |
| bcftools_15x | D6_15 | 610 | 0.8607 | 0.8748 | 0.8677 |
| bcftools_15x | I16_PLUS | 104 | 0.8654 | 0.8911 | 0.8780 |
| bcftools_15x | I1_5 | 4941 | 0.8498 | 0.9384 | 0.8919 |
| bcftools_15x | I6_15 | 557 | 0.8689 | 0.9057 | 0.8869 |
| nindels_15x | homopolymer_ge7 | 4908 | 0.8201 | 0.8260 | 0.8230 |
| nindels_15x | not_in_TR_or_HP | 3876 | 0.9727 | 0.8679 | 0.9173 |
| nindels_15x | tandem_repeat_51to200 | 653 | 0.7994 | 0.7543 | 0.7762 |
| nindels_15x | tandem_repeat_le50 | 1449 | 0.8130 | 0.7475 | 0.7788 |
| nindels_15x | D16_PLUS | 136 | 0.7941 | 0.8640 | 0.8276 |
| nindels_15x | D1_5 | 5066 | 0.8723 | 0.7795 | 0.8233 |
| nindels_15x | D6_15 | 610 | 0.7721 | 0.7734 | 0.7728 |
| nindels_15x | I16_PLUS | 104 | 0.8269 | 0.9149 | 0.8687 |
| nindels_15x | I1_5 | 4941 | 0.8646 | 0.9003 | 0.8821 |
| nindels_15x | I6_15 | 557 | 0.7504 | 0.7656 | 0.7579 |
| bcftools_30x | homopolymer_ge7 | 4908 | 0.8566 | 0.9144 | 0.8845 |
| bcftools_30x | not_in_TR_or_HP | 3876 | 0.9856 | 0.9940 | 0.9898 |
| bcftools_30x | tandem_repeat_51to200 | 653 | 0.9280 | 0.9373 | 0.9327 |
| bcftools_30x | tandem_repeat_le50 | 1449 | 0.9455 | 0.9493 | 0.9474 |
| bcftools_30x | D16_PLUS | 136 | 0.9338 | 0.9845 | 0.9585 |
| bcftools_30x | D1_5 | 5066 | 0.9066 | 0.9482 | 0.9269 |
| bcftools_30x | D6_15 | 610 | 0.9344 | 0.9199 | 0.9271 |
| bcftools_30x | I16_PLUS | 104 | 0.9038 | 0.9038 | 0.9038 |
| bcftools_30x | I1_5 | 4941 | 0.9184 | 0.9554 | 0.9365 |
| bcftools_30x | I6_15 | 557 | 0.9084 | 0.9204 | 0.9144 |
| nindels_30x | homopolymer_ge7 | 4908 | 0.8924 | 0.7815 | 0.8333 |
| nindels_30x | not_in_TR_or_HP | 3876 | 0.9912 | 0.9665 | 0.9787 |
| nindels_30x | tandem_repeat_51to200 | 653 | 0.8545 | 0.7480 | 0.7977 |
| nindels_30x | tandem_repeat_le50 | 1449 | 0.8461 | 0.7275 | 0.7823 |
| nindels_30x | D16_PLUS | 136 | 0.8529 | 0.8593 | 0.8561 |
| nindels_30x | D1_5 | 5066 | 0.9224 | 0.8074 | 0.8611 |
| nindels_30x | D6_15 | 610 | 0.8180 | 0.7463 | 0.7805 |
| nindels_30x | I16_PLUS | 104 | 0.8846 | 0.8679 | 0.8762 |
| nindels_30x | I1_5 | 4941 | 0.9168 | 0.8904 | 0.9034 |
| nindels_30x | I6_15 | 557 | 0.7756 | 0.7200 | 0.7468 |
| bcftools_67x | homopolymer_ge7 | 4908 | 0.9091 | 0.9376 | 0.9231 |
| bcftools_67x | not_in_TR_or_HP | 3876 | 0.9894 | 0.9948 | 0.9921 |
| bcftools_67x | tandem_repeat_51to200 | 653 | 0.9556 | 0.9609 | 0.9582 |
| bcftools_67x | tandem_repeat_le50 | 1449 | 0.9620 | 0.9698 | 0.9659 |
| bcftools_67x | D16_PLUS | 136 | 0.9632 | 0.9850 | 0.9740 |
| bcftools_67x | D1_5 | 5066 | 0.9410 | 0.9618 | 0.9513 |
| bcftools_67x | D6_15 | 610 | 0.9525 | 0.9447 | 0.9486 |
| bcftools_67x | I16_PLUS | 104 | 0.9231 | 0.9320 | 0.9275 |
| bcftools_67x | I1_5 | 4941 | 0.9474 | 0.9665 | 0.9568 |
| bcftools_67x | I6_15 | 557 | 0.9425 | 0.9584 | 0.9504 |
| nindels_67x | homopolymer_ge7 | 4908 | 0.9093 | 0.6904 | 0.7849 |
| nindels_67x | not_in_TR_or_HP | 3876 | 0.9910 | 0.9559 | 0.9731 |
| nindels_67x | tandem_repeat_51to200 | 653 | 0.8469 | 0.6844 | 0.7570 |
| nindels_67x | tandem_repeat_le50 | 1449 | 0.8323 | 0.6634 | 0.7383 |
| nindels_67x | D16_PLUS | 136 | 0.8750 | 0.8500 | 0.8623 |
| nindels_67x | D1_5 | 5066 | 0.9313 | 0.7181 | 0.8109 |
| nindels_67x | D6_15 | 610 | 0.8066 | 0.6819 | 0.7390 |
| nindels_67x | I16_PLUS | 104 | 0.8846 | 0.8070 | 0.8440 |
| nindels_67x | I1_5 | 4941 | 0.9249 | 0.8554 | 0.8888 |
| nindels_67x | I6_15 | 557 | 0.7648 | 0.6625 | 0.7100 |

## Error taxonomy: nindels_15x


**False negatives by cause**

| Cause | Count |
|---|---|
| never a candidate (discovery miss) | 502 |
| multi-allelic site (truth 1/2), alleles genotyped separately | 474 |
| genotype wrong (allele found, GT differs) | 188 |
| different allele/representation called nearby | 102 |
| candidate found but filtered (LowProb/UNCERTAIN) | 24 |

**False positives by cause and repeat context**

| Cause \| context | Count |
|---|---|
| no truth variant nearby (pure FP) \| other | 478 |
| allele/representation mismatch near a truth variant \| HP>=7 | 450 |
| no truth variant nearby (pure FP) \| HP>=7 | 267 |
| allele/representation mismatch near a truth variant \| TR<=50 | 201 |
| no truth variant nearby (pure FP) \| TR<=50 | 170 |
| genotype error (allele right, GT wrong) \| HP>=7 | 131 |
| allele/representation mismatch near a truth variant \| other | 100 |
| allele/representation mismatch near a truth variant \| TR51-200 | 86 |
| no truth variant nearby (pure FP) \| TR51-200 | 70 |
| genotype error (allele right, GT wrong) \| TR<=50 | 22 |
| genotype error (allele right, GT wrong) \| TR51-200 | 14 |
| genotype error (allele right, GT wrong) \| other | 12 |

**FP called genotype**

- no truth variant nearby, called 0/1: 692
- allele/representation mismatch near a truth variant, called 0/1: 413
- allele/representation mismatch near a truth variant, called 1/1: 380
- no truth variant nearby, called 1/1: 292
- genotype error, called 0/1: 116
- genotype error, called 1/1: 63
- allele/representation mismatch near a truth variant, called 2/1: 44
- no truth variant nearby, called 2/1: 1

**Recall: clustered (another truth variant within 10 bp) vs isolated**

- isolated: 9190/10352 = 0.8878
- clustered: 520/648 = 0.8025

**Median features: TP calls vs pure-FP calls**

| Feature | TP | pure FP |
|---|---|---|
| NALT | 8.00 | 2.00 |
| NREF | 4.00 | 8.00 |
| F | 0.62 | 0.22 |
| HQ | 74.67 | 15.28 |
| QMIN | 0.92 | 0.93 |

Pure FPs with NALT <= 2: 623/985

**Examples**

- *genotype wrong (allele found, GT differs)*
  - `chr20:217716 ATG>A truth GT 0/1 [TR51-200] -> ATG>A PASS F=0.25;NALT=2;NREF=6;NOTH=0;CR=0.36336826329242866;RPRI=0.185816040041`
  - `chr20:314905 CT>C truth GT 0/1 [HP>=7] -> CT>C PASS F=1;NALT=3;NREF=0;NOTH=1;CR=0.43012059070764175;RPRI=0.231185439371513`
  - `chr20:454236 CAA>C truth GT 0/1 [HP>=7] -> CAA>C UNCERTAIN F=1;NALT=2;NREF=0;NOTH=1;CR=0.4371581736513161;RPRI=0.6438829574040897`
  - `chr20:596536 G>GAAAACA truth GT 0/1 [HP>=7] -> G>GAAAACA PASS F=1;NALT=8;NREF=0;NOTH=1;CR=0.2187733899982377;RPRI=0.1093866950012681`
- *multi-allelic site (truth 1/2), alleles genotyped separately*
  - `chr20:232033 A>AAAATAAAT,AAAATAAATAAAT truth GT 1/2 [TR<=50] -> A>AAAATAAAT PASS F=1;NALT=4;NREF=0;NOTH=6;CR=0.28809519845678094;RPRI=0.144231795889787; A>AAAATAAATAAAT PASS F=1;NALT=5;NREF=0;NOTH=5;CR=0.28809519845678094;RPRI=0.144049572603514`
  - `chr20:348484 CGTGTGT>CGTGT,C truth GT 2/1 [TR<=50] -> CGT>C PASS F=1;NALT=4;NREF=0;NOTH=5;CR=0.368640649723186;RPRI=0.1845045215229898;; CGTGTGT>C PASS F=1;NALT=5;NREF=0;NOTH=4;CR=0.368640649723186;RPRI=0.1843222982367173;`
  - `chr20:511703 TACAC>CACAC,T truth GT 2/1 [TR<=50] -> TACAC>T PASS F=0.1875;NALT=3;NREF=13;NOTH=0;CR=0.21112491269261563;RPRI=0.107298289`
  - `chr20:571051 C>CT,CTTT truth GT 1/2 [HP>=7] -> C>CT PASS F=1;NALT=4;NREF=0;NOTH=3;CR=0.5443119885471114;RPRI=0.2723401909349525`
- *never a candidate (discovery miss)*
  - `chr20:255818 GATATAT>G truth GT 0/1 [TR<=50] -> no candidate`
  - `chr20:356342 CA>C truth GT 0/1 [HP>=7] -> no candidate`
  - `chr20:364763 C>CAA truth GT 0/1 [HP>=7] -> no candidate`
  - `chr20:365798 CTT>CT,C truth GT 2/1 [HP>=7] -> no candidate`
- *different allele/representation called nearby*
  - `chr20:367246 AC>A truth GT 0/1 [TR<=50] -> AAAC>A PASS F=0.6363636363636364;NALT=7;NREF=4;NOTH=0;CR=0.10839129895466793;RPRI=`
  - `chr20:367250 AAC>A truth GT 0/1 [TR<=50] -> AAAC>A PASS F=0.6363636363636364;NALT=7;NREF=4;NOTH=0;CR=0.10839129895466793;RPRI=`
  - `chr20:392557 A>AT truth GT 0/1 [HP>=7] -> A>AT PASS F=0.5714285714285714;NALT=4;NREF=3;NOTH=3;CR=0.3961044420673366;RPRI=0; AT>A PASS F=0.4;NALT=2;NREF=3;NOTH=5;CR=0.3928607445472544;RPRI=0.22763251491406`
  - `chr20:454295 G>GA truth GT 0/1 [HP>=7] -> G>GA PASS F=0.8888888888888888;NALT=8;NREF=1;NOTH=0;CR=0.2870592687291093;RPRI=0`
- *candidate found but filtered (LowProb/UNCERTAIN)*
  - `chr20:2224954 CTT>CT,C truth GT 2/1 [HP>=7] -> CT>C UNCERTAIN F=1;NALT=2;NREF=0;NOTH=0;CR=0.3947689395392092;RPRI=0.6226883403480363`
  - `chr20:3860373 A>AT,ATT truth GT 1/2 [HP>=7] -> A>AT UNCERTAIN F=1;NALT=2;NREF=0;NOTH=1;CR=0.5729435787201076;RPRI=0.7117756599384855`
  - `chr20:11543211 TACACACAC>TAC,T truth GT 2/1 [TR51-200] -> TACACAC>T UNCERTAIN F=1;NALT=3;NREF=0;NOTH=0;CR=0.3039767144659411;RPRI=0.1681135012506630`
  - `chr20:11769171 CTT>CT,C truth GT 2/1 [HP>=7] -> CT>C UNCERTAIN F=1;NALT=2;NREF=0;NOTH=1;CR=0.5234539240409751;RPRI=0.6870308325989193`
- *no truth variant nearby (pure FP) | HP>=7*
  - `chr20:92968 GCT>G called 0/1 truth ./. NALT=2 NREF=2 F=0.5 HQ=24.25370567397715`
  - `chr20:218313 GT>G called 0/1 truth ./. NALT=3 NREF=15 F=0.16666666666666666 HQ=12.034083515026282`
  - `chr20:736736 CA>C called 1/1 truth ./. NALT=11 NREF=0 F=1 HQ=126.84766398803845`
- *no truth variant nearby (pure FP) | other*
  - `chr20:214157 GT>G called 0/1 truth ./. NALT=2 NREF=12 F=0.14285714285714285 HQ=11.22985344896`
  - `chr20:291411 TC>T called 0/1 truth ./. NALT=2 NREF=7 F=0.2222222222222222 HQ=11.693366179495683`
  - `chr20:307400 TG>T called 0/1 truth ./. NALT=2 NREF=9 F=0.18181818181818182 HQ=10.525805240029342`
- *genotype error (allele right, GT wrong) | TR51-200*
  - `chr20:217716 ATG>A called 1/1 truth 0/1 NALT=2 NREF=6 F=0.25 HQ=27.189189696091773`
  - `chr20:10507904 T>TTA called 1/1 truth 0/1 NALT=10 NREF=0 F=1 HQ=107.02680841333262`
  - `chr20:13269934 A>AGGATGGATGGAT called 1/1 truth 0/1 NALT=3 NREF=1 F=0.75 HQ=58.453428044622456`
- *allele/representation mismatch near a truth variant | TR<=50*
  - `chr20:232033 A>AAAATAAAT called 1/1 truth ./. NALT=4 NREF=0 F=1 HQ=159.51820552796164`
  - `chr20:232033 A>AAAATAAATAAAT called 1/1 truth ./. NALT=5 NREF=0 F=1 HQ=119.6837424713915`
  - `chr20:367244 AAAC>A called 0/1 truth ./. NALT=7 NREF=4 F=0.6363636363636364 HQ=72.39246104683593`
- *no truth variant nearby (pure FP) | TR51-200*
  - `chr20:271443 TAC>T called 0/1 truth ./. NALT=3 NREF=9 F=0.25 HQ=19.257348446017204`
  - `chr20:2011231 CCTCTCTCT>C called 1/1 truth ./. NALT=10 NREF=0 F=1 HQ=225.3552454124325`
  - `chr20:2011231 CCTCTCTCTCTCTCTCT>C called 0/1 truth ./. NALT=3 NREF=0 F=1 HQ=45.78977560396197`
- *genotype error (allele right, GT wrong) | HP>=7*
  - `chr20:314905 CT>C called 1/1 truth 0/1 NALT=3 NREF=0 F=1 HQ=29.1069486719569`
  - `chr20:454295 G>GA called 1/1 truth 0/1 NALT=8 NREF=1 F=0.8888888888888888 HQ=66.847556119984`
  - `chr20:596536 G>GAAAACA called 1/1 truth 0/1 NALT=8 NREF=0 F=1 HQ=137.14279345144385`
- *no truth variant nearby (pure FP) | TR<=50*
  - `chr20:348484 CGT>C called 1/1 truth ./. NALT=4 NREF=0 F=1 HQ=59.069095022529154`
  - `chr20:348484 CGTGTGT>C called 0/1 truth ./. NALT=5 NREF=0 F=1 HQ=73.8675599286779`
  - `chr20:685696 AAATAATAAT>A called 1/1 truth ./. NALT=8 NREF=0 F=1 HQ=177.26885122861114`
- *allele/representation mismatch near a truth variant | HP>=7*
  - `chr20:385413 TA>T called 0/1 truth ./. NALT=2 NREF=7 F=0.2222222222222222 HQ=16.53147731202216`
  - `chr20:392557 AT>ATT,A called 2/1 truth ./. NALT=4 NREF=3 F=0.5714285714285714 HQ=32.122280465839225`
  - `chr20:458330 TA>T called 0/1 truth ./. NALT=2 NREF=8 F=0.2 HQ=28.04805384416511`
- *allele/representation mismatch near a truth variant | other*
  - `chr20:992033 GGAC>G called 0/1 truth ./. NALT=2 NREF=6 F=0.25 HQ=11.11528177118987`
  - `chr20:2501317 GA>G called 0/1 truth ./. NALT=2 NREF=12 F=0.14285714285714285 HQ=9.835710879777935`
  - `chr20:3696456 G>GGCT called 1/1 truth ./. NALT=2 NREF=17 F=0.10526315789473684 HQ=110.21029166374547`
- *genotype error (allele right, GT wrong) | other*
  - `chr20:1373848 GAGAC>G called 1/1 truth 0/1 NALT=9 NREF=0 F=1 HQ=138.28365231588404`
  - `chr20:2874203 GTCTCCCAGGCCCCTC>G called 1/1 truth 0/1 NALT=8 NREF=0 F=1 HQ=128.94866385112954`
  - `chr20:9025537 A>AAAAC called 1/1 truth 0/1 NALT=8 NREF=0 F=1 HQ=120.90383393135853`
- *allele/representation mismatch near a truth variant | TR51-200*
  - `chr20:1850385 CGT>C called 1/1 truth ./. NALT=4 NREF=0 F=1 HQ=36.30658582470923`
  - `chr20:2945582 GGTGTGTGTGT>GGT,G called 2/1 truth ./. NALT=3 NREF=6 F=0.3333333333333333 HQ=72.77588755125623`
  - `chr20:3310333 A>AAAAGAAAG called 1/1 truth ./. NALT=2 NREF=0 F=1 HQ=169.83454815875194`
- *genotype error (allele right, GT wrong) | TR<=50*
  - `chr20:1948230 TATTC>T called 1/1 truth 0/1 NALT=9 NREF=0 F=1 HQ=126.01917253578767`
  - `chr20:4411737 T>TTG called 0/1 truth 1/1 NALT=11 NREF=1 F=0.9166666666666666 HQ=103.01267905333711`
  - `chr20:5454164 AATATAT>A called 1/1 truth 0/1 NALT=10 NREF=1 F=0.9090909090909091 HQ=211.8334767465526`

## Error taxonomy: nindels_30x


**False negatives by cause**

| Cause | Count |
|---|---|
| multi-allelic site (truth 1/2), alleles genotyped separately | 484 |
| different allele/representation called nearby | 153 |
| genotype wrong (allele found, GT differs) | 114 |
| never a candidate (discovery miss) | 105 |
| candidate found but filtered (LowProb/UNCERTAIN) | 2 |

**False positives by cause and repeat context**

| Cause \| context | Count |
|---|---|
| allele/representation mismatch near a truth variant \| HP>=7 | 725 |
| no truth variant nearby (pure FP) \| HP>=7 | 400 |
| allele/representation mismatch near a truth variant \| TR<=50 | 255 |
| no truth variant nearby (pure FP) \| TR<=50 | 188 |
| allele/representation mismatch near a truth variant \| TR51-200 | 102 |
| genotype error (allele right, GT wrong) \| HP>=7 | 100 |
| no truth variant nearby (pure FP) \| other | 80 |
| no truth variant nearby (pure FP) \| TR51-200 | 78 |
| allele/representation mismatch near a truth variant \| other | 71 |
| genotype error (allele right, GT wrong) \| TR<=50 | 12 |
| genotype error (allele right, GT wrong) \| TR51-200 | 8 |
| genotype error (allele right, GT wrong) \| other | 4 |

**FP called genotype**

- allele/representation mismatch near a truth variant, called 0/1: 631
- allele/representation mismatch near a truth variant, called 1/1: 422
- no truth variant nearby, called 1/1: 384
- no truth variant nearby, called 0/1: 359
- allele/representation mismatch near a truth variant, called 2/1: 100
- genotype error, called 0/1: 92
- genotype error, called 1/1: 32
- no truth variant nearby, called 2/1: 3

**Recall: clustered (another truth variant within 10 bp) vs isolated**

- isolated: 9605/10352 = 0.9278
- clustered: 537/648 = 0.8287

**Median features: TP calls vs pure-FP calls**

| Feature | TP | pure FP |
|---|---|---|
| NALT | 15.00 | 5.00 |
| NREF | 9.00 | 0.00 |
| F | 0.58 | 1.00 |
| HQ | 143.40 | 89.07 |
| QMIN | 0.92 | 0.84 |

Pure FPs with NALT <= 2: 218/746

**Examples**

- *genotype wrong (allele found, GT differs)*
  - `chr20:217716 ATG>A truth GT 0/1 [TR51-200] -> ATG>A PASS F=0.35;NALT=7;NREF=13;NOTH=0;CR=0.31961826329242865;RPRI=0.15990113792`
  - `chr20:454236 CAA>C truth GT 0/1 [HP>=7] -> CAA>C PASS F=1;NALT=6;NREF=0;NOTH=2;CR=0.4788248403179828;RPRI=0.2394124407196778`
  - `chr20:454295 G>GA truth GT 0/1 [HP>=7] -> G>GA PASS F=0.8666666666666667;NALT=13;NREF=2;NOTH=0;CR=0.29261482428466484;RPRI`
  - `chr20:1122295 CTTTT>C truth GT 1/1 [HP>=7] -> CTTTT>C PASS F=0.8666666666666667;NALT=13;NREF=2;NOTH=1;CR=0.5599922581805352;RPRI=`
- *multi-allelic site (truth 1/2), alleles genotyped separately*
  - `chr20:232033 A>AAAATAAAT,AAAATAAATAAAT truth GT 1/2 [TR<=50] -> A>AAAATAAAT PASS F=1;NALT=6;NREF=0;NOTH=11;CR=0.28809519845678094;RPRI=0.14404761978907; A>AAAATAAATAAAT PASS F=1;NALT=9;NREF=0;NOTH=8;CR=0.28809519845678094;RPRI=0.144047599228412`
  - `chr20:348484 CGTGTGT>CGTGT,C truth GT 2/1 [TR<=50] -> CGT>C PASS F=1;NALT=6;NREF=0;NOTH=6;CR=0.3408628719454082;RPRI=0.1704314565333906; CGTGTGT>C PASS F=1;NALT=6;NREF=0;NOTH=6;CR=0.3408628719454082;RPRI=0.1704314565333906`
  - `chr20:365798 CTT>CT,C truth GT 2/1 [HP>=7] -> CTT>C PASS F=1;NALT=3;NREF=0;NOTH=1;CR=0.48763776212656507;RPRI=0.259944025080975`
  - `chr20:511703 TACAC>CACAC,T truth GT 2/1 [TR<=50] -> TACAC>T PASS F=0.3939393939393939;NALT=13;NREF=20;NOTH=0;CR=0.19502642784413077;RPR`
- *different allele/representation called nearby*
  - `chr20:367246 AC>A truth GT 0/1 [TR<=50] -> AAAC>A PASS F=0.45;NALT=9;NREF=11;NOTH=2;CR=0.11975493531830429;RPRI=0.05987747704; AAAC>A PASS F=0.14285714285714285;NALT=2;NREF=12;NOTH=9;CR=0.10932479580603646;RPR`
  - `chr20:367250 AAC>A truth GT 0/1 [TR<=50] -> AAAC>A PASS F=0.45;NALT=9;NREF=11;NOTH=2;CR=0.11975493531830429;RPRI=0.05987747704; AAAC>A PASS F=0.14285714285714285;NALT=2;NREF=12;NOTH=9;CR=0.10932479580603646;RPR`
  - `chr20:499958 CA>C truth GT 0/1 [HP>=7] -> CA>C PASS F=0.6190476190476191;NALT=13;NREF=8;NOTH=5;CR=0.4447384505149794;RPRI=; CAA>C PASS F=0.2;NALT=2;NREF=8;NOTH=16;CR=0.43781690864807266;RPRI=0.219965659073`
  - `chr20:653259 CG>C truth GT 0/1 [HP>=7] -> CG>C PASS F=0.39285714285714285;NALT=11;NREF=17;NOTH=2;CR=0.3150409482496533;RPR; CGG>C PASS F=0.10526315789473684;NALT=2;NREF=17;NOTH=11;CR=0.3150409482496533;RPR`
- *never a candidate (discovery miss)*
  - `chr20:527001 C>CAGACAGAGAGAGAGAG truth GT 0/1 [TR51-200] -> no candidate`
  - `chr20:1404738 CTTT>C truth GT 0/1 [HP>=7] -> no candidate`
  - `chr20:1849802 AAAAT>A truth GT 0/1 [TR<=50] -> no candidate`
  - `chr20:1914466 AG>A truth GT 1/1 [other] -> no candidate`
- *candidate found but filtered (LowProb/UNCERTAIN)*
  - `chr20:3101547 CTTTTTTTT>CTTT,C truth GT 2/1 [HP>=7] -> CTTTTTTTT>C UNCERTAIN F=1;NALT=2;NREF=0;NOTH=1;CR=0.4296594755469182;RPRI=0.6401336083518908`
  - `chr20:23554017 CA>CAA,C truth GT 2/1 [HP>=7] -> C>CA UNCERTAIN F=1;NALT=2;NREF=0;NOTH=0;CR=0.6742299682032735;RPRI=0.7624188546800684`
- *no truth variant nearby (pure FP) | HP>=7*
  - `chr20:92968 GCT>G called 0/1 truth ./. NALT=2 NREF=4 F=0.3333333333333333 HQ=40.75797720680305`
  - `chr20:444130 TAA>T called 0/1 truth ./. NALT=2 NREF=11 F=0.15384615384615385 HQ=11.092779901475069`
  - `chr20:736735 TC>T called 0/1 truth ./. NALT=2 NREF=3 F=0.4 HQ=10.677855668112215`
- *allele/representation mismatch near a truth variant | HP>=7*
  - `chr20:152467 AAGT>A called 0/1 truth ./. NALT=2 NREF=25 F=0.07407407407407407 HQ=44.995140076364805`
  - `chr20:365798 CTT>C called 1/1 truth ./. NALT=3 NREF=0 F=1 HQ=29.263237514086242`
  - `chr20:385413 TA>T called 0/1 truth ./. NALT=2 NREF=8 F=0.2 HQ=14.43916291405435`
- *genotype error (allele right, GT wrong) | TR51-200*
  - `chr20:217716 ATG>A called 1/1 truth 0/1 NALT=7 NREF=13 F=0.35 HQ=87.50731841679578`
  - `chr20:10507904 T>TTA called 1/1 truth 0/1 NALT=19 NREF=0 F=1 HQ=210.66927268749157`
  - `chr20:18488953 A>ATG called 1/1 truth 0/1 NALT=12 NREF=2 F=0.8571428571428571 HQ=215.619389882296`
- *allele/representation mismatch near a truth variant | TR<=50*
  - `chr20:232033 A>AAAATAAAT called 1/1 truth ./. NALT=6 NREF=0 F=1 HQ=271.9733577737402`
  - `chr20:232033 A>AAAATAAATAAAT called 1/1 truth ./. NALT=9 NREF=0 F=1 HQ=208.94344353179457`
  - `chr20:319668 CAA>C called 0/1 truth ./. NALT=2 NREF=30 F=0.0625 HQ=44.95238386176061`
- *no truth variant nearby (pure FP) | TR51-200*
  - `chr20:271443 TAC>T called 0/1 truth ./. NALT=3 NREF=14 F=0.17647058823529413 HQ=22.279940742246822`
  - `chr20:2011231 CCTCTCTCT>C called 1/1 truth ./. NALT=11 NREF=0 F=1 HQ=369.6604746303136`
  - `chr20:2011231 CCTCTCTCTCTCTCTCT>C called 1/1 truth ./. NALT=11 NREF=0 F=1 HQ=173.98602147346287`
- *no truth variant nearby (pure FP) | other*
  - `chr20:336556 A>ATTTT called 0/1 truth ./. NALT=2 NREF=23 F=0.08 HQ=13.115196667346368`
  - `chr20:525486 AC>A called 0/1 truth ./. NALT=2 NREF=14 F=0.125 HQ=9.851752203758963`
  - `chr20:676477 TCA>T called 0/1 truth ./. NALT=2 NREF=18 F=0.1 HQ=10.24619703151643`
- *no truth variant nearby (pure FP) | TR<=50*
  - `chr20:348484 CGT>C called 1/1 truth ./. NALT=6 NREF=0 F=1 HQ=81.70202796950718`
  - `chr20:348484 CGTGTGT>C called 0/1 truth ./. NALT=6 NREF=0 F=1 HQ=87.82039647301093`
  - `chr20:685696 AAATAATAAT>A called 1/1 truth ./. NALT=12 NREF=0 F=1 HQ=290.06510792984756`
- *genotype error (allele right, GT wrong) | HP>=7*
  - `chr20:454236 CAA>C called 1/1 truth 0/1 NALT=6 NREF=0 F=1 HQ=70.7052232338718`
  - `chr20:454295 G>GA called 1/1 truth 0/1 NALT=13 NREF=2 F=0.8666666666666667 HQ=99.20079102218102`
  - `chr20:887636 CT>C called 0/1 truth 1/1 NALT=3 NREF=1 F=0.75 HQ=28.588303504267913`
- *allele/representation mismatch near a truth variant | TR51-200*
  - `chr20:957300 ATG>A called 0/1 truth ./. NALT=2 NREF=25 F=0.07407407407407407 HQ=13.616187660561195`
  - `chr20:2945582 GGTGTGTGTGT>GGT,G called 2/1 truth ./. NALT=9 NREF=11 F=0.45 HQ=161.85597806068182`
  - `chr20:3310333 A>AAAAGAAAG called 1/1 truth ./. NALT=5 NREF=0 F=1 HQ=377.63050487562725`
- *allele/representation mismatch near a truth variant | other*
  - `chr20:1012670 A>AG called 0/1 truth ./. NALT=2 NREF=19 F=0.09523809523809523 HQ=23.478562747469084`
  - `chr20:2674873 AGT>A called 0/1 truth ./. NALT=2 NREF=13 F=0.13333333333333333 HQ=23.56192741170628`
  - `chr20:3079560 A>ATAT called 0/1 truth ./. NALT=2 NREF=22 F=0.08333333333333333 HQ=53.42672775584737`
- *genotype error (allele right, GT wrong) | TR<=50*
  - `chr20:4411737 T>TTG called 0/1 truth 1/1 NALT=15 NREF=3 F=0.8333333333333334 HQ=145.14078881986018`
  - `chr20:5454164 AATATAT>A called 1/1 truth 0/1 NALT=12 NREF=1 F=0.9230769230769231 HQ=294.1390946872243`
  - `chr20:9553557 A>AGT called 1/1 truth 0/1 NALT=16 NREF=11 F=0.5925925925925926 HQ=152.44434121580332`
- *genotype error (allele right, GT wrong) | other*
  - `chr20:19711132 A>AAC called 0/1 truth 1/1 NALT=17 NREF=1 F=0.9444444444444444 HQ=82.3275527086239`
  - `chr20:22513143 TCTCTCTCC>T called 0/1 truth 1/1 NALT=9 NREF=13 F=0.4090909090909091 HQ=172.09093494068287`
  - `chr20:37236093 T>TTGTGGGTGTGTG called 0/1 truth 1/1 NALT=12 NREF=1 F=0.9230769230769231 HQ=199.81644325080117`

## Error taxonomy: nindels_67x


**False negatives by cause**

| Cause | Count |
|---|---|
| multi-allelic site (truth 1/2), alleles genotyped separately | 481 |
| different allele/representation called nearby | 242 |
| genotype wrong (allele found, GT differs) | 69 |
| never a candidate (discovery miss) | 31 |

**False positives by cause and repeat context**

| Cause \| context | Count |
|---|---|
| allele/representation mismatch near a truth variant \| HP>=7 | 1256 |
| no truth variant nearby (pure FP) \| HP>=7 | 689 |
| allele/representation mismatch near a truth variant \| TR<=50 | 361 |
| no truth variant nearby (pure FP) \| TR<=50 | 237 |
| allele/representation mismatch near a truth variant \| other | 157 |
| allele/representation mismatch near a truth variant \| TR51-200 | 152 |
| no truth variant nearby (pure FP) \| TR51-200 | 93 |
| genotype error (allele right, GT wrong) \| HP>=7 | 58 |
| no truth variant nearby (pure FP) \| other | 56 |
| genotype error (allele right, GT wrong) \| TR51-200 | 10 |
| genotype error (allele right, GT wrong) \| TR<=50 | 7 |
| genotype error (allele right, GT wrong) \| other | 2 |

**FP called genotype**

- allele/representation mismatch near a truth variant, called 0/1: 1164
- no truth variant nearby, called 0/1: 588
- allele/representation mismatch near a truth variant, called 1/1: 553
- no truth variant nearby, called 1/1: 463
- allele/representation mismatch near a truth variant, called 2/1: 209
- genotype error, called 0/1: 54
- no truth variant nearby, called 2/1: 24
- genotype error, called 1/1: 23

**Recall: clustered (another truth variant within 10 bp) vs isolated**

- isolated: 9632/10352 = 0.9304
- clustered: 545/648 = 0.8410

**Median features: TP calls vs pure-FP calls**

| Feature | TP | pure FP |
|---|---|---|
| NALT | 32.00 | 5.00 |
| NREF | 23.00 | 1.00 |
| F | 0.55 | 0.95 |
| HQ | 317.28 | 146.82 |
| QMIN | 0.91 | 0.79 |

Pure FPs with NALT <= 2: 337/1075

**Examples**

- *genotype wrong (allele found, GT differs)*
  - `chr20:217716 ATG>A truth GT 0/1 [TR51-200] -> ATG>A PASS F=0.42857142857142855;NALT=18;NREF=24;NOTH=0;CR=0.30533254900671436;RP`
  - `chr20:887636 CT>C truth GT 1/1 [HP>=7] -> CT>C PASS F=0.75;NALT=6;NREF=2;NOTH=8;CR=0.46569580501495844;RPRI=0.232850381680`
  - `chr20:1138518 C>CA truth GT 1/1 [HP>=7] -> C>CA PASS F=0.75;NALT=15;NREF=5;NOTH=4;CR=0.5118917596697001;RPRI=0.255945879834`
  - `chr20:2542038 G>GA truth GT 1/1 [HP>=7] -> G>GA PASS F=0.9818181818181818;NALT=54;NREF=1;NOTH=7;CR=0.2216611847967109;RPRI=; G>GAA LowProb F=0.6666666666666666;NALT=2;NREF=1;NOTH=59;CR=0.2216611847967109;RPRI=`
- *multi-allelic site (truth 1/2), alleles genotyped separately*
  - `chr20:232033 A>AAAATAAAT,AAAATAAATAAAT truth GT 1/2 [TR<=50] -> A>AAAATAAAT PASS F=1;NALT=14;NREF=0;NOTH=22;CR=0.29503964290122536;RPRI=0.1475198214506; A>AAAATAAATAAAT PASS F=1;NALT=19;NREF=0;NOTH=17;CR=0.29503964290122536;RPRI=0.1475198214506`
  - `chr20:348484 CGTGTGT>CGTGT,C truth GT 2/1 [TR<=50] -> CGT>C PASS F=1;NALT=11;NREF=0;NOTH=14;CR=0.3275295386120749;RPRI=0.16376476930603; CGTGTGT>C PASS F=1;NALT=14;NREF=0;NOTH=11;CR=0.3275295386120749;RPRI=0.16376476930603`
  - `chr20:365798 CTT>CT,C truth GT 2/1 [HP>=7] -> CT>C PASS F=1;NALT=5;NREF=0;NOTH=7;CR=0.5293044287932317;RPRI=0.2646541877717401; CTT>C PASS F=1;NALT=6;NREF=0;NOTH=6;CR=0.5293044287932317;RPRI=0.2646522349573023`
  - `chr20:511703 TACAC>CACAC,T truth GT 2/1 [TR<=50] -> TACAC>T PASS F=0.5142857142857142;NALT=36;NREF=34;NOTH=1;CR=0.2010016732559959;RPRI`
- *different allele/representation called nearby*
  - `chr20:316724 GTT>G truth GT 0/1 [HP>=7] -> GT>G PASS F=0.08;NALT=2;NREF=23;NOTH=15;CR=0.4587943148124999;RPRI=0.22939719495; GTT>G PASS F=0.3783783783783784;NALT=14;NREF=23;NOTH=3;CR=0.45264214230694844;RPR`
  - `chr20:367246 AC>A truth GT 0/1 [TR<=50] -> AAAC>A PASS F=0.48936170212765956;NALT=23;NREF=24;NOTH=3;CR=0.11339129895466793;RP; AAAC>A PASS F=0.10344827586206896;NALT=3;NREF=26;NOTH=23;CR=0.10326292289633747;RP`
  - `chr20:367250 AAC>A truth GT 0/1 [TR<=50] -> AAAC>A PASS F=0.48936170212765956;NALT=23;NREF=24;NOTH=3;CR=0.11339129895466793;RP; AAAC>A PASS F=0.10344827586206896;NALT=3;NREF=26;NOTH=23;CR=0.10326292289633747;RP`
  - `chr20:499958 CA>C truth GT 0/1 [HP>=7] -> AG>A LowProb F=0.05263157894736842;NALT=2;NREF=36;NOTH=36;CR=0.13566070966823549;RP; GC>G LowProb F=0.05405405405405406;NALT=2;NREF=35;NOTH=36;CR=0.13760917518252758;RP`
- *never a candidate (discovery miss)*
  - `chr20:1404738 CTTT>C truth GT 0/1 [HP>=7] -> no candidate`
  - `chr20:1914466 AG>A truth GT 1/1 [other] -> no candidate`
  - `chr20:1914707 CT>C truth GT 1/1 [other] -> no candidate`
  - `chr20:1915067 T>TAA truth GT 1/1 [other] -> no candidate`
- *no truth variant nearby (pure FP) | HP>=7*
  - `chr20:92967 TGC>T called 0/1 truth ./. NALT=2 NREF=10 F=0.16666666666666666 HQ=28.763800251348016`
  - `chr20:92968 GCT>G called 0/1 truth ./. NALT=3 NREF=9 F=0.25 HQ=64.91670796801418`
  - `chr20:152447 TAGC>T called 0/1 truth ./. NALT=3 NREF=14 F=0.17647058823529413 HQ=33.35242703083259`
- *allele/representation mismatch near a truth variant | HP>=7*
  - `chr20:151662 CT>C called 0/1 truth ./. NALT=2 NREF=1 F=0.6666666666666666 HQ=9.68462166253721`
  - `chr20:240075 C>CAA called 1/1 truth ./. NALT=4 NREF=0 F=1 HQ=143.2504682826077`
  - `chr20:314924 TG>T called 0/1 truth ./. NALT=3 NREF=21 F=0.125 HQ=27.8069860433784`
- *genotype error (allele right, GT wrong) | TR51-200*
  - `chr20:217716 ATG>A called 1/1 truth 0/1 NALT=18 NREF=24 F=0.42857142857142855 HQ=194.9873812657304`
  - `chr20:10507904 T>TTA called 1/1 truth 0/1 NALT=35 NREF=1 F=0.9722222222222222 HQ=395.16450002071497`
  - `chr20:16801236 T>TTG called 0/1 truth 1/1 NALT=10 NREF=4 F=0.7142857142857143 HQ=120.92990948263673`
- *allele/representation mismatch near a truth variant | TR<=50*
  - `chr20:232033 A>AAAATAAAT called 1/1 truth ./. NALT=14 NREF=0 F=1 HQ=577.3326278784127`
  - `chr20:232033 A>AAAATAAATAAAT called 1/1 truth ./. NALT=19 NREF=0 F=1 HQ=450.9537169040026`
  - `chr20:319668 CAA>C called 0/1 truth ./. NALT=2 NREF=51 F=0.03773584905660377 HQ=51.32600662086111`
- *no truth variant nearby (pure FP) | TR<=50*
  - `chr20:348484 CGT>C called 1/1 truth ./. NALT=11 NREF=0 F=1 HQ=163.6232470038858`
  - `chr20:348484 CGTGTGT>C called 0/1 truth ./. NALT=14 NREF=0 F=1 HQ=207.1052065032933`
  - `chr20:685696 AAATAATAAT>A called 1/1 truth ./. NALT=27 NREF=0 F=1 HQ=837.9012780926106`
- *allele/representation mismatch near a truth variant | other*
  - `chr20:611880 GA>G called 0/1 truth ./. NALT=2 NREF=60 F=0.03225806451612903 HQ=11.926755947063839`
  - `chr20:1012670 A>AG called 0/1 truth ./. NALT=2 NREF=45 F=0.0425531914893617 HQ=16.45180216478252`
  - `chr20:1362011 TG>T called 0/1 truth ./. NALT=2 NREF=55 F=0.03508771929824561 HQ=17.579218780969086`
- *genotype error (allele right, GT wrong) | HP>=7*
  - `chr20:887636 CT>C called 0/1 truth 1/1 NALT=6 NREF=2 F=0.75 HQ=72.86555990533355`
  - `chr20:1138518 C>CA called 0/1 truth 1/1 NALT=15 NREF=5 F=0.75 HQ=118.8957190297765`
  - `chr20:2542038 G>GA called 0/1 truth 1/1 NALT=54 NREF=1 F=0.9818181818181818 HQ=286.26300565390375`
- *no truth variant nearby (pure FP) | TR51-200*
  - `chr20:2011231 CCTCTCTCT>C called 1/1 truth ./. NALT=25 NREF=0 F=1 HQ=915.7641071938913`
  - `chr20:2011231 CCTCTCTCTCTCTCTCT>C called 1/1 truth ./. NALT=31 NREF=0 F=1 HQ=498.54588307890367`
  - `chr20:4779162 TTGTGCG>T called 0/1 truth ./. NALT=2 NREF=13 F=0.13333333333333333 HQ=301.1240480177976`
- *allele/representation mismatch near a truth variant | TR51-200*
  - `chr20:2249002 CTT>C called 0/1 truth ./. NALT=2 NREF=30 F=0.0625 HQ=124.22497889371562`
  - `chr20:2249004 T>TTC called 0/1 truth ./. NALT=2 NREF=31 F=0.06060606060606061 HQ=14.169981185321376`
  - `chr20:2945582 GGTGTGTGTGT>GGT,G called 2/1 truth ./. NALT=24 NREF=24 F=0.5 HQ=383.5165737073278`
- *no truth variant nearby (pure FP) | other*
  - `chr20:4374492 CA>C called 0/1 truth ./. NALT=7 NREF=42 F=0.14285714285714285 HQ=44.54178567920849`
  - `chr20:4374503 CTACCTCCAACCTT>C called 0/1 truth ./. NALT=3 NREF=36 F=0.07692307692307693 HQ=36.74718330053361`
  - `chr20:7508395 G>GA called 0/1 truth ./. NALT=13 NREF=66 F=0.16455696202531644 HQ=55.96247500393963`
- *genotype error (allele right, GT wrong) | TR<=50*
  - `chr20:5454164 AATATAT>A called 1/1 truth 0/1 NALT=23 NREF=1 F=0.9583333333333334 HQ=631.3292268302607`
  - `chr20:10476142 C>CTG called 1/1 truth 0/1 NALT=29 NREF=0 F=1 HQ=348.1402439049523`
  - `chr20:10476143 T>TGC called 1/1 truth 0/1 NALT=22 NREF=5 F=0.8148148148148148 HQ=277.1129285109305`
- *genotype error (allele right, GT wrong) | other*
  - `chr20:22513143 TCTCTCTCC>T called 0/1 truth 1/1 NALT=13 NREF=24 F=0.35135135135135137 HQ=293.55138102595015`
  - `chr20:63701286 A>AC called 1/1 truth 0/1 NALT=35 NREF=41 F=0.4605263157894737 HQ=304.60559228727266`
