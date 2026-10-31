## Benchmark: HG004 chr20, GIAB v4.2.1 confident regions, hap.py (xcmp)

| Caller | Wall (s) | CPU (s) | Peak RSS (MB) | SNP F1 | INDEL F1 | INDEL Precision | INDEL Recall | INDEL TP | INDEL FP | INDEL FN | GT concordance |
|---|---|---|---|---|---|---|---|---|---|---|---|
| nindels_30x | 29.2 (n=3) | 160.7 | 502 | N/A (no SNP calls) | 0.8743 | 0.8334 | 0.9194 | 10113 | 2023 | 887 | 0.9879 (124 GT err) |
| expA_multi_30x | 28.9 (n=3) | 185.0 | 549 | N/A (no SNP calls) | 0.9438 | 0.9412 | 0.9464 | 10410 | 674 | 590 | 0.9863 (145 GT err) |
| expB_str_30x | 36.3 (n=3) | 236.7 | 498 | N/A (no SNP calls) | 0.8871 | 0.8646 | 0.9108 | 10019 | 1573 | 981 | 0.9927 (75 GT err) |
| expC_frag_30x | 18.6 (n=3) | 125.4 | 512 | N/A (no SNP calls) | 0.8438 | 0.8794 | 0.8110 | 8921 | 1227 | 2079 | 0.9825 (175 GT err) |
| expC2_fragend_30x | 21.9 (n=3) | 130.1 | 543 | N/A (no SNP calls) | 0.8455 | 0.8727 | 0.8199 | 9019 | 1319 | 1981 | 0.9829 (171 GT err) |
| expAB_30x | 39.1 (n=3) | 216.0 | 519 | N/A (no SNP calls) | 0.9397 | 0.9515 | 0.9282 | 10210 | 539 | 790 | 0.9913 (92 GT err) |
| expABC2_30x | 25.5 (n=3) | 138.3 | 538 | N/A (no SNP calls) | 0.8697 | 0.9566 | 0.7973 | 8770 | 410 | 2230 | 0.9903 (97 GT err) |
| bcftools_30x | 660.5 (n=1) | 581.0 | 110 | 0.9911 | 0.9330 | 0.9483 | 0.9181 | 10099 | 578 | 901 | 0.9747 (262 GT err) |
| nindels_67x | 64.9 (n=3) | 367.5 | 960 | N/A (no SNP calls) | 0.8390 | 0.7678 | 0.9248 | 10173 | 3078 | 827 | 0.9925 (77 GT err) |
| expA_multi_67x | 65.6 (n=1) | 417.5 | 1044 | N/A (no SNP calls) | 0.9516 | 0.9445 | 0.9589 | 10548 | 644 | 452 | 0.9870 (139 GT err) |
| expB_str_67x | 115.3 (n=1) | 538.0 | 893 | N/A (no SNP calls) | 0.8697 | 0.8212 | 0.9244 | 10168 | 2217 | 832 | 0.9954 (47 GT err) |
| expC_frag_67x | 65.7 (n=1) | 325.2 | 929 | N/A (no SNP calls) | 0.8575 | 0.8227 | 0.8954 | 9849 | 2124 | 1151 | 0.9899 (103 GT err) |
| expC2_fragend_67x | 63.7 (n=1) | 336.6 | 948 | N/A (no SNP calls) | 0.8521 | 0.8087 | 0.9003 | 9903 | 2344 | 1097 | 0.9921 (80 GT err) |
| expAB_67x | 70.5 (n=1) | 490.9 | 1000 | N/A (no SNP calls) | 0.9548 | 0.9504 | 0.9594 | 10553 | 572 | 447 | 0.9917 (89 GT err) |
| expABC2_67x | 67.8 (n=1) | 359.1 | 1063 | N/A (no SNP calls) | 0.9394 | 0.9548 | 0.9244 | 10168 | 498 | 832 | 0.9931 (72 GT err) |
| bcftools_67x | 1019.8 (n=1) | 1019.2 | 111 | 0.9926 | 0.9548 | 0.9627 | 0.9470 | 10417 | 426 | 583 | 0.9770 (245 GT err) |

## INDEL (PASS) by length class, genotype and repeat context

| Caller | Stratum | Truth | Recall | Precision | F1 |
|---|---|---|---|---|---|
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
| expA_multi_30x | homopolymer_ge7 | 4908 | 0.9159 | 0.9136 | 0.9147 |
| expA_multi_30x | not_in_TR_or_HP | 3876 | 0.9894 | 0.9841 | 0.9867 |
| expA_multi_30x | tandem_repeat_51to200 | 653 | 0.9479 | 0.9430 | 0.9454 |
| expA_multi_30x | tandem_repeat_le50 | 1449 | 0.9337 | 0.9227 | 0.9282 |
| expA_multi_30x | D16_PLUS | 136 | 0.9485 | 0.9556 | 0.9520 |
| expA_multi_30x | D1_5 | 5066 | 0.9406 | 0.9277 | 0.9341 |
| expA_multi_30x | D6_15 | 610 | 0.9459 | 0.9263 | 0.9360 |
| expA_multi_30x | I16_PLUS | 104 | 0.9712 | 0.9712 | 0.9712 |
| expA_multi_30x | I1_5 | 4941 | 0.9429 | 0.9559 | 0.9493 |
| expA_multi_30x | I6_15 | 557 | 0.9354 | 0.9507 | 0.9430 |
| expB_str_30x | homopolymer_ge7 | 4908 | 0.8771 | 0.8394 | 0.8578 |
| expB_str_30x | not_in_TR_or_HP | 3876 | 0.9912 | 0.9675 | 0.9792 |
| expB_str_30x | tandem_repeat_51to200 | 653 | 0.8407 | 0.7565 | 0.7964 |
| expB_str_30x | tandem_repeat_le50 | 1449 | 0.8385 | 0.7422 | 0.7874 |
| expB_str_30x | D16_PLUS | 136 | 0.8529 | 0.8593 | 0.8561 |
| expB_str_30x | D1_5 | 5066 | 0.9084 | 0.8625 | 0.8849 |
| expB_str_30x | D6_15 | 610 | 0.8197 | 0.7504 | 0.7835 |
| expB_str_30x | I16_PLUS | 104 | 0.8846 | 0.8679 | 0.8762 |
| expB_str_30x | I1_5 | 4941 | 0.9065 | 0.9033 | 0.9049 |
| expB_str_30x | I6_15 | 557 | 0.7684 | 0.7179 | 0.7423 |
| expC_frag_30x | homopolymer_ge7 | 4908 | 0.7209 | 0.8497 | 0.7800 |
| expC_frag_30x | not_in_TR_or_HP | 3876 | 0.9494 | 0.9800 | 0.9645 |
| expC_frag_30x | tandem_repeat_51to200 | 653 | 0.7703 | 0.7788 | 0.7745 |
| expC_frag_30x | tandem_repeat_le50 | 1449 | 0.7619 | 0.7510 | 0.7564 |
| expC_frag_30x | D16_PLUS | 136 | 0.7353 | 0.8860 | 0.8036 |
| expC_frag_30x | D1_5 | 5066 | 0.8032 | 0.8728 | 0.8366 |
| expC_frag_30x | D6_15 | 610 | 0.7295 | 0.7884 | 0.7578 |
| expC_frag_30x | I16_PLUS | 104 | 0.7981 | 0.8646 | 0.8300 |
| expC_frag_30x | I1_5 | 4941 | 0.8061 | 0.9160 | 0.8576 |
| expC_frag_30x | I6_15 | 557 | 0.7163 | 0.7575 | 0.7364 |
| expC2_fragend_30x | homopolymer_ge7 | 4908 | 0.7492 | 0.8381 | 0.7912 |
| expC2_fragend_30x | not_in_TR_or_HP | 3876 | 0.9407 | 0.9772 | 0.9586 |
| expC2_fragend_30x | tandem_repeat_51to200 | 653 | 0.7611 | 0.7743 | 0.7677 |
| expC2_fragend_30x | tandem_repeat_le50 | 1449 | 0.7585 | 0.7517 | 0.7551 |
| expC2_fragend_30x | D16_PLUS | 136 | 0.6985 | 0.8649 | 0.7728 |
| expC2_fragend_30x | D1_5 | 5066 | 0.8139 | 0.8588 | 0.8357 |
| expC2_fragend_30x | D6_15 | 610 | 0.7311 | 0.7915 | 0.7601 |
| expC2_fragend_30x | I16_PLUS | 104 | 0.7981 | 0.8842 | 0.8389 |
| expC2_fragend_30x | I1_5 | 4941 | 0.8179 | 0.9152 | 0.8638 |
| expC2_fragend_30x | I6_15 | 557 | 0.7092 | 0.7590 | 0.7332 |
| expAB_30x | homopolymer_ge7 | 4908 | 0.8780 | 0.9308 | 0.9036 |
| expAB_30x | not_in_TR_or_HP | 3876 | 0.9894 | 0.9846 | 0.9870 |
| expAB_30x | tandem_repeat_51to200 | 653 | 0.9342 | 0.9476 | 0.9408 |
| expAB_30x | tandem_repeat_le50 | 1449 | 0.9296 | 0.9340 | 0.9318 |
| expAB_30x | D16_PLUS | 136 | 0.9485 | 0.9627 | 0.9556 |
| expAB_30x | D1_5 | 5066 | 0.9187 | 0.9382 | 0.9283 |
| expAB_30x | D6_15 | 610 | 0.9508 | 0.9370 | 0.9439 |
| expAB_30x | I16_PLUS | 104 | 0.9712 | 0.9712 | 0.9712 |
| expAB_30x | I1_5 | 4941 | 0.9235 | 0.9665 | 0.9445 |
| expAB_30x | I6_15 | 557 | 0.9354 | 0.9578 | 0.9465 |
| expABC2_30x | homopolymer_ge7 | 4908 | 0.6948 | 0.9428 | 0.8000 |
| expABC2_30x | not_in_TR_or_HP | 3876 | 0.9340 | 0.9829 | 0.9578 |
| expABC2_30x | tandem_repeat_51to200 | 653 | 0.7810 | 0.9394 | 0.8529 |
| expABC2_30x | tandem_repeat_le50 | 1449 | 0.7805 | 0.9298 | 0.8487 |
| expABC2_30x | D16_PLUS | 136 | 0.7647 | 0.9541 | 0.8490 |
| expABC2_30x | D1_5 | 5066 | 0.7837 | 0.9485 | 0.8582 |
| expABC2_30x | D6_15 | 610 | 0.8131 | 0.9424 | 0.8730 |
| expABC2_30x | I16_PLUS | 104 | 0.8558 | 0.9468 | 0.8990 |
| expABC2_30x | I1_5 | 4941 | 0.7911 | 0.9685 | 0.8709 |
| expABC2_30x | I6_15 | 557 | 0.8276 | 0.9526 | 0.8857 |
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
| expA_multi_67x | homopolymer_ge7 | 4908 | 0.9399 | 0.9216 | 0.9307 |
| expA_multi_67x | not_in_TR_or_HP | 3876 | 0.9879 | 0.9810 | 0.9844 |
| expA_multi_67x | tandem_repeat_51to200 | 653 | 0.9541 | 0.9394 | 0.9467 |
| expA_multi_67x | tandem_repeat_le50 | 1449 | 0.9462 | 0.9315 | 0.9388 |
| expA_multi_67x | D16_PLUS | 136 | 0.9779 | 0.9638 | 0.9708 |
| expA_multi_67x | D1_5 | 5066 | 0.9510 | 0.9301 | 0.9405 |
| expA_multi_67x | D6_15 | 610 | 0.9525 | 0.9326 | 0.9424 |
| expA_multi_67x | I16_PLUS | 104 | 0.9904 | 0.9717 | 0.9810 |
| expA_multi_67x | I1_5 | 4941 | 0.9587 | 0.9593 | 0.9590 |
| expA_multi_67x | I6_15 | 557 | 0.9641 | 0.9555 | 0.9598 |
| expB_str_67x | homopolymer_ge7 | 4908 | 0.9085 | 0.7826 | 0.8408 |
| expB_str_67x | not_in_TR_or_HP | 3876 | 0.9910 | 0.9559 | 0.9731 |
| expB_str_67x | tandem_repeat_51to200 | 653 | 0.8499 | 0.7025 | 0.7692 |
| expB_str_67x | tandem_repeat_le50 | 1449 | 0.8302 | 0.6918 | 0.7547 |
| expB_str_67x | D16_PLUS | 136 | 0.8750 | 0.8500 | 0.8623 |
| expB_str_67x | D1_5 | 5066 | 0.9287 | 0.8147 | 0.8680 |
| expB_str_67x | D6_15 | 610 | 0.8082 | 0.6852 | 0.7416 |
| expB_str_67x | I16_PLUS | 104 | 0.8846 | 0.8070 | 0.8440 |
| expB_str_67x | I1_5 | 4941 | 0.9233 | 0.8687 | 0.8952 |
| expB_str_67x | I6_15 | 557 | 0.7630 | 0.6630 | 0.7095 |
| expC_frag_67x | homopolymer_ge7 | 4908 | 0.8517 | 0.7644 | 0.8057 |
| expC_frag_67x | not_in_TR_or_HP | 3876 | 0.9837 | 0.9739 | 0.9788 |
| expC_frag_67x | tandem_repeat_51to200 | 653 | 0.8346 | 0.7290 | 0.7782 |
| expC_frag_67x | tandem_repeat_le50 | 1449 | 0.8316 | 0.7022 | 0.7614 |
| expC_frag_67x | D16_PLUS | 136 | 0.7941 | 0.8438 | 0.8182 |
| expC_frag_67x | D1_5 | 5066 | 0.8968 | 0.7930 | 0.8417 |
| expC_frag_67x | D6_15 | 610 | 0.7967 | 0.7424 | 0.7686 |
| expC_frag_67x | I16_PLUS | 104 | 0.8558 | 0.8241 | 0.8396 |
| expC_frag_67x | I1_5 | 4941 | 0.8935 | 0.8847 | 0.8891 |
| expC_frag_67x | I6_15 | 557 | 0.7540 | 0.6988 | 0.7254 |
| expC2_fragend_67x | homopolymer_ge7 | 4908 | 0.8643 | 0.7433 | 0.7992 |
| expC2_fragend_67x | not_in_TR_or_HP | 3876 | 0.9830 | 0.9709 | 0.9769 |
| expC2_fragend_67x | tandem_repeat_51to200 | 653 | 0.8331 | 0.7190 | 0.7718 |
| expC2_fragend_67x | tandem_repeat_le50 | 1449 | 0.8288 | 0.6956 | 0.7564 |
| expC2_fragend_67x | D16_PLUS | 136 | 0.8088 | 0.8462 | 0.8271 |
| expC2_fragend_67x | D1_5 | 5066 | 0.9031 | 0.7685 | 0.8304 |
| expC2_fragend_67x | D6_15 | 610 | 0.7967 | 0.7360 | 0.7652 |
| expC2_fragend_67x | I16_PLUS | 104 | 0.8558 | 0.8241 | 0.8396 |
| expC2_fragend_67x | I1_5 | 4941 | 0.8988 | 0.8824 | 0.8906 |
| expC2_fragend_67x | I6_15 | 557 | 0.7522 | 0.6960 | 0.7230 |
| expAB_67x | homopolymer_ge7 | 4908 | 0.9415 | 0.9341 | 0.9378 |
| expAB_67x | not_in_TR_or_HP | 3876 | 0.9879 | 0.9805 | 0.9842 |
| expAB_67x | tandem_repeat_51to200 | 653 | 0.9525 | 0.9380 | 0.9452 |
| expAB_67x | tandem_repeat_le50 | 1449 | 0.9448 | 0.9331 | 0.9389 |
| expAB_67x | D16_PLUS | 136 | 0.9779 | 0.9638 | 0.9708 |
| expAB_67x | D1_5 | 5066 | 0.9489 | 0.9355 | 0.9421 |
| expAB_67x | D6_15 | 610 | 0.9541 | 0.9327 | 0.9433 |
| expAB_67x | I16_PLUS | 104 | 0.9904 | 0.9717 | 0.9810 |
| expAB_67x | I1_5 | 4941 | 0.9613 | 0.9668 | 0.9641 |
| expAB_67x | I6_15 | 557 | 0.9659 | 0.9624 | 0.9642 |
| expABC2_67x | homopolymer_ge7 | 4908 | 0.8810 | 0.9404 | 0.9097 |
| expABC2_67x | not_in_TR_or_HP | 3876 | 0.9801 | 0.9847 | 0.9824 |
| expABC2_67x | tandem_repeat_51to200 | 653 | 0.9173 | 0.9371 | 0.9271 |
| expABC2_67x | tandem_repeat_le50 | 1449 | 0.9248 | 0.9346 | 0.9297 |
| expABC2_67x | D16_PLUS | 136 | 0.9118 | 0.9612 | 0.9358 |
| expABC2_67x | D1_5 | 5066 | 0.9177 | 0.9425 | 0.9299 |
| expABC2_67x | D6_15 | 610 | 0.9180 | 0.9410 | 0.9294 |
| expABC2_67x | I16_PLUS | 104 | 0.9615 | 0.9524 | 0.9569 |
| expABC2_67x | I1_5 | 4941 | 0.9205 | 0.9691 | 0.9441 |
| expABC2_67x | I6_15 | 557 | 0.9300 | 0.9647 | 0.9470 |
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

## Error taxonomy: expA_multi_30x


**False negatives by cause**

| Cause | Count |
|---|---|
| different allele/representation called nearby | 165 |
| multi-allelic site (truth 1/2), alleles genotyped separately | 148 |
| genotype wrong (allele found, GT differs) | 130 |
| never a candidate (discovery miss) | 116 |
| candidate found but filtered (LowProb/UNCERTAIN) | 1 |

**False positives by cause and repeat context**

| Cause \| context | Count |
|---|---|
| allele/representation mismatch near a truth variant \| HP>=7 | 288 |
| genotype error (allele right, GT wrong) \| HP>=7 | 129 |
| allele/representation mismatch near a truth variant \| TR<=50 | 111 |
| no truth variant nearby (pure FP) \| other | 38 |
| allele/representation mismatch near a truth variant \| TR51-200 | 35 |
| allele/representation mismatch near a truth variant \| other | 29 |
| no truth variant nearby (pure FP) \| HP>=7 | 24 |
| genotype error (allele right, GT wrong) \| TR<=50 | 9 |
| genotype error (allele right, GT wrong) \| other | 4 |
| genotype error (allele right, GT wrong) \| TR51-200 | 3 |
| no truth variant nearby (pure FP) \| TR51-200 | 3 |
| no truth variant nearby (pure FP) \| TR<=50 | 1 |

**FP called genotype**

- allele/representation mismatch near a truth variant, called 0/1: 373
- genotype error, called 0/1: 127
- no truth variant nearby, called 0/1: 66
- allele/representation mismatch near a truth variant, called 1/1: 55
- allele/representation mismatch near a truth variant, called 2/1: 35
- genotype error, called 1/1: 18

**Recall: clustered (another truth variant within 10 bp) vs isolated**

- isolated: 9922/10352 = 0.9585
- clustered: 518/648 = 0.7994

**Median features: TP calls vs pure-FP calls**

| Feature | TP | pure FP |
|---|---|---|
| NALT | n/a | n/a |
| NREF | n/a | n/a |
| F | n/a | n/a |
| HQ | n/a | n/a |
| QMIN | n/a | n/a |

**Examples**

- *genotype wrong (allele found, GT differs)*
  - `chr20:217716 ATG>A truth GT 0/1 [TR51-200] -> ATG>A PASS K=1;DP=20`
  - `chr20:454236 CAA>C truth GT 0/1 [HP>=7] -> CAA>C PASS K=1;DP=8`
  - `chr20:454295 G>GA truth GT 0/1 [HP>=7] -> G>GA PASS K=1;DP=15`
  - `chr20:1122295 CTTTT>C truth GT 1/1 [HP>=7] -> CTTTT>C PASS K=1;DP=16`
- *multi-allelic site (truth 1/2), alleles genotyped separately*
  - `chr20:365798 CTT>CT,C truth GT 2/1 [HP>=7] -> CTT>C PASS K=1;DP=4`
  - `chr20:511703 TACAC>CACAC,T truth GT 2/1 [TR<=50] -> TACAC>T PASS K=1;DP=31`
  - `chr20:586045 CT>CTT,C truth GT 2/1 [HP>=7] -> CT>C PASS K=1;DP=9`
  - `chr20:942277 A>AT,T truth GT 1/2 [TR<=50] -> A>AT PASS K=2;DP=26`
- *different allele/representation called nearby*
  - `chr20:367246 AC>A truth GT 0/1 [TR<=50] -> AAAC>A PASS K=2;DP=23`
  - `chr20:367250 AAC>A truth GT 0/1 [TR<=50] -> AAAC>A PASS K=2;DP=23`
  - `chr20:499958 CA>C truth GT 0/1 [HP>=7] -> CAA>C PASS K=2;DP=26`
  - `chr20:887636 CT>C truth GT 1/1 [HP>=7] -> CT>C PASS K=1;DP=7`
- *never a candidate (discovery miss)*
  - `chr20:527001 C>CAGACAGAGAGAGAGAG truth GT 0/1 [TR51-200] -> no candidate`
  - `chr20:1404738 CTTT>C truth GT 0/1 [HP>=7] -> no candidate`
  - `chr20:1849802 AAAAT>A truth GT 0/1 [TR<=50] -> no candidate`
  - `chr20:1914466 AG>A truth GT 1/1 [other] -> no candidate`
- *candidate found but filtered (LowProb/UNCERTAIN)*
  - `chr20:13366867 T>TACAC truth GT 0/1 [TR<=50] -> C>CACA LowProb K=1;DP=15`
- *genotype error (allele right, GT wrong) | TR51-200*
  - `chr20:217716 ATG>A called 1/1 truth 0/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:21798505 C>CTG called 0/1 truth 1/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:35321252 C>CTG called 1/1 truth 0/1 NALT=None NREF=None F=None HQ=None`
- *no truth variant nearby (pure FP) | TR51-200*
  - `chr20:271443 TAC>T called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:8490860 GCATATATATATATATATGTACACATATCTATATGTACACATATATATGTGTA>G called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:40465791 TTG>T called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *no truth variant nearby (pure FP) | other*
  - `chr20:336556 A>ATTTT called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:4374492 CA>C called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:4374503 CTACCTCCAACCTT>C called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *allele/representation mismatch near a truth variant | HP>=7*
  - `chr20:365798 CTT>C called 1/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:458330 TA>T called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:499958 CAA>C called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *allele/representation mismatch near a truth variant | TR<=50*
  - `chr20:367249 AAAC>A called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:511703 TACAC>T called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:942277 A>AT called 1/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *genotype error (allele right, GT wrong) | HP>=7*
  - `chr20:454236 CAA>C called 1/1 truth 0/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:454295 G>GA called 1/1 truth 0/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:887636 CT>C called 0/1 truth 1/1 NALT=None NREF=None F=None HQ=None`
- *no truth variant nearby (pure FP) | HP>=7*
  - `chr20:1700552 TA>T called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:4208933 TG>T called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:5972193 CT>C called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *allele/representation mismatch near a truth variant | TR51-200*
  - `chr20:3982103 GTC>G called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:3982105 C>CTCTCTCTCTC called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:6124250 G>GCA,GCACA called 2/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *genotype error (allele right, GT wrong) | TR<=50*
  - `chr20:4411737 T>TTG called 0/1 truth 1/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:9553557 A>AGT called 1/1 truth 0/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:35530219 A>AATT called 0/1 truth 1/1 NALT=None NREF=None F=None HQ=None`
- *allele/representation mismatch near a truth variant | other*
  - `chr20:4619359 TC>T called 1/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:9051911 T>TA called 1/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:13929964 TA>T called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *genotype error (allele right, GT wrong) | other*
  - `chr20:19711132 A>AAC called 0/1 truth 1/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:22513143 TCTCTCTCC>T called 0/1 truth 1/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:37236093 T>TTGTGGGTGTGTG called 0/1 truth 1/1 NALT=None NREF=None F=None HQ=None`
- *no truth variant nearby (pure FP) | TR<=50*
  - `chr20:26079727 ACTAATCACTGGCTGCAACT>A called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`

## Error taxonomy: expB_str_30x


**False negatives by cause**

| Cause | Count |
|---|---|
| multi-allelic site (truth 1/2), alleles genotyped separately | 496 |
| different allele/representation called nearby | 116 |
| never a candidate (discovery miss) | 105 |
| genotype wrong (allele found, GT differs) | 73 |
| candidate found but filtered (LowProb/UNCERTAIN) | 7 |

**False positives by cause and repeat context**

| Cause \| context | Count |
|---|---|
| allele/representation mismatch near a truth variant \| HP>=7 | 478 |
| no truth variant nearby (pure FP) \| HP>=7 | 296 |
| allele/representation mismatch near a truth variant \| TR<=50 | 234 |
| no truth variant nearby (pure FP) \| TR<=50 | 175 |
| allele/representation mismatch near a truth variant \| TR51-200 | 93 |
| no truth variant nearby (pure FP) \| other | 77 |
| no truth variant nearby (pure FP) \| TR51-200 | 75 |
| allele/representation mismatch near a truth variant \| other | 70 |
| genotype error (allele right, GT wrong) \| HP>=7 | 53 |
| genotype error (allele right, GT wrong) \| TR<=50 | 10 |
| genotype error (allele right, GT wrong) \| TR51-200 | 9 |
| genotype error (allele right, GT wrong) \| other | 3 |

**FP called genotype**

- allele/representation mismatch near a truth variant, called 1/1: 441
- no truth variant nearby, called 1/1: 377
- allele/representation mismatch near a truth variant, called 0/1: 366
- no truth variant nearby, called 0/1: 242
- allele/representation mismatch near a truth variant, called 2/1: 68
- genotype error, called 1/1: 46
- genotype error, called 0/1: 29
- no truth variant nearby, called 2/1: 4

**Recall: clustered (another truth variant within 10 bp) vs isolated**

- isolated: 9672/10352 = 0.9343
- clustered: 531/648 = 0.8194

**Median features: TP calls vs pure-FP calls**

| Feature | TP | pure FP |
|---|---|---|
| NALT | 15.00 | 7.00 |
| NREF | 9.00 | 0.00 |
| F | 0.58 | 1.00 |
| HQ | 110.43 | 75.64 |
| QMIN | 0.92 | 0.84 |

Pure FPs with NALT <= 2: 144/623

**Examples**

- *multi-allelic site (truth 1/2), alleles genotyped separately*
  - `chr20:152450 CAAA>CAAAA,C truth GT 2/1 [HP>=7] -> C>CA LowProb F=1;NALT=3;NREF=0;NOTH=14;CR=0.4772186005384677;RPRI=0.254734444286926; CAAA>C PASS F=1;NALT=10;NREF=0;NOTH=7;CR=0.4616893532706749;RPRI=0.230844676635337`
  - `chr20:232033 A>AAAATAAAT,AAAATAAATAAAT truth GT 1/2 [TR<=50] -> A>AAAATAAAT PASS F=1;NALT=6;NREF=0;NOTH=11;CR=0.28809519845678094;RPRI=0.14404761978907; A>AAAATAAATAAAT PASS F=1;NALT=9;NREF=0;NOTH=8;CR=0.28809519845678094;RPRI=0.144047599228412`
  - `chr20:348484 CGTGTGT>CGTGT,C truth GT 2/1 [TR<=50] -> CGT>C PASS F=1;NALT=6;NREF=0;NOTH=6;CR=0.3408628719454082;RPRI=0.1704314565333906; CGTGTGT>C PASS F=1;NALT=6;NREF=0;NOTH=6;CR=0.3408628719454082;RPRI=0.1704314565333906`
  - `chr20:365798 CTT>CT,C truth GT 2/1 [HP>=7] -> CTT>C PASS F=1;NALT=3;NREF=0;NOTH=1;CR=0.48763776212656507;RPRI=0.259944025080975`
- *genotype wrong (allele found, GT differs)*
  - `chr20:217716 ATG>A truth GT 0/1 [TR51-200] -> ATG>A PASS F=0.35;NALT=7;NREF=13;NOTH=0;CR=0.31961826329242865;RPRI=0.15990113792`
  - `chr20:454236 CAA>C truth GT 0/1 [HP>=7] -> CAA>C PASS F=1;NALT=6;NREF=0;NOTH=2;CR=0.4788248403179828;RPRI=0.2394124407196778`
  - `chr20:454295 G>GA truth GT 0/1 [HP>=7] -> G>GA PASS F=0.8666666666666667;NALT=13;NREF=2;NOTH=0;CR=0.29261482428466484;RPRI`
  - `chr20:649701 CA>C truth GT 0/1 [HP>=7] -> CA>C PASS F=0.35;NALT=7;NREF=13;NOTH=0;CR=0.41882219109579155;RPRI=0.20950310183`
- *different allele/representation called nearby*
  - `chr20:367246 AC>A truth GT 0/1 [TR<=50] -> AAAC>A LowProb F=0.45;NALT=9;NREF=11;NOTH=2;CR=0.11975493531830429;RPRI=0.05987747704; AAAC>A PASS F=0.14285714285714285;NALT=2;NREF=12;NOTH=9;CR=0.10932479580603646;RPR`
  - `chr20:367250 AAC>A truth GT 0/1 [TR<=50] -> AAAC>A LowProb F=0.45;NALT=9;NREF=11;NOTH=2;CR=0.11975493531830429;RPRI=0.05987747704; AAAC>A PASS F=0.14285714285714285;NALT=2;NREF=12;NOTH=9;CR=0.10932479580603646;RPR`
  - `chr20:499958 CA>C truth GT 0/1 [HP>=7] -> CA>C PASS F=0.6190476190476191;NALT=13;NREF=8;NOTH=5;CR=0.4447384505149794;RPRI=; CAA>C PASS F=0.2;NALT=2;NREF=8;NOTH=16;CR=0.43781690864807266;RPRI=0.219965659073`
  - `chr20:887636 CT>C truth GT 1/1 [HP>=7] -> CT>C PASS F=0.75;NALT=3;NREF=1;NOTH=3;CR=0.4947136621578156;RPRI=0.3755096907743`
- *never a candidate (discovery miss)*
  - `chr20:527001 C>CAGACAGAGAGAGAGAG truth GT 0/1 [TR51-200] -> no candidate`
  - `chr20:1404738 CTTT>C truth GT 0/1 [HP>=7] -> no candidate`
  - `chr20:1849802 AAAAT>A truth GT 0/1 [TR<=50] -> no candidate`
  - `chr20:1914466 AG>A truth GT 1/1 [other] -> no candidate`
- *candidate found but filtered (LowProb/UNCERTAIN)*
  - `chr20:3101547 CTTTTTTTT>CTTT,C truth GT 2/1 [HP>=7] -> CTTTTTTTT>C UNCERTAIN F=1;NALT=2;NREF=0;NOTH=1;CR=0.4296594755469182;RPRI=0.6401336083518908`
  - `chr20:4208484 CAA>C truth GT 0/1 [HP>=7] -> CA>C LowProb F=0.4;NALT=2;NREF=3;NOTH=3;CR=0.5813877621265651;RPRI=0.32189602370372; CAA>C LowProb F=0.4;NALT=2;NREF=3;NOTH=3;CR=0.5813877621265651;RPRI=0.32189602370372`
  - `chr20:13366867 T>TACAC truth GT 0/1 [TR<=50] -> C>CACA LowProb F=0.15384615384615385;NALT=2;NREF=11;NOTH=3;CR=0.05382614221358484;RPR`
  - `chr20:16944080 CA>C truth GT 0/1 [HP>=7] -> C>CA LowProb F=0.2222222222222222;NALT=2;NREF=7;NOTH=4;CR=0.47602637505431555;RPRI=; CA>C LowProb F=0.36363636363636365;NALT=4;NREF=7;NOTH=2;CR=0.4692852181881585;RPRI=`
- *no truth variant nearby (pure FP) | HP>=7*
  - `chr20:92968 GCT>G called 0/1 truth ./. NALT=2 NREF=4 F=0.3333333333333333 HQ=21.3271122905428`
  - `chr20:736736 CA>C called 1/1 truth ./. NALT=14 NREF=0 F=1 HQ=107.79971965386869`
  - `chr20:736736 CAA>C called 1/1 truth ./. NALT=21 NREF=0 F=1 HQ=130.36764723673866`
- *allele/representation mismatch near a truth variant | HP>=7*
  - `chr20:152450 CAAA>C called 0/1 truth ./. NALT=10 NREF=0 F=1 HQ=94.47319983496834`
  - `chr20:152467 AAGT>A called 0/1 truth ./. NALT=2 NREF=25 F=0.07407407407407407 HQ=24.951765070727788`
  - `chr20:365798 CTT>C called 1/1 truth ./. NALT=3 NREF=0 F=1 HQ=15.721706760708926`
- *genotype error (allele right, GT wrong) | TR51-200*
  - `chr20:217716 ATG>A called 1/1 truth 0/1 NALT=7 NREF=13 F=0.35 HQ=70.13553191129836`
  - `chr20:18488953 A>ATG called 1/1 truth 0/1 NALT=12 NREF=2 F=0.8571428571428571 HQ=118.54233840318568`
  - `chr20:21798505 C>CTG called 0/1 truth 1/1 NALT=12 NREF=5 F=0.7058823529411765 HQ=49.196487935442406`
- *allele/representation mismatch near a truth variant | TR<=50*
  - `chr20:232033 A>AAAATAAAT called 1/1 truth ./. NALT=6 NREF=0 F=1 HQ=233.26578032823483`
  - `chr20:232033 A>AAAATAAATAAAT called 1/1 truth ./. NALT=9 NREF=0 F=1 HQ=190.3963411436746`
  - `chr20:367249 AAAC>A called 0/1 truth ./. NALT=2 NREF=12 F=0.14285714285714285 HQ=32.93539789994429`
- *no truth variant nearby (pure FP) | other*
  - `chr20:336556 A>ATTTT called 0/1 truth ./. NALT=2 NREF=23 F=0.08 HQ=13.112103604661684`
  - `chr20:525486 AC>A called 0/1 truth ./. NALT=2 NREF=14 F=0.125 HQ=9.851752203758963`
  - `chr20:676477 TCA>T called 0/1 truth ./. NALT=2 NREF=18 F=0.1 HQ=10.246197031516317`
- *no truth variant nearby (pure FP) | TR<=50*
  - `chr20:348484 CGT>C called 1/1 truth ./. NALT=6 NREF=0 F=1 HQ=48.29047094875415`
  - `chr20:348484 CGTGTGT>C called 0/1 truth ./. NALT=6 NREF=0 F=1 HQ=66.14979742487202`
  - `chr20:685696 AAATAATAAT>A called 1/1 truth ./. NALT=12 NREF=0 F=1 HQ=284.71517459542525`
- *genotype error (allele right, GT wrong) | HP>=7*
  - `chr20:454236 CAA>C called 1/1 truth 0/1 NALT=6 NREF=0 F=1 HQ=40.06420972834334`
  - `chr20:454295 G>GA called 1/1 truth 0/1 NALT=13 NREF=2 F=0.8666666666666667 HQ=60.17710085310671`
  - `chr20:649701 CA>C called 1/1 truth 0/1 NALT=7 NREF=13 F=0.35 HQ=58.01055927094308`
- *allele/representation mismatch near a truth variant | TR51-200*
  - `chr20:957300 ATG>A called 0/1 truth ./. NALT=2 NREF=25 F=0.07407407407407407 HQ=12.511869896646203`
  - `chr20:2945582 GGTGTGTGTGT>GGT,G called 2/1 truth ./. NALT=9 NREF=11 F=0.45 HQ=158.18074326428302`
  - `chr20:3310333 A>AAAAGAAAG called 1/1 truth ./. NALT=5 NREF=0 F=1 HQ=317.6417349769573`
- *allele/representation mismatch near a truth variant | other*
  - `chr20:1012670 A>AG called 0/1 truth ./. NALT=2 NREF=19 F=0.09523809523809523 HQ=23.478562747468914`
  - `chr20:3079560 A>ATAT called 0/1 truth ./. NALT=2 NREF=22 F=0.08333333333333333 HQ=53.42672775584731`
  - `chr20:3188828 C>CTTT called 1/1 truth ./. NALT=2 NREF=0 F=1 HQ=131.61009333865184`
- *no truth variant nearby (pure FP) | TR51-200*
  - `chr20:2011231 CCTCTCTCT>C called 1/1 truth ./. NALT=11 NREF=0 F=1 HQ=366.12678194705927`
  - `chr20:2011231 CCTCTCTCTCTCTCTCT>C called 1/1 truth ./. NALT=11 NREF=0 F=1 HQ=188.05655501391834`
  - `chr20:4779167 CGTGTGT>C called 1/1 truth ./. NALT=6 NREF=0 F=1 HQ=276.1817286393695`
- *genotype error (allele right, GT wrong) | TR<=50*
  - `chr20:4411737 T>TTG called 0/1 truth 1/1 NALT=15 NREF=3 F=0.8333333333333334 HQ=86.20770891409506`
  - `chr20:5454164 AATATAT>A called 1/1 truth 0/1 NALT=12 NREF=1 F=0.9230769230769231 HQ=212.6975999211677`
  - `chr20:9553557 A>AGT called 1/1 truth 0/1 NALT=16 NREF=11 F=0.5925925925925926 HQ=113.41515530366468`
- *genotype error (allele right, GT wrong) | other*
  - `chr20:19711132 A>AAC called 0/1 truth 1/1 NALT=17 NREF=1 F=0.9444444444444444 HQ=37.21491202429797`
  - `chr20:37236093 T>TTGTGGGTGTGTG called 0/1 truth 1/1 NALT=12 NREF=1 F=0.9230769230769231 HQ=198.78066437009585`
  - `chr20:62420350 GCA>G called 1/1 truth 0/1 NALT=11 NREF=8 F=0.5789473684210527 HQ=45.76519985287564`

## Error taxonomy: expC_frag_30x


**False negatives by cause**

| Cause | Count |
|---|---|
| genotype wrong (allele found, GT differs) | 507 |
| multi-allelic site (truth 1/2), alleles genotyped separately | 366 |
| candidate found but filtered (LowProb/UNCERTAIN) | 121 |
| never a candidate (discovery miss) | 114 |
| different allele/representation called nearby | 88 |

**False positives by cause and repeat context**

| Cause \| context | Count |
|---|---|
| allele/representation mismatch near a truth variant \| HP>=7 | 309 |
| no truth variant nearby (pure FP) \| HP>=7 | 211 |
| allele/representation mismatch near a truth variant \| TR<=50 | 190 |
| no truth variant nearby (pure FP) \| TR<=50 | 149 |
| genotype error (allele right, GT wrong) \| HP>=7 | 108 |
| allele/representation mismatch near a truth variant \| TR51-200 | 69 |
| no truth variant nearby (pure FP) \| TR51-200 | 55 |
| no truth variant nearby (pure FP) \| other | 37 |
| allele/representation mismatch near a truth variant \| other | 32 |
| genotype error (allele right, GT wrong) \| TR<=50 | 27 |
| genotype error (allele right, GT wrong) \| TR51-200 | 20 |
| genotype error (allele right, GT wrong) \| other | 20 |

**FP called genotype**

- no truth variant nearby, called 1/1: 284
- allele/representation mismatch near a truth variant, called 0/1: 282
- allele/representation mismatch near a truth variant, called 1/1: 279
- no truth variant nearby, called 0/1: 167
- genotype error, called 0/1: 107
- genotype error, called 1/1: 68
- allele/representation mismatch near a truth variant, called 2/1: 39
- no truth variant nearby, called 2/1: 1

**Recall: clustered (another truth variant within 10 bp) vs isolated**

- isolated: 9294/10352 = 0.8978
- clustered: 510/648 = 0.7870

**Median features: TP calls vs pure-FP calls**

| Feature | TP | pure FP |
|---|---|---|
| NALT | 13.00 | 8.00 |
| NREF | 9.00 | 0.00 |
| F | 0.58 | 1.00 |
| HQ | 103.08 | 67.75 |
| QMIN | 0.93 | 0.85 |

Pure FPs with NALT <= 2: 82/452

**Examples**

- *genotype wrong (allele found, GT differs)*
  - `chr20:92969 CTT>C truth GT 1/1 [HP>=7] -> CTT>C LowProb F=1;NALT=6;NREF=0;NOTH=2;CR=0.500026912053441;RPRI=0.25001347658740697`
  - `chr20:126062 C>CT truth GT 1/1 [HP>=7] -> C>CT LowProb F=1;NALT=23;NREF=0;NOTH=0;CR=0.39569967100534;RPRI=0.19784983550267;TA`
  - `chr20:217716 ATG>A truth GT 0/1 [TR51-200] -> ATG>A PASS F=0.375;NALT=6;NREF=10;NOTH=0;CR=0.33211826329242866;RPRI=0.1666514764`
  - `chr20:240075 C>CAAA truth GT 1/1 [HP>=7] -> C>CAAA PASS F=1;NALT=12;NREF=0;NOTH=1;CR=0.4973683769954572;RPRI=0.248684188497728`
- *candidate found but filtered (LowProb/UNCERTAIN)*
  - `chr20:290494 G>GC truth GT 1/1 [other] -> G>GC LowProb F=1;NALT=10;NREF=0;NOTH=0;CR=0.12340895967450088;RPRI=0.06170447983725`
  - `chr20:365798 CTT>CT,C truth GT 2/1 [HP>=7] -> CTT>C UNCERTAIN F=1;NALT=2;NREF=0;NOTH=1;CR=0.42513776212656507;RPRI=0.637872751641714`
  - `chr20:367246 AC>A truth GT 0/1 [TR<=50] -> AAAC>A LowProb F=0.5;NALT=9;NREF=9;NOTH=2;CR=0.12089129895466792;RPRI=0.0604456518720; AAAC>A LowProb F=0.18181818181818182;NALT=2;NREF=9;NOTH=9;CR=0.11095523058864516;RPRI`
  - `chr20:367250 AAC>A truth GT 0/1 [TR<=50] -> AAAC>A LowProb F=0.5;NALT=9;NREF=9;NOTH=2;CR=0.12089129895466792;RPRI=0.0604456518720; AAAC>A LowProb F=0.18181818181818182;NALT=2;NREF=9;NOTH=9;CR=0.11095523058864516;RPRI`
- *multi-allelic site (truth 1/2), alleles genotyped separately*
  - `chr20:348484 CGTGTGT>CGTGT,C truth GT 2/1 [TR<=50] -> CGT>C PASS F=1;NALT=6;NREF=0;NOTH=5;CR=0.3484386295211658;RPRI=0.1742193353212693; CGTGTGT>C PASS F=1;NALT=5;NREF=0;NOTH=6;CR=0.3484386295211658;RPRI=0.1742212881357072`
  - `chr20:511703 TACAC>CACAC,T truth GT 2/1 [TR<=50] -> TACAC>T PASS F=0.4230769230769231;NALT=11;NREF=15;NOTH=0;CR=0.18949029730800024;RPR`
  - `chr20:571051 C>CT,CTTT truth GT 1/2 [HP>=7] -> C>CT PASS F=0.875;NALT=7;NREF=1;NOTH=4;CR=0.5413357980709209;RPRI=0.270667910592; C>CTTT LowProb F=0.6666666666666666;NALT=2;NREF=1;NOTH=9;CR=0.5413357980709209;RPRI=0`
  - `chr20:685696 AAATAATAATAAT>AAAT,A truth GT 2/1 [TR<=50] -> AAATAATAAT>A PASS F=1;NALT=10;NREF=0;NOTH=5;CR=0.3587215551410115;RPRI=0.179360777570505; AAATAATAATAAT>A PASS F=1;NALT=5;NREF=0;NOTH=10;CR=0.3587215551410115;RPRI=0.179362750945630`
- *different allele/representation called nearby*
  - `chr20:499958 CA>C truth GT 0/1 [HP>=7] -> CA>C PASS F=0.6;NALT=12;NREF=8;NOTH=3;CR=0.44975517292300615;RPRI=0.224877586461; CAA>C PASS F=0.2;NALT=2;NREF=8;NOTH=13;CR=0.4428336310560994;RPRI=0.2224740202775`
  - `chr20:1031180 TA>T truth GT 0/1 [HP>=7] -> TA>T PASS F=0.5;NALT=8;NREF=8;NOTH=4;CR=0.3288611456159177;RPRI=0.16443063328095; TAA>T PASS F=0.2727272727272727;NALT=3;NREF=8;NOTH=9;CR=0.3163611456159177;RPRI=0`
  - `chr20:2597654 ATTTT>A truth GT 0/1 [HP>=7] -> AT>A PASS F=0.16666666666666666;NALT=2;NREF=10;NOTH=15;CR=0.4388870328829897;RPR; ATTTT>A PASS F=0.5238095238095238;NALT=11;NREF=10;NOTH=6;CR=0.43348978200302;RPRI=0`
  - `chr20:2813920 C>CAA truth GT 0/1 [HP>=7] -> C>CAA PASS F=0.3333333333333333;NALT=7;NREF=14;NOTH=6;CR=0.34551929155028693;RPRI; CA>C PASS F=0.17647058823529413;NALT=3;NREF=14;NOTH=10;CR=0.34551929155028693;RP`
- *never a candidate (discovery miss)*
  - `chr20:527001 C>CAGACAGAGAGAGAGAG truth GT 0/1 [TR51-200] -> no candidate`
  - `chr20:1404738 CTTT>C truth GT 0/1 [HP>=7] -> no candidate`
  - `chr20:1849802 AAAAT>A truth GT 0/1 [TR<=50] -> no candidate`
  - `chr20:1914466 AG>A truth GT 1/1 [other] -> no candidate`
- *genotype error (allele right, GT wrong) | TR51-200*
  - `chr20:217716 ATG>A called 1/1 truth 0/1 NALT=6 NREF=10 F=0.375 HQ=10.273894107785338`
  - `chr20:10507904 T>TTA called 1/1 truth 0/1 NALT=15 NREF=0 F=1 HQ=160.77308847223787`
  - `chr20:11498279 A>AGGAAGAAAGAGAGGGAGGAAAAGAGAGAGGAAGGGAGGG called 1/1 truth 0/1 NALT=5 NREF=13 F=0.2777777777777778 HQ=32.45933489323086`
- *genotype error (allele right, GT wrong) | HP>=7*
  - `chr20:240075 C>CAAA called 0/1 truth 1/1 NALT=12 NREF=0 F=1 HQ=12.605009281055423`
  - `chr20:240477 G>GA called 0/1 truth 1/1 NALT=14 NREF=1 F=0.9333333333333333 HQ=33.5932354891637`
  - `chr20:1065654 C>CT called 0/1 truth 1/1 NALT=13 NREF=1 F=0.9285714285714286 HQ=13.314328610528293`
- *no truth variant nearby (pure FP) | TR<=50*
  - `chr20:348484 CGT>C called 1/1 truth ./. NALT=6 NREF=0 F=1 HQ=59.06917788176776`
  - `chr20:348484 CGTGTGT>C called 0/1 truth ./. NALT=5 NREF=0 F=1 HQ=28.589456175165992`
  - `chr20:685696 AAATAATAAT>A called 1/1 truth ./. NALT=10 NREF=0 F=1 HQ=177.25942593297205`
- *allele/representation mismatch near a truth variant | HP>=7*
  - `chr20:499958 CAA>CA,C called 2/1 truth ./. NALT=12 NREF=8 F=0.6 HQ=17.59318640608859`
  - `chr20:571051 C>CT called 1/1 truth ./. NALT=7 NREF=1 F=0.875 HQ=16.601034214225734`
  - `chr20:787529 TG>T called 0/1 truth ./. NALT=2 NREF=16 F=0.1111111111111111 HQ=18.851557469219756`
- *allele/representation mismatch near a truth variant | TR<=50*
  - `chr20:511703 TACAC>T called 0/1 truth ./. NALT=11 NREF=15 F=0.4230769230769231 HQ=136.26400736191738`
  - `chr20:734610 G>GGT called 1/1 truth ./. NALT=7 NREF=0 F=1 HQ=82.81052873469099`
  - `chr20:734610 G>GGTGTGTGT called 0/1 truth ./. NALT=10 NREF=0 F=1 HQ=88.38606121842673`
- *allele/representation mismatch near a truth variant | TR51-200*
  - `chr20:672388 CGT>C called 0/1 truth ./. NALT=2 NREF=13 F=0.13333333333333333 HQ=10.125237969989598`
  - `chr20:3310333 A>AAAAGAAAG called 1/1 truth ./. NALT=4 NREF=0 F=1 HQ=265.2534845535227`
  - `chr20:3310333 A>AAAAGAAAGAAAGAAAG called 0/1 truth ./. NALT=12 NREF=0 F=1 HQ=135.05090203214718`
- *no truth variant nearby (pure FP) | HP>=7*
  - `chr20:736736 CA>C called 1/1 truth ./. NALT=12 NREF=0 F=1 HQ=118.01519660762108`
  - `chr20:736736 CAA>C called 1/1 truth ./. NALT=14 NREF=0 F=1 HQ=126.38477738565155`
  - `chr20:1700552 TA>T called 0/1 truth ./. NALT=3 NREF=13 F=0.1875 HQ=12.008635388879611`
- *no truth variant nearby (pure FP) | TR51-200*
  - `chr20:2011231 CCTCTCTCT>C called 1/1 truth ./. NALT=10 NREF=0 F=1 HQ=305.4728964078706`
  - `chr20:2011231 CCTCTCTCTCTCTCTCT>C called 1/1 truth ./. NALT=9 NREF=0 F=1 HQ=77.79438314335516`
  - `chr20:4779167 CGTGTGT>C called 1/1 truth ./. NALT=5 NREF=0 F=1 HQ=171.52319972908623`
- *genotype error (allele right, GT wrong) | TR<=50*
  - `chr20:2260781 AAAAC>A called 0/1 truth 1/1 NALT=23 NREF=0 F=1 HQ=14.373078906267315`
  - `chr20:4411737 T>TTG called 0/1 truth 1/1 NALT=15 NREF=2 F=0.8823529411764706 HQ=145.83375454229548`
  - `chr20:5454164 AATATAT>A called 1/1 truth 0/1 NALT=11 NREF=1 F=0.9166666666666666 HQ=176.81080245484281`
- *allele/representation mismatch near a truth variant | other*
  - `chr20:2674873 AGT>A called 0/1 truth ./. NALT=2 NREF=12 F=0.14285714285714285 HQ=9.561321330225951`
  - `chr20:3188828 C>CTTT called 1/1 truth ./. NALT=2 NREF=0 F=1 HQ=99.0162731452632`
  - `chr20:4619359 TC>T called 1/1 truth ./. NALT=2 NREF=35 F=0.05405405405405406 HQ=332.5458513587737`
- *genotype error (allele right, GT wrong) | other*
  - `chr20:3273603 CA>C called 1/1 truth 0/1 NALT=18 NREF=9 F=0.6666666666666666 HQ=39.16222261261348`
  - `chr20:8495564 A>AT called 1/1 truth 0/1 NALT=12 NREF=4 F=0.75 HQ=76.62369132292645`
  - `chr20:10631606 AT>A called 1/1 truth 0/1 NALT=7 NREF=5 F=0.5833333333333334 HQ=39.14553129562921`
- *no truth variant nearby (pure FP) | other*
  - `chr20:4374492 CA>C called 0/1 truth ./. NALT=5 NREF=15 F=0.25 HQ=39.134160585044526`
  - `chr20:4374503 CTACCTCCAACCTT>C called 0/1 truth ./. NALT=2 NREF=12 F=0.14285714285714285 HQ=22.6358644216337`
  - `chr20:25887339 AG>A called 0/1 truth ./. NALT=4 NREF=23 F=0.14814814814814814 HQ=20.38825257926294`

## Error taxonomy: expC2_fragend_30x


**False negatives by cause**

| Cause | Count |
|---|---|
| genotype wrong (allele found, GT differs) | 492 |
| multi-allelic site (truth 1/2), alleles genotyped separately | 369 |
| candidate found but filtered (LowProb/UNCERTAIN) | 117 |
| never a candidate (discovery miss) | 111 |
| different allele/representation called nearby | 93 |

**False positives by cause and repeat context**

| Cause \| context | Count |
|---|---|
| allele/representation mismatch near a truth variant \| HP>=7 | 396 |
| no truth variant nearby (pure FP) \| HP>=7 | 221 |
| allele/representation mismatch near a truth variant \| TR<=50 | 188 |
| no truth variant nearby (pure FP) \| TR<=50 | 149 |
| genotype error (allele right, GT wrong) \| HP>=7 | 95 |
| allele/representation mismatch near a truth variant \| TR51-200 | 66 |
| no truth variant nearby (pure FP) \| TR51-200 | 57 |
| no truth variant nearby (pure FP) \| other | 36 |
| allele/representation mismatch near a truth variant \| other | 35 |
| genotype error (allele right, GT wrong) \| other | 28 |
| genotype error (allele right, GT wrong) \| TR<=50 | 25 |
| genotype error (allele right, GT wrong) \| TR51-200 | 23 |

**FP called genotype**

- allele/representation mismatch near a truth variant, called 0/1: 357
- no truth variant nearby, called 1/1: 291
- allele/representation mismatch near a truth variant, called 1/1: 285
- no truth variant nearby, called 0/1: 171
- genotype error, called 0/1: 96
- genotype error, called 1/1: 75
- allele/representation mismatch near a truth variant, called 2/1: 43
- no truth variant nearby, called 2/1: 1

**Recall: clustered (another truth variant within 10 bp) vs isolated**

- isolated: 9310/10352 = 0.8993
- clustered: 508/648 = 0.7840

**Median features: TP calls vs pure-FP calls**

| Feature | TP | pure FP |
|---|---|---|
| NALT | 13.00 | 8.00 |
| NREF | 9.00 | 0.00 |
| F | 0.58 | 1.00 |
| HQ | 108.10 | 77.81 |
| QMIN | 0.93 | 0.84 |

Pure FPs with NALT <= 2: 85/463

**Examples**

- *genotype wrong (allele found, GT differs)*
  - `chr20:126062 C>CT truth GT 1/1 [HP>=7] -> C>CT LowProb F=1;NALT=23;NREF=0;NOTH=0;CR=0.39569967100534;RPRI=0.19784983550267;TA`
  - `chr20:217716 ATG>A truth GT 0/1 [TR51-200] -> ATG>A PASS F=0.375;NALT=6;NREF=10;NOTH=0;CR=0.33211826329242866;RPRI=0.1666514764`
  - `chr20:240075 C>CAAA truth GT 1/1 [HP>=7] -> C>CAAA PASS F=1;NALT=12;NREF=0;NOTH=1;CR=0.4973683769954572;RPRI=0.248684188497728`
  - `chr20:566178 AAAAAAT>A truth GT 1/1 [HP>=7] -> AAAAAAT>A LowProb F=1;NALT=18;NREF=0;NOTH=0;CR=0.3012729609356477;RPRI=0.150636480467823`
- *candidate found but filtered (LowProb/UNCERTAIN)*
  - `chr20:290494 G>GC truth GT 1/1 [other] -> G>GC LowProb F=1;NALT=10;NREF=0;NOTH=0;CR=0.12340895967450088;RPRI=0.06170447983725`
  - `chr20:365798 CTT>CT,C truth GT 2/1 [HP>=7] -> CTT>C UNCERTAIN F=1;NALT=2;NREF=0;NOTH=1;CR=0.5084710954598984;RPRI=0.6795394183083808`
  - `chr20:586045 CT>CTT,C truth GT 2/1 [HP>=7] -> CT>C LowProb F=1;NALT=3;NREF=0;NOTH=6;CR=0.6068173838569234;RPRI=0.3195338359461542`
  - `chr20:787512 CT>C truth GT 1/1 [HP>=7] -> CT>C LowProb F=1;NALT=10;NREF=0;NOTH=2;CR=0.4838248403179828;RPRI=0.241912420158991`
- *multi-allelic site (truth 1/2), alleles genotyped separately*
  - `chr20:348484 CGTGTGT>CGTGT,C truth GT 2/1 [TR<=50] -> CGT>C PASS F=1;NALT=6;NREF=0;NOTH=5;CR=0.3484386295211658;RPRI=0.1742193353212693; CGTGTGT>C PASS F=1;NALT=5;NREF=0;NOTH=6;CR=0.3484386295211658;RPRI=0.1742212881357072`
  - `chr20:511703 TACAC>CACAC,T truth GT 2/1 [TR<=50] -> TACAC>T PASS F=0.4230769230769231;NALT=11;NREF=15;NOTH=0;CR=0.18949029730800024;RPR`
  - `chr20:571051 C>CT,CTTT truth GT 1/2 [HP>=7] -> C>CT PASS F=0.875;NALT=7;NREF=1;NOTH=4;CR=0.5413357980709209;RPRI=0.270667910592; C>CTTT LowProb F=0.6666666666666666;NALT=2;NREF=1;NOTH=9;CR=0.5413357980709209;RPRI=0`
  - `chr20:685696 AAATAATAATAAT>AAAT,A truth GT 2/1 [TR<=50] -> AAATAATAAT>A PASS F=1;NALT=10;NREF=0;NOTH=5;CR=0.3587215551410115;RPRI=0.179360777570505; AAATAATAATAAT>A PASS F=1;NALT=5;NREF=0;NOTH=10;CR=0.3587215551410115;RPRI=0.179362750945630`
- *different allele/representation called nearby*
  - `chr20:367246 AC>A truth GT 0/1 [TR<=50] -> AAAC>A LowProb F=0.5;NALT=9;NREF=9;NOTH=2;CR=0.12089129895466792;RPRI=0.0604456518720; AAAC>A PASS F=0.18181818181818182;NALT=2;NREF=9;NOTH=9;CR=0.09845523058864517;RPRI`
  - `chr20:367250 AAC>A truth GT 0/1 [TR<=50] -> AAAC>A LowProb F=0.5;NALT=9;NREF=9;NOTH=2;CR=0.12089129895466792;RPRI=0.0604456518720; AAAC>A PASS F=0.18181818181818182;NALT=2;NREF=9;NOTH=9;CR=0.09845523058864517;RPRI`
  - `chr20:499958 CA>C truth GT 0/1 [HP>=7] -> CA>C PASS F=0.6;NALT=12;NREF=8;NOTH=3;CR=0.44975517292300615;RPRI=0.224877586461; CAA>C PASS F=0.2;NALT=2;NREF=8;NOTH=13;CR=0.4428336310560994;RPRI=0.2224740202775`
  - `chr20:1031180 TA>T truth GT 0/1 [HP>=7] -> TA>T PASS F=0.5;NALT=8;NREF=8;NOTH=4;CR=0.3288611456159177;RPRI=0.16443063328095; TAA>T PASS F=0.2727272727272727;NALT=3;NREF=8;NOTH=9;CR=0.3163611456159177;RPRI=0`
- *never a candidate (discovery miss)*
  - `chr20:527001 C>CAGACAGAGAGAGAGAG truth GT 0/1 [TR51-200] -> no candidate`
  - `chr20:1404738 CTTT>C truth GT 0/1 [HP>=7] -> no candidate`
  - `chr20:1849802 AAAAT>A truth GT 0/1 [TR<=50] -> no candidate`
  - `chr20:1914466 AG>A truth GT 1/1 [other] -> no candidate`
- *no truth variant nearby (pure FP) | HP>=7*
  - `chr20:92968 GCT>G called 1/1 truth ./. NALT=2 NREF=3 F=0.4 HQ=47.601198377543824`
  - `chr20:736735 TC>T called 0/1 truth ./. NALT=2 NREF=2 F=0.5 HQ=9.26132044553367`
  - `chr20:736736 CA>C called 1/1 truth ./. NALT=11 NREF=0 F=1 HQ=118.42202463764636`
- *genotype error (allele right, GT wrong) | TR51-200*
  - `chr20:217716 ATG>A called 1/1 truth 0/1 NALT=6 NREF=10 F=0.375 HQ=20.28217599363785`
  - `chr20:10507904 T>TTA called 1/1 truth 0/1 NALT=14 NREF=0 F=1 HQ=158.46083572225325`
  - `chr20:10507919 T>TAC called 1/1 truth 0/1 NALT=10 NREF=15 F=0.4 HQ=124.43981718647262`
- *genotype error (allele right, GT wrong) | HP>=7*
  - `chr20:240075 C>CAAA called 0/1 truth 1/1 NALT=12 NREF=0 F=1 HQ=12.605009281055423`
  - `chr20:1122295 CTTTT>C called 0/1 truth 1/1 NALT=11 NREF=1 F=0.9166666666666666 HQ=108.63895804692868`
  - `chr20:1138518 C>CA called 0/1 truth 1/1 NALT=7 NREF=3 F=0.7 HQ=45.697963314081335`
- *allele/representation mismatch near a truth variant | TR<=50*
  - `chr20:319668 CAA>C called 0/1 truth ./. NALT=2 NREF=24 F=0.07692307692307693 HQ=43.65754438631882`
  - `chr20:367249 AAAC>A called 0/1 truth ./. NALT=2 NREF=9 F=0.18181818181818182 HQ=110.36699803376086`
  - `chr20:511703 TACAC>T called 0/1 truth ./. NALT=11 NREF=15 F=0.4230769230769231 HQ=136.264132324064`
- *no truth variant nearby (pure FP) | TR<=50*
  - `chr20:348484 CGT>C called 1/1 truth ./. NALT=6 NREF=0 F=1 HQ=36.26445605415182`
  - `chr20:348484 CGTGTGT>C called 0/1 truth ./. NALT=5 NREF=0 F=1 HQ=73.86744067274412`
  - `chr20:685696 AAATAATAAT>A called 1/1 truth ./. NALT=10 NREF=0 F=1 HQ=177.25942593297205`
- *allele/representation mismatch near a truth variant | HP>=7*
  - `chr20:499958 CAA>CA,C called 2/1 truth ./. NALT=12 NREF=8 F=0.6 HQ=44.53543870614999`
  - `chr20:560085 AC>A called 0/1 truth ./. NALT=2 NREF=1 F=0.6666666666666666 HQ=17.438918597385282`
  - `chr20:571051 C>CT called 1/1 truth ./. NALT=7 NREF=1 F=0.875 HQ=16.601034214225734`
- *allele/representation mismatch near a truth variant | TR51-200*
  - `chr20:672388 CGT>C called 0/1 truth ./. NALT=2 NREF=13 F=0.13333333333333333 HQ=10.125237969989598`
  - `chr20:3310333 A>AAAAGAAAG called 1/1 truth ./. NALT=4 NREF=0 F=1 HQ=265.25347294434306`
  - `chr20:3310333 A>AAAAGAAAGAAAGAAAG called 0/1 truth ./. NALT=12 NREF=0 F=1 HQ=135.05090203214718`
- *no truth variant nearby (pure FP) | TR51-200*
  - `chr20:2011231 CCTCTCTCT>C called 1/1 truth ./. NALT=10 NREF=0 F=1 HQ=305.4728964078706`
  - `chr20:2011231 CCTCTCTCTCTCTCTCT>C called 1/1 truth ./. NALT=9 NREF=0 F=1 HQ=77.96252776608253`
  - `chr20:4779167 CGTGTGT>C called 1/1 truth ./. NALT=5 NREF=0 F=1 HQ=171.52319972908623`
- *genotype error (allele right, GT wrong) | TR<=50*
  - `chr20:2260781 AAAAC>A called 0/1 truth 1/1 NALT=23 NREF=0 F=1 HQ=14.373078906267315`
  - `chr20:4411737 T>TTG called 0/1 truth 1/1 NALT=15 NREF=2 F=0.8823529411764706 HQ=145.83375454229548`
  - `chr20:5454164 AATATAT>A called 1/1 truth 0/1 NALT=11 NREF=1 F=0.9166666666666666 HQ=217.36012484478601`
- *allele/representation mismatch near a truth variant | other*
  - `chr20:3188828 C>CTTT called 1/1 truth ./. NALT=2 NREF=0 F=1 HQ=99.0162731452632`
  - `chr20:4619359 TC>T called 1/1 truth ./. NALT=2 NREF=35 F=0.05405405405405406 HQ=332.54585160358283`
  - `chr20:9049617 TG>T called 0/1 truth ./. NALT=2 NREF=0 F=1 HQ=27.406260015743328`
- *no truth variant nearby (pure FP) | other*
  - `chr20:4374492 CA>C called 0/1 truth ./. NALT=5 NREF=15 F=0.25 HQ=39.134160585044526`
  - `chr20:4374503 CTACCTCCAACCTT>C called 0/1 truth ./. NALT=2 NREF=11 F=0.15384615384615385 HQ=24.003630015080773`
  - `chr20:22457785 G>GTTTTT called 0/1 truth ./. NALT=2 NREF=17 F=0.10526315789473684 HQ=44.03547311742591`
- *genotype error (allele right, GT wrong) | other*
  - `chr20:5206755 A>AC called 1/1 truth 0/1 NALT=11 NREF=13 F=0.4583333333333333 HQ=26.130730136124118`
  - `chr20:10631606 AT>A called 1/1 truth 0/1 NALT=7 NREF=5 F=0.5833333333333334 HQ=39.14553129562921`
  - `chr20:11273742 G>GCC called 0/1 truth 1/1 NALT=33 NREF=0 F=1 HQ=13.13347253281688`

## Error taxonomy: expAB_30x


**False negatives by cause**

| Cause | Count |
|---|---|
| multi-allelic site (truth 1/2), alleles genotyped separately | 158 |
| different allele/representation called nearby | 120 |
| never a candidate (discovery miss) | 116 |
| genotype wrong (allele found, GT differs) | 95 |
| candidate found but filtered (LowProb/UNCERTAIN) | 17 |

**False positives by cause and repeat context**

| Cause \| context | Count |
|---|---|
| allele/representation mismatch near a truth variant \| HP>=7 | 252 |
| allele/representation mismatch near a truth variant \| TR<=50 | 94 |
| genotype error (allele right, GT wrong) \| HP>=7 | 76 |
| no truth variant nearby (pure FP) \| other | 36 |
| allele/representation mismatch near a truth variant \| TR51-200 | 31 |
| allele/representation mismatch near a truth variant \| other | 28 |
| genotype error (allele right, GT wrong) \| TR<=50 | 8 |
| genotype error (allele right, GT wrong) \| TR51-200 | 5 |
| no truth variant nearby (pure FP) \| HP>=7 | 4 |
| genotype error (allele right, GT wrong) \| other | 3 |
| no truth variant nearby (pure FP) \| TR51-200 | 1 |
| no truth variant nearby (pure FP) \| TR<=50 | 1 |

**FP called genotype**

- allele/representation mismatch near a truth variant, called 0/1: 303
- allele/representation mismatch near a truth variant, called 1/1: 74
- genotype error, called 0/1: 58
- no truth variant nearby, called 0/1: 42
- genotype error, called 1/1: 34
- allele/representation mismatch near a truth variant, called 2/1: 28

**Recall: clustered (another truth variant within 10 bp) vs isolated**

- isolated: 9985/10352 = 0.9645
- clustered: 509/648 = 0.7855

**Median features: TP calls vs pure-FP calls**

| Feature | TP | pure FP |
|---|---|---|
| NALT | n/a | n/a |
| NREF | n/a | n/a |
| F | n/a | n/a |
| HQ | n/a | n/a |
| QMIN | n/a | n/a |

**Examples**

- *genotype wrong (allele found, GT differs)*
  - `chr20:217716 ATG>A truth GT 0/1 [TR51-200] -> ATG>A PASS K=1;DP=20`
  - `chr20:454236 CAA>C truth GT 0/1 [HP>=7] -> CAA>C PASS K=1;DP=8`
  - `chr20:454295 G>GA truth GT 0/1 [HP>=7] -> G>GA PASS K=1;DP=15`
  - `chr20:649701 CA>C truth GT 0/1 [HP>=7] -> CA>C PASS K=1;DP=20`
- *multi-allelic site (truth 1/2), alleles genotyped separately*
  - `chr20:365798 CTT>CT,C truth GT 2/1 [HP>=7] -> CTT>C PASS K=1;DP=4`
  - `chr20:511703 TACAC>CACAC,T truth GT 2/1 [TR<=50] -> TACAC>T PASS K=1;DP=31`
  - `chr20:586045 CT>CTT,C truth GT 2/1 [HP>=7] -> CT>C PASS K=1;DP=9`
  - `chr20:942277 A>AT,T truth GT 1/2 [TR<=50] -> A>AT PASS K=2;DP=26`
- *different allele/representation called nearby*
  - `chr20:367246 AC>A truth GT 0/1 [TR<=50] -> AAAC>A PASS K=2;DP=23`
  - `chr20:367250 AAC>A truth GT 0/1 [TR<=50] -> AAAC>A PASS K=2;DP=23`
  - `chr20:499958 CA>C truth GT 0/1 [HP>=7] -> CAA>C PASS K=2;DP=26`
  - `chr20:887636 CT>C truth GT 1/1 [HP>=7] -> CT>C PASS K=1;DP=7`
- *never a candidate (discovery miss)*
  - `chr20:527001 C>CAGACAGAGAGAGAGAG truth GT 0/1 [TR51-200] -> no candidate`
  - `chr20:1404738 CTTT>C truth GT 0/1 [HP>=7] -> no candidate`
  - `chr20:1849802 AAAAT>A truth GT 0/1 [TR<=50] -> no candidate`
  - `chr20:1914466 AG>A truth GT 1/1 [other] -> no candidate`
- *candidate found but filtered (LowProb/UNCERTAIN)*
  - `chr20:4208484 CAA>C truth GT 0/1 [HP>=7] -> CA>C LowProb K=2;DP=8`
  - `chr20:11769171 CTT>CT,C truth GT 2/1 [HP>=7] -> CT>C LowProb K=1;DP=5`
  - `chr20:13366867 T>TACAC truth GT 0/1 [TR<=50] -> C>CACA LowProb K=1;DP=15`
  - `chr20:23554017 CA>CAA,C truth GT 2/1 [HP>=7] -> C>CA LowProb K=1;DP=2`
- *genotype error (allele right, GT wrong) | TR51-200*
  - `chr20:217716 ATG>A called 1/1 truth 0/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:21798505 C>CTG called 0/1 truth 1/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:35321252 C>CTG called 1/1 truth 0/1 NALT=None NREF=None F=None HQ=None`
- *allele/representation mismatch near a truth variant | HP>=7*
  - `chr20:365798 CTT>C called 1/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:458330 TA>T called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:499958 CAA>C called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *allele/representation mismatch near a truth variant | TR<=50*
  - `chr20:367249 AAAC>A called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:511703 TACAC>T called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:942277 A>AT called 1/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *genotype error (allele right, GT wrong) | HP>=7*
  - `chr20:454236 CAA>C called 1/1 truth 0/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:454295 G>GA called 1/1 truth 0/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:649701 CA>C called 1/1 truth 0/1 NALT=None NREF=None F=None HQ=None`
- *allele/representation mismatch near a truth variant | TR51-200*
  - `chr20:3982103 GTC>G called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:3982105 C>CTCTCTCTCTC called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:6124250 G>GCA,GCACA called 2/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *no truth variant nearby (pure FP) | HP>=7*
  - `chr20:4208933 TG>T called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:25926126 CA>C called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:26077020 GA>G called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *no truth variant nearby (pure FP) | other*
  - `chr20:4374492 CA>C called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:4374503 CTACCTCCAACCTT>C called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:10699738 A>AC called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *genotype error (allele right, GT wrong) | TR<=50*
  - `chr20:4411737 T>TTG called 0/1 truth 1/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:9553557 A>AGT called 1/1 truth 0/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:40040821 G>GTCTC called 1/1 truth 0/1 NALT=None NREF=None F=None HQ=None`
- *allele/representation mismatch near a truth variant | other*
  - `chr20:4619359 TC>T called 1/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:9051911 T>TA called 1/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:13929964 TA>T called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *no truth variant nearby (pure FP) | TR51-200*
  - `chr20:8490860 GCATATATATATATATATGTACACATATCTATATGTACACATATATATGTGTA>G called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *genotype error (allele right, GT wrong) | other*
  - `chr20:19711132 A>AAC called 0/1 truth 1/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:37236093 T>TTGTGGGTGTGTG called 0/1 truth 1/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:62420350 GCA>G called 1/1 truth 0/1 NALT=None NREF=None F=None HQ=None`
- *no truth variant nearby (pure FP) | TR<=50*
  - `chr20:26079727 ACTAATCACTGGCTGCAACT>A called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`

## Error taxonomy: expABC2_30x


**False negatives by cause**

| Cause | Count |
|---|---|
| genotype wrong (allele found, GT differs) | 458 |
| candidate found but filtered (LowProb/UNCERTAIN) | 271 |
| never a candidate (discovery miss) | 130 |
| multi-allelic site (truth 1/2), alleles genotyped separately | 126 |
| different allele/representation called nearby | 86 |

**False positives by cause and repeat context**

| Cause \| context | Count |
|---|---|
| allele/representation mismatch near a truth variant \| HP>=7 | 157 |
| allele/representation mismatch near a truth variant \| TR<=50 | 82 |
| genotype error (allele right, GT wrong) \| HP>=7 | 54 |
| no truth variant nearby (pure FP) \| other | 29 |
| allele/representation mismatch near a truth variant \| TR51-200 | 26 |
| genotype error (allele right, GT wrong) \| other | 23 |
| allele/representation mismatch near a truth variant \| other | 17 |
| genotype error (allele right, GT wrong) \| TR51-200 | 10 |
| genotype error (allele right, GT wrong) \| TR<=50 | 10 |
| no truth variant nearby (pure FP) \| HP>=7 | 2 |

**FP called genotype**

- allele/representation mismatch near a truth variant, called 0/1: 199
- allele/representation mismatch near a truth variant, called 1/1: 70
- genotype error, called 1/1: 62
- genotype error, called 0/1: 35
- no truth variant nearby, called 0/1: 31
- allele/representation mismatch near a truth variant, called 2/1: 13

**Recall: clustered (another truth variant within 10 bp) vs isolated**

- isolated: 9444/10352 = 0.9123
- clustered: 485/648 = 0.7485

**Median features: TP calls vs pure-FP calls**

| Feature | TP | pure FP |
|---|---|---|
| NALT | n/a | n/a |
| NREF | n/a | n/a |
| F | n/a | n/a |
| HQ | n/a | n/a |
| QMIN | n/a | n/a |

**Examples**

- *different allele/representation called nearby*
  - `chr20:92969 CTT>C truth GT 1/1 [HP>=7] -> GCTT>GC,GT PASS K=2;DP=12`
  - `chr20:499958 CA>C truth GT 0/1 [HP>=7] -> CAA>C PASS K=2;DP=23`
  - `chr20:1409097 G>GTGT truth GT 0/1 [other] -> ACTGC>A PASS K=2;DP=17`
  - `chr20:3609710 A>ATGTG truth GT 0/1 [TR<=50] -> A>ATGTG PASS K=2;DP=26`
- *genotype wrong (allele found, GT differs)*
  - `chr20:126062 C>CT truth GT 1/1 [HP>=7] -> C>CT LowProb K=1;DP=23`
  - `chr20:217716 ATG>A truth GT 0/1 [TR51-200] -> ATG>A PASS K=1;DP=16`
  - `chr20:240075 C>CAAA truth GT 1/1 [HP>=7] -> C>CAAA LowProb K=1;DP=13`
  - `chr20:566178 AAAAAAT>A truth GT 1/1 [HP>=7] -> AAAAAAT>A LowProb K=1;DP=18`
- *candidate found but filtered (LowProb/UNCERTAIN)*
  - `chr20:232033 A>AAAATAAAT,AAAATAAATAAAT truth GT 1/2 [TR<=50] -> A>AAAATAAATAAAT LowProb K=2;DP=17`
  - `chr20:290494 G>GC truth GT 1/1 [other] -> G>GC LowProb K=1;DP=10`
  - `chr20:365798 CTT>CT,C truth GT 2/1 [HP>=7] -> CTT>C LowProb K=1;DP=3`
  - `chr20:367246 AC>A truth GT 0/1 [TR<=50] -> AAAC>A LowProb K=2;DP=21`
- *multi-allelic site (truth 1/2), alleles genotyped separately*
  - `chr20:511703 TACAC>CACAC,T truth GT 2/1 [TR<=50] -> TACAC>T PASS K=1;DP=26`
  - `chr20:2314530 TACACAC>TACACACACAC,T truth GT 2/1 [TR<=50] -> TACACAC>T PASS K=2;DP=23`
  - `chr20:2647120 TTTTATTTA>TTTTATTTATTTATTTA,T truth GT 2/1 [TR<=50] -> TTTTATTTA>T,TTTTATTTATTTATTTA PASS K=2;DP=15`
  - `chr20:2754979 C>CA,G truth GT 2/1 [HP>=7] -> C>CA PASS K=1;DP=28`
- *never a candidate (discovery miss)*
  - `chr20:527001 C>CAGACAGAGAGAGAGAG truth GT 0/1 [TR51-200] -> no candidate`
  - `chr20:1404738 CTTT>C truth GT 0/1 [HP>=7] -> no candidate`
  - `chr20:1849802 AAAAT>A truth GT 0/1 [TR<=50] -> no candidate`
  - `chr20:1914466 AG>A truth GT 1/1 [other] -> no candidate`
- *allele/representation mismatch near a truth variant | HP>=7*
  - `chr20:92969 CTT>C,T called 2/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:499958 CAA>C called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:2523976 TCTC>T called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *genotype error (allele right, GT wrong) | TR51-200*
  - `chr20:217716 ATG>A called 1/1 truth 0/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:11498279 A>AGGAAGAAAGAGAGGGAGGAAAAGAGAGAGGAAGGGAGGG called 1/1 truth 0/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:13311243 T>TGATA called 1/1 truth 0/1 NALT=None NREF=None F=None HQ=None`
- *allele/representation mismatch near a truth variant | TR<=50*
  - `chr20:511703 TACAC>T called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:2314530 TACACAC>T called 1/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:2647120 TTTTATTTA>T called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *genotype error (allele right, GT wrong) | HP>=7*
  - `chr20:649701 CA>C called 1/1 truth 0/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:1122295 CTTTT>C called 0/1 truth 1/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:1138518 C>CA called 0/1 truth 1/1 NALT=None NREF=None F=None HQ=None`
- *allele/representation mismatch near a truth variant | TR51-200*
  - `chr20:3982103 GTC>G called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:3982105 C>CTCTCTCTCTC called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:4779167 CGTGTGTGT>C called 1/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *no truth variant nearby (pure FP) | other*
  - `chr20:4374492 CA>C called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:4374503 CTACCTCCAACCTT>C called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:22457785 G>GTTTTT called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *allele/representation mismatch near a truth variant | other*
  - `chr20:4619359 TC>T called 1/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:9051911 T>TA called 1/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:18410617 AC>A called 1/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *genotype error (allele right, GT wrong) | other*
  - `chr20:5206755 A>AC called 1/1 truth 0/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:10631606 AT>A called 1/1 truth 0/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:11561257 ATATAAGAATATCTCTTATATATGTACT>A called 0/1 truth 1/1 NALT=None NREF=None F=None HQ=None`
- *genotype error (allele right, GT wrong) | TR<=50*
  - `chr20:9553557 A>AGT called 1/1 truth 0/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:19356852 T>TTAATAATAATAATAATAA called 1/1 truth 0/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:40040821 G>GTCTC called 1/1 truth 0/1 NALT=None NREF=None F=None HQ=None`
- *no truth variant nearby (pure FP) | HP>=7*
  - `chr20:25926126 CA>C called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:26077020 GA>G called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`

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

## Error taxonomy: expA_multi_67x


**False negatives by cause**

| Cause | Count |
|---|---|
| different allele/representation called nearby | 159 |
| genotype wrong (allele found, GT differs) | 126 |
| multi-allelic site (truth 1/2), alleles genotyped separately | 113 |
| never a candidate (discovery miss) | 47 |

**False positives by cause and repeat context**

| Cause \| context | Count |
|---|---|
| allele/representation mismatch near a truth variant \| HP>=7 | 272 |
| genotype error (allele right, GT wrong) \| HP>=7 | 123 |
| allele/representation mismatch near a truth variant \| TR<=50 | 101 |
| allele/representation mismatch near a truth variant \| other | 43 |
| allele/representation mismatch near a truth variant \| TR51-200 | 38 |
| no truth variant nearby (pure FP) \| other | 34 |
| no truth variant nearby (pure FP) \| HP>=7 | 14 |
| genotype error (allele right, GT wrong) \| other | 7 |
| genotype error (allele right, GT wrong) \| TR51-200 | 5 |
| genotype error (allele right, GT wrong) \| TR<=50 | 4 |
| no truth variant nearby (pure FP) \| TR<=50 | 2 |
| no truth variant nearby (pure FP) \| TR51-200 | 1 |

**FP called genotype**

- allele/representation mismatch near a truth variant, called 0/1: 366
- genotype error, called 0/1: 128
- no truth variant nearby, called 0/1: 51
- allele/representation mismatch near a truth variant, called 1/1: 50
- allele/representation mismatch near a truth variant, called 2/1: 38
- genotype error, called 1/1: 11

**Recall: clustered (another truth variant within 10 bp) vs isolated**

- isolated: 10055/10352 = 0.9713
- clustered: 500/648 = 0.7716

**Median features: TP calls vs pure-FP calls**

| Feature | TP | pure FP |
|---|---|---|
| NALT | n/a | n/a |
| NREF | n/a | n/a |
| F | n/a | n/a |
| HQ | n/a | n/a |
| QMIN | n/a | n/a |

**Examples**

- *genotype wrong (allele found, GT differs)*
  - `chr20:217716 ATG>A truth GT 0/1 [TR51-200] -> ATG>A PASS K=1;DP=42`
  - `chr20:887636 CT>C truth GT 1/1 [HP>=7] -> CT>C PASS K=2;DP=38`
  - `chr20:1122295 CTTTT>C truth GT 1/1 [HP>=7] -> CTTTT>C PASS K=2;DP=45`
  - `chr20:1138518 C>CA truth GT 1/1 [HP>=7] -> C>CA PASS K=1;DP=24`
- *different allele/representation called nearby*
  - `chr20:240075 C>CAAA truth GT 1/1 [HP>=7] -> C>CAAA,CAA PASS K=2;DP=26`
  - `chr20:367246 AC>A truth GT 0/1 [TR<=50] -> AAAC>A PASS K=2;DP=52`
  - `chr20:367250 AAC>A truth GT 0/1 [TR<=50] -> AAAC>A PASS K=2;DP=52`
  - `chr20:499958 CA>C truth GT 0/1 [HP>=7] -> CAA>C PASS K=4;DP=73`
- *multi-allelic site (truth 1/2), alleles genotyped separately*
  - `chr20:511703 TACAC>CACAC,T truth GT 2/1 [TR<=50] -> TACAC>T PASS K=1;DP=68`
  - `chr20:942277 A>AT,T truth GT 1/2 [TR<=50] -> A>AT PASS K=3;DP=58`
  - `chr20:1903011 CAAAAAAAAA>CAAAAAAAAAAAAA,C truth GT 2/1 [HP>=7] -> CAAAAAAAAA>C PASS K=1;DP=40`
  - `chr20:2647120 TTTTATTTA>TTTTATTTATTTATTTA,T truth GT 2/1 [TR<=50] -> TTTTATTTA>T,TTTTATTTATTTATTTA PASS K=2;DP=40`
- *never a candidate (discovery miss)*
  - `chr20:1404738 CTTT>C truth GT 0/1 [HP>=7] -> no candidate`
  - `chr20:1914466 AG>A truth GT 1/1 [other] -> no candidate`
  - `chr20:1914707 CT>C truth GT 1/1 [other] -> no candidate`
  - `chr20:1915067 T>TAA truth GT 1/1 [other] -> no candidate`
- *genotype error (allele right, GT wrong) | TR51-200*
  - `chr20:217716 ATG>A called 1/1 truth 0/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:16801236 T>TTG called 0/1 truth 1/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:21798505 C>CTG called 0/1 truth 1/1 NALT=None NREF=None F=None HQ=None`
- *allele/representation mismatch near a truth variant | HP>=7*
  - `chr20:240075 C>CAA,CAAA called 2/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:458330 TA>T called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:499958 CAA>C called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *allele/representation mismatch near a truth variant | TR<=50*
  - `chr20:367249 AAAC>A called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:511703 TACAC>T called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:942277 A>AT called 1/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *genotype error (allele right, GT wrong) | HP>=7*
  - `chr20:887636 CT>C called 0/1 truth 1/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:1122295 CTTTT>C called 0/1 truth 1/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:1138518 C>CA called 0/1 truth 1/1 NALT=None NREF=None F=None HQ=None`
- *no truth variant nearby (pure FP) | HP>=7*
  - `chr20:3832562 CA>C called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:15931646 TA>T called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:25926126 CA>C called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *allele/representation mismatch near a truth variant | TR51-200*
  - `chr20:3982103 GTC>G called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:3982105 C>CTCTCTCTCTC called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:6124250 G>GCA,GCACA called 2/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *no truth variant nearby (pure FP) | other*
  - `chr20:4374492 CA>C called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:4374503 CTACCTCCAACCTT>C called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:7508395 G>GA called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *allele/representation mismatch near a truth variant | other*
  - `chr20:4513315 AG>A called 1/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:4619359 TC>T called 1/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:9051911 T>TA called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *genotype error (allele right, GT wrong) | other*
  - `chr20:9051906 G>GA called 0/1 truth 1/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:12759430 CG>C called 0/1 truth 1/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:22513143 TCTCTCTCC>T called 0/1 truth 1/1 NALT=None NREF=None F=None HQ=None`
- *genotype error (allele right, GT wrong) | TR<=50*
  - `chr20:20490861 A>AAC called 0/1 truth 1/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:40040821 G>GTCTC called 1/1 truth 0/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:49366207 TAC>T called 0/1 truth 1/1 NALT=None NREF=None F=None HQ=None`
- *no truth variant nearby (pure FP) | TR51-200*
  - `chr20:25895715 AT>A called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *no truth variant nearby (pure FP) | TR<=50*
  - `chr20:26079727 ACTAATCACTGGCTGCAACT>A called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:55234925 GAT>G called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`

## Error taxonomy: expB_str_67x


**False negatives by cause**

| Cause | Count |
|---|---|
| multi-allelic site (truth 1/2), alleles genotyped separately | 493 |
| different allele/representation called nearby | 215 |
| genotype wrong (allele found, GT differs) | 41 |
| never a candidate (discovery miss) | 31 |
| candidate found but filtered (LowProb/UNCERTAIN) | 4 |

**False positives by cause and repeat context**

| Cause \| context | Count |
|---|---|
| allele/representation mismatch near a truth variant \| HP>=7 | 758 |
| no truth variant nearby (pure FP) \| HP>=7 | 459 |
| allele/representation mismatch near a truth variant \| TR<=50 | 309 |
| no truth variant nearby (pure FP) \| TR<=50 | 211 |
| allele/representation mismatch near a truth variant \| other | 152 |
| allele/representation mismatch near a truth variant \| TR51-200 | 136 |
| no truth variant nearby (pure FP) \| TR51-200 | 89 |
| no truth variant nearby (pure FP) \| other | 56 |
| genotype error (allele right, GT wrong) \| HP>=7 | 25 |
| genotype error (allele right, GT wrong) \| TR51-200 | 10 |
| genotype error (allele right, GT wrong) \| TR<=50 | 9 |
| genotype error (allele right, GT wrong) \| other | 3 |

**FP called genotype**

- allele/representation mismatch near a truth variant, called 0/1: 598
- allele/representation mismatch near a truth variant, called 1/1: 577
- no truth variant nearby, called 1/1: 453
- no truth variant nearby, called 0/1: 342
- allele/representation mismatch near a truth variant, called 2/1: 180
- genotype error, called 1/1: 34
- no truth variant nearby, called 2/1: 20
- genotype error, called 0/1: 13

**Recall: clustered (another truth variant within 10 bp) vs isolated**

- isolated: 9678/10352 = 0.9349
- clustered: 538/648 = 0.8302

**Median features: TP calls vs pure-FP calls**

| Feature | TP | pure FP |
|---|---|---|
| NALT | 32.00 | 11.00 |
| NREF | 23.00 | 0.00 |
| F | 0.55 | 1.00 |
| HQ | 244.88 | 146.29 |
| QMIN | 0.91 | 0.80 |

Pure FPs with NALT <= 2: 168/815

**Examples**

- *genotype wrong (allele found, GT differs)*
  - `chr20:217716 ATG>A truth GT 0/1 [TR51-200] -> ATG>A PASS F=0.42857142857142855;NALT=18;NREF=24;NOTH=0;CR=0.30533254900671436;RP`
  - `chr20:649701 CA>C truth GT 0/1 [HP>=7] -> CA>C PASS F=0.5;NALT=19;NREF=19;NOTH=4;CR=0.44679838157198204;RPRI=0.22339919078`
  - `chr20:3106037 TG>T truth GT 0/1 [HP>=7] -> TG>T PASS F=0.5454545454545454;NALT=36;NREF=30;NOTH=1;CR=0.27816189992254053;RPR`
  - `chr20:3270136 C>CT truth GT 1/1 [HP>=7] -> C>CT PASS F=0.7368421052631579;NALT=14;NREF=5;NOTH=0;CR=0.5227490775817135;RPRI=`
- *multi-allelic site (truth 1/2), alleles genotyped separately*
  - `chr20:232033 A>AAAATAAAT,AAAATAAATAAAT truth GT 1/2 [TR<=50] -> A>AAAATAAAT PASS F=1;NALT=14;NREF=0;NOTH=22;CR=0.29503964290122536;RPRI=0.1475198214506; A>AAAATAAATAAAT PASS F=1;NALT=19;NREF=0;NOTH=17;CR=0.29503964290122536;RPRI=0.1475198214506`
  - `chr20:348484 CGTGTGT>CGTGT,C truth GT 2/1 [TR<=50] -> CGT>C PASS F=1;NALT=11;NREF=0;NOTH=14;CR=0.3275295386120749;RPRI=0.16376476930603; CGTGTGT>C PASS F=1;NALT=14;NREF=0;NOTH=11;CR=0.3275295386120749;RPRI=0.16376476930603`
  - `chr20:365798 CTT>CT,C truth GT 2/1 [HP>=7] -> CT>C PASS F=1;NALT=5;NREF=0;NOTH=7;CR=0.5293044287932317;RPRI=0.2646541877717401; CTT>C PASS F=1;NALT=6;NREF=0;NOTH=6;CR=0.5293044287932317;RPRI=0.2646522349573023`
  - `chr20:511703 TACAC>CACAC,T truth GT 2/1 [TR<=50] -> TACAC>T PASS F=0.5142857142857142;NALT=36;NREF=34;NOTH=1;CR=0.2010016732559959;RPRI`
- *different allele/representation called nearby*
  - `chr20:316724 GTT>G truth GT 0/1 [HP>=7] -> GT>G PASS F=0.08;NALT=2;NREF=23;NOTH=15;CR=0.4587943148124999;RPRI=0.22939719495; GTT>G PASS F=0.3783783783783784;NALT=14;NREF=23;NOTH=3;CR=0.45264214230694844;RPR`
  - `chr20:367246 AC>A truth GT 0/1 [TR<=50] -> AAAC>A LowProb F=0.48936170212765956;NALT=23;NREF=24;NOTH=3;CR=0.11339129895466793;RP; AAAC>A PASS F=0.10344827586206896;NALT=3;NREF=26;NOTH=23;CR=0.10326292289633747;RP`
  - `chr20:367250 AAC>A truth GT 0/1 [TR<=50] -> AAAC>A LowProb F=0.48936170212765956;NALT=23;NREF=24;NOTH=3;CR=0.11339129895466793;RP; AAAC>A PASS F=0.10344827586206896;NALT=3;NREF=26;NOTH=23;CR=0.10326292289633747;RP`
  - `chr20:499958 CA>C truth GT 0/1 [HP>=7] -> AG>A LowProb F=0.05263157894736842;NALT=2;NREF=36;NOTH=36;CR=0.13566070966823549;RP; GC>G LowProb F=0.05405405405405406;NALT=2;NREF=35;NOTH=36;CR=0.13760917518252758;RP`
- *never a candidate (discovery miss)*
  - `chr20:1404738 CTTT>C truth GT 0/1 [HP>=7] -> no candidate`
  - `chr20:1914466 AG>A truth GT 1/1 [other] -> no candidate`
  - `chr20:1914707 CT>C truth GT 1/1 [other] -> no candidate`
  - `chr20:1915067 T>TAA truth GT 1/1 [other] -> no candidate`
- *candidate found but filtered (LowProb/UNCERTAIN)*
  - `chr20:32208648 CAAA>CAAAA,C truth GT 2/1 [HP>=7] -> C>CA LowProb F=0.5;NALT=2;NREF=2;NOTH=3;CR=0.5563556539272072;RPRI=0.33816829806178`
  - `chr20:53745828 C>CA truth GT 0/1 [HP>=7] -> AC>A LowProb F=0.031746031746031744;NALT=2;NREF=61;NOTH=1;CR=0.1683116445783832;RPR`
  - `chr20:59524356 T>TA truth GT 0/1 [HP>=7] -> T>TAA LowProb F=0.0625;NALT=2;NREF=30;NOTH=23;CR=0.3769647309549455;RPRI=0.188482365`
  - `chr20:59524369 G>GA truth GT 0/1 [HP>=7] -> A>AAG LowProb F=0.40384615384615385;NALT=21;NREF=31;NOTH=1;CR=0.3066571180728725;RPR`
- *no truth variant nearby (pure FP) | HP>=7*
  - `chr20:92967 TGC>T called 0/1 truth ./. NALT=2 NREF=10 F=0.16666666666666666 HQ=10.793388486883373`
  - `chr20:92968 GCT>G called 0/1 truth ./. NALT=3 NREF=9 F=0.25 HQ=24.385772036459684`
  - `chr20:152450 CAAAA>C called 0/1 truth ./. NALT=2 NREF=0 F=1 HQ=107.96599371151706`
- *genotype error (allele right, GT wrong) | TR51-200*
  - `chr20:217716 ATG>A called 1/1 truth 0/1 NALT=18 NREF=24 F=0.42857142857142855 HQ=143.58869919142944`
  - `chr20:16801236 T>TTG called 0/1 truth 1/1 NALT=10 NREF=4 F=0.7142857142857143 HQ=35.78872643764328`
  - `chr20:18488953 A>ATG called 1/1 truth 0/1 NALT=30 NREF=2 F=0.9375 HQ=276.9337815447602`
- *allele/representation mismatch near a truth variant | TR<=50*
  - `chr20:232033 A>AAAATAAAT called 1/1 truth ./. NALT=14 NREF=0 F=1 HQ=497.9742852082418`
  - `chr20:232033 A>AAAATAAATAAAT called 1/1 truth ./. NALT=19 NREF=0 F=1 HQ=406.6161792333198`
  - `chr20:367249 AAAC>A called 0/1 truth ./. NALT=3 NREF=26 F=0.10344827586206896 HQ=71.70705411297877`
- *allele/representation mismatch near a truth variant | HP>=7*
  - `chr20:240075 C>CAA called 1/1 truth ./. NALT=4 NREF=0 F=1 HQ=122.4707730810319`
  - `chr20:316724 GTT>GT,G called 2/1 truth ./. NALT=2 NREF=23 F=0.08 HQ=13.437775349710591`
  - `chr20:368302 TTG>T called 0/1 truth ./. NALT=4 NREF=31 F=0.11428571428571428 HQ=15.977844485590595`
- *no truth variant nearby (pure FP) | TR<=50*
  - `chr20:348484 CGT>C called 1/1 truth ./. NALT=11 NREF=0 F=1 HQ=99.06368319226885`
  - `chr20:348484 CGTGTGT>C called 0/1 truth ./. NALT=14 NREF=0 F=1 HQ=153.73056339773768`
  - `chr20:685696 AAATAATAAT>A called 1/1 truth ./. NALT=27 NREF=0 F=1 HQ=822.7399334601173`
- *allele/representation mismatch near a truth variant | other*
  - `chr20:611880 GA>G called 0/1 truth ./. NALT=2 NREF=60 F=0.03225806451612903 HQ=11.927025091363475`
  - `chr20:1012670 A>AG called 0/1 truth ./. NALT=2 NREF=45 F=0.0425531914893617 HQ=16.45180216478218`
  - `chr20:1362011 TG>T called 0/1 truth ./. NALT=2 NREF=55 F=0.03508771929824561 HQ=15.359432219768223`
- *genotype error (allele right, GT wrong) | HP>=7*
  - `chr20:649701 CA>C called 1/1 truth 0/1 NALT=19 NREF=19 F=0.5 HQ=120.98339770637244`
  - `chr20:3106037 TG>T called 1/1 truth 0/1 NALT=36 NREF=30 F=0.5454545454545454 HQ=250.64054890320062`
  - `chr20:3270136 C>CT called 0/1 truth 1/1 NALT=14 NREF=5 F=0.7368421052631579 HQ=37.532543715236216`
- *allele/representation mismatch near a truth variant | TR51-200*
  - `chr20:672388 CGT>C called 0/1 truth ./. NALT=2 NREF=29 F=0.06451612903225806 HQ=12.491642067620091`
  - `chr20:2249002 CTT>C called 0/1 truth ./. NALT=2 NREF=30 F=0.0625 HQ=249.8805514867372`
  - `chr20:2249004 T>TTC called 0/1 truth ./. NALT=2 NREF=31 F=0.06060606060606061 HQ=42.829649029331904`
- *no truth variant nearby (pure FP) | TR51-200*
  - `chr20:2011231 CCTCTCTCT>C called 1/1 truth ./. NALT=25 NREF=0 F=1 HQ=907.506380483889`
  - `chr20:2011231 CCTCTCTCTCTCTCTCT>C called 1/1 truth ./. NALT=31 NREF=0 F=1 HQ=522.3248595285688`
  - `chr20:4779162 TTGTGCG>T called 0/1 truth ./. NALT=2 NREF=13 F=0.13333333333333333 HQ=189.81944752449022`
- *no truth variant nearby (pure FP) | other*
  - `chr20:4374492 CA>C called 0/1 truth ./. NALT=7 NREF=42 F=0.14285714285714285 HQ=44.54178567920849`
  - `chr20:4374503 CTACCTCCAACCTT>C called 0/1 truth ./. NALT=3 NREF=36 F=0.07692307692307693 HQ=36.74718330053361`
  - `chr20:7508395 G>GA called 0/1 truth ./. NALT=13 NREF=66 F=0.16455696202531644 HQ=54.65798647948225`
- *genotype error (allele right, GT wrong) | TR<=50*
  - `chr20:5454164 AATATAT>A called 1/1 truth 0/1 NALT=23 NREF=1 F=0.9583333333333334 HQ=449.3046759699655`
  - `chr20:5454204 TATA>T called 1/1 truth 0/1 NALT=20 NREF=49 F=0.2898550724637681 HQ=377.4577162470696`
  - `chr20:9553557 A>AGT called 1/1 truth 0/1 NALT=31 NREF=29 F=0.5166666666666667 HQ=233.49482890086074`
- *genotype error (allele right, GT wrong) | other*
  - `chr20:19711132 A>AAC called 0/1 truth 1/1 NALT=41 NREF=1 F=0.9761904761904762 HQ=105.99346769347812`
  - `chr20:62420350 GCA>G called 1/1 truth 0/1 NALT=19 NREF=25 F=0.4318181818181818 HQ=88.91922850655078`
  - `chr20:63701286 A>AC called 1/1 truth 0/1 NALT=35 NREF=41 F=0.4605263157894737 HQ=304.60559228727266`

## Error taxonomy: expC_frag_67x


**False negatives by cause**

| Cause | Count |
|---|---|
| multi-allelic site (truth 1/2), alleles genotyped separately | 438 |
| genotype wrong (allele found, GT differs) | 208 |
| different allele/representation called nearby | 168 |
| candidate found but filtered (LowProb/UNCERTAIN) | 70 |
| never a candidate (discovery miss) | 33 |

**False positives by cause and repeat context**

| Cause \| context | Count |
|---|---|
| allele/representation mismatch near a truth variant \| HP>=7 | 793 |
| no truth variant nearby (pure FP) \| HP>=7 | 420 |
| allele/representation mismatch near a truth variant \| TR<=50 | 289 |
| no truth variant nearby (pure FP) \| TR<=50 | 206 |
| allele/representation mismatch near a truth variant \| TR51-200 | 113 |
| no truth variant nearby (pure FP) \| TR51-200 | 80 |
| genotype error (allele right, GT wrong) \| HP>=7 | 77 |
| allele/representation mismatch near a truth variant \| other | 70 |
| no truth variant nearby (pure FP) \| other | 50 |
| genotype error (allele right, GT wrong) \| TR<=50 | 13 |
| genotype error (allele right, GT wrong) \| TR51-200 | 10 |
| genotype error (allele right, GT wrong) \| other | 3 |

**FP called genotype**

- allele/representation mismatch near a truth variant, called 0/1: 690
- allele/representation mismatch near a truth variant, called 1/1: 447
- no truth variant nearby, called 1/1: 403
- no truth variant nearby, called 0/1: 341
- allele/representation mismatch near a truth variant, called 2/1: 128
- genotype error, called 0/1: 71
- genotype error, called 1/1: 32
- no truth variant nearby, called 2/1: 12

**Recall: clustered (another truth variant within 10 bp) vs isolated**

- isolated: 9552/10352 = 0.9227
- clustered: 531/648 = 0.8194

**Median features: TP calls vs pure-FP calls**

| Feature | TP | pure FP |
|---|---|---|
| NALT | 28.00 | 11.00 |
| NREF | 20.00 | 0.00 |
| F | 0.54 | 1.00 |
| HQ | 245.46 | 132.15 |
| QMIN | 0.92 | 0.79 |

Pure FPs with NALT <= 2: 213/756

**Examples**

- *genotype wrong (allele found, GT differs)*
  - `chr20:92969 CTT>C truth GT 1/1 [HP>=7] -> CTT>C LowProb F=1;NALT=15;NREF=0;NOTH=4;CR=0.52469796468502;RPRI=0.26234898234251;TA`
  - `chr20:217716 ATG>A truth GT 0/1 [TR51-200] -> ATG>A PASS F=0.4444444444444444;NALT=16;NREF=20;NOTH=0;CR=0.31128492995909535;RPR`
  - `chr20:240075 C>CAAA truth GT 1/1 [HP>=7] -> C>CAA LowProb F=1;NALT=4;NREF=0;NOTH=20;CR=0.49816965904673927;RPRI=0.24926902618476; C>CAAA PASS F=1;NALT=19;NREF=0;NOTH=5;CR=0.49816965904673927;RPRI=0.24908482952336`
  - `chr20:887636 CT>C truth GT 1/1 [HP>=7] -> CT>C LowProb F=0.75;NALT=6;NREF=2;NOTH=8;CR=0.46569580501495844;RPRI=0.232850381680`
- *multi-allelic site (truth 1/2), alleles genotyped separately*
  - `chr20:232033 A>AAAATAAAT,AAAATAAATAAAT truth GT 1/2 [TR<=50] -> A>AAAATAAAT PASS F=1;NALT=13;NREF=0;NOTH=21;CR=0.2954481396332515;RPRI=0.14772406981662; A>AAAATAAATAAAT PASS F=1;NALT=18;NREF=0;NOTH=16;CR=0.2954481396332515;RPRI=0.14772406981662`
  - `chr20:348484 CGTGTGT>CGTGT,C truth GT 2/1 [TR<=50] -> CGT>C PASS F=1;NALT=11;NREF=0;NOTH=12;CR=0.33361649513381403;RPRI=0.1668082475669; CGTGTGT>C PASS F=1;NALT=12;NREF=0;NOTH=11;CR=0.33361649513381403;RPRI=0.1668082475669`
  - `chr20:511703 TACAC>CACAC,T truth GT 2/1 [TR<=50] -> TACAC>T PASS F=0.55;NALT=33;NREF=27;NOTH=1;CR=0.19626835531556644;RPRI=0.0981341776`
  - `chr20:571051 C>CT,CTTT truth GT 1/2 [HP>=7] -> C>CT PASS F=0.9285714285714286;NALT=13;NREF=1;NOTH=9;CR=0.534995218360776;RPRI=0; C>CTTT PASS F=0.8571428571428571;NALT=6;NREF=1;NOTH=16;CR=0.534995218360776;RPRI=0`
- *candidate found but filtered (LowProb/UNCERTAIN)*
  - `chr20:290494 G>GC truth GT 1/1 [other] -> G>GC LowProb F=1;NALT=23;NREF=0;NOTH=0;CR=0.1418872205440661;RPRI=0.070943610272033`
  - `chr20:787512 CT>C truth GT 1/1 [HP>=7] -> CT>C LowProb F=1;NALT=14;NREF=0;NOTH=6;CR=0.4838248403179828;RPRI=0.241912420158991`
  - `chr20:1311900 C>CTT,CTTTTTT truth GT 1/2 [HP>=7] -> C>CT LowProb F=1;NALT=2;NREF=0;NOTH=20;CR=0.4650257196527999;RPRI=0.657816730404831; C>CTT LowProb F=1;NALT=12;NREF=0;NOTH=10;CR=0.4650257196527999;RPRI=0.23251285982639`
  - `chr20:1882863 CTTTTTTTTT>C truth GT 0/1 [HP>=7] -> CTTTTTTTTT>C LowProb F=0.76;NALT=19;NREF=6;NOTH=4;CR=0.4692768152441251;RPRI=0.234638407622; CTTTTTTTTTT>C LowProb F=0.25;NALT=2;NREF=6;NOTH=21;CR=0.46627051276420733;RPRI=0.23726716477`
- *different allele/representation called nearby*
  - `chr20:316724 GTT>G truth GT 0/1 [HP>=7] -> GT>G PASS F=0.10526315789473684;NALT=2;NREF=17;NOTH=14;CR=0.43190037541856047;RP; GTT>G PASS F=0.43333333333333335;NALT=13;NREF=17;NOTH=3;CR=0.425748202913009;RPRI`
  - `chr20:367246 AC>A truth GT 0/1 [TR<=50] -> AAAC>A PASS F=0.5;NALT=20;NREF=20;NOTH=3;CR=0.11420525244304001;RPRI=0.05710262622; AAAC>A PASS F=0.13043478260869565;NALT=3;NREF=20;NOTH=20;CR=0.10426918407701725;RP`
  - `chr20:367250 AAC>A truth GT 0/1 [TR<=50] -> AAAC>A PASS F=0.5;NALT=20;NREF=20;NOTH=3;CR=0.11420525244304001;RPRI=0.05710262622; AAAC>A PASS F=0.13043478260869565;NALT=3;NREF=20;NOTH=20;CR=0.10426918407701725;RP`
  - `chr20:499958 CA>C truth GT 0/1 [HP>=7] -> CA>C PASS F=0.5294117647058824;NALT=27;NREF=24;NOTH=6;CR=0.4413646313516866;RPRI; CAA>C PASS F=0.1111111111111111;NALT=3;NREF=24;NOTH=30;CR=0.43444308948477983;RPR`
- *never a candidate (discovery miss)*
  - `chr20:1404738 CTTT>C truth GT 0/1 [HP>=7] -> no candidate`
  - `chr20:1914466 AG>A truth GT 1/1 [other] -> no candidate`
  - `chr20:1914707 CT>C truth GT 1/1 [other] -> no candidate`
  - `chr20:1915067 T>TAA truth GT 1/1 [other] -> no candidate`
- *no truth variant nearby (pure FP) | HP>=7*
  - `chr20:152447 TAGC>T called 0/1 truth ./. NALT=3 NREF=11 F=0.21428571428571427 HQ=36.21712055353919`
  - `chr20:152450 CAAAA>C called 0/1 truth ./. NALT=2 NREF=0 F=1 HQ=110.6358534446531`
  - `chr20:736736 CA>C called 1/1 truth ./. NALT=26 NREF=0 F=1 HQ=279.79366989857317`
- *genotype error (allele right, GT wrong) | TR51-200*
  - `chr20:217716 ATG>A called 1/1 truth 0/1 NALT=16 NREF=20 F=0.4444444444444444 HQ=135.32607933926352`
  - `chr20:10507904 T>TTA called 1/1 truth 0/1 NALT=27 NREF=1 F=0.9642857142857143 HQ=304.53837759961584`
  - `chr20:18488953 A>ATG called 1/1 truth 0/1 NALT=27 NREF=2 F=0.9310344827586207 HQ=363.0805487784121`
- *allele/representation mismatch near a truth variant | TR<=50*
  - `chr20:232033 A>AAAATAAAT called 1/1 truth ./. NALT=13 NREF=0 F=1 HQ=143.37122807362243`
  - `chr20:232033 A>AAAATAAATAAAT called 1/1 truth ./. NALT=18 NREF=0 F=1 HQ=377.16074594724284`
  - `chr20:367244 AAAC>A called 0/1 truth ./. NALT=20 NREF=20 F=0.5 HQ=68.99391527042974`
- *genotype error (allele right, GT wrong) | HP>=7*
  - `chr20:240075 C>CAAA called 0/1 truth 1/1 NALT=19 NREF=0 F=1 HQ=12.605009281055423`
  - `chr20:1065654 C>CT called 0/1 truth 1/1 NALT=34 NREF=3 F=0.918918918918919 HQ=22.446423445923756`
  - `chr20:1091299 C>CT called 0/1 truth 1/1 NALT=22 NREF=1 F=0.9565217391304348 HQ=20.319526032692536`
- *allele/representation mismatch near a truth variant | HP>=7*
  - `chr20:316724 GTT>GT,G called 2/1 truth ./. NALT=2 NREF=17 F=0.10526315789473684 HQ=16.02862004237943`
  - `chr20:368302 TTG>T called 0/1 truth ./. NALT=2 NREF=28 F=0.06666666666666667 HQ=11.506532197444812`
  - `chr20:499958 CAA>CA,C called 2/1 truth ./. NALT=27 NREF=24 F=0.5294117647058824 HQ=172.2451325565362`
- *no truth variant nearby (pure FP) | TR<=50*
  - `chr20:348484 CGT>C called 1/1 truth ./. NALT=11 NREF=0 F=1 HQ=136.38214241367393`
  - `chr20:348484 CGTGTGT>C called 0/1 truth ./. NALT=12 NREF=0 F=1 HQ=57.214760165256564`
  - `chr20:685696 AAATAATAAT>A called 1/1 truth ./. NALT=23 NREF=0 F=1 HQ=418.9509191239294`
- *no truth variant nearby (pure FP) | TR51-200*
  - `chr20:2011231 CCTCTCTCT>C called 1/1 truth ./. NALT=23 NREF=0 F=1 HQ=755.3464759755009`
  - `chr20:2011231 CCTCTCTCTCTCTCTCT>C called 1/1 truth ./. NALT=25 NREF=0 F=1 HQ=240.83885796308186`
  - `chr20:4779167 CGTGTGT>C called 1/1 truth ./. NALT=12 NREF=0 F=1 HQ=379.23817329261783`
- *genotype error (allele right, GT wrong) | other*
  - `chr20:2499288 AC>A called 1/1 truth 0/1 NALT=18 NREF=13 F=0.5806451612903226 HQ=14.842942556216784`
  - `chr20:62491856 C>CCCCCACCCCAGAGTGCCCTCTCTGCTCTCCCAGCT called 0/1 truth 1/1 NALT=24 NREF=0 F=1 HQ=16.52356079239982`
  - `chr20:63701286 A>AC called 1/1 truth 0/1 NALT=27 NREF=33 F=0.45 HQ=233.61542165158153`
- *no truth variant nearby (pure FP) | other*
  - `chr20:2783142 G>GTTTTTCA called 0/1 truth ./. NALT=2 NREF=27 F=0.06896551724137931 HQ=12.826983344258224`
  - `chr20:2783152 CGATTAGG>C called 0/1 truth ./. NALT=2 NREF=29 F=0.06451612903225806 HQ=12.131570061498792`
  - `chr20:4374492 CA>C called 0/1 truth ./. NALT=7 NREF=36 F=0.16279069767441862 HQ=50.086727720491126`
- *allele/representation mismatch near a truth variant | other*
  - `chr20:3188828 C>CTTT called 1/1 truth ./. NALT=4 NREF=0 F=1 HQ=241.56108957520337`
  - `chr20:4125629 A>AC called 1/1 truth ./. NALT=2 NREF=0 F=1 HQ=56.74491845533703`
  - `chr20:4619359 TC>T called 1/1 truth ./. NALT=2 NREF=57 F=0.03389830508474576 HQ=525.0988926408149`
- *allele/representation mismatch near a truth variant | TR51-200*
  - `chr20:3310333 A>AAG called 0/1 truth ./. NALT=3 NREF=13 F=0.1875 HQ=91.974890453645`
  - `chr20:3310333 A>AAAAGAAAG called 1/1 truth ./. NALT=12 NREF=0 F=1 HQ=525.7516535185596`
  - `chr20:3310333 A>AAAAGAAAGAAAGAAAG called 1/1 truth ./. NALT=19 NREF=0 F=1 HQ=267.73353926069706`
- *genotype error (allele right, GT wrong) | TR<=50*
  - `chr20:5454164 AATATAT>A called 1/1 truth 0/1 NALT=21 NREF=1 F=0.9545454545454546 HQ=335.9276328654146`
  - `chr20:9553557 A>AGT called 1/1 truth 0/1 NALT=28 NREF=23 F=0.5490196078431373 HQ=262.0270840772897`
  - `chr20:10476142 C>CTG called 1/1 truth 0/1 NALT=25 NREF=0 F=1 HQ=298.2444361806263`

## Error taxonomy: expC2_fragend_67x


**False negatives by cause**

| Cause | Count |
|---|---|
| multi-allelic site (truth 1/2), alleles genotyped separately | 446 |
| different allele/representation called nearby | 182 |
| genotype wrong (allele found, GT differs) | 179 |
| candidate found but filtered (LowProb/UNCERTAIN) | 59 |
| never a candidate (discovery miss) | 32 |

**False positives by cause and repeat context**

| Cause \| context | Count |
|---|---|
| allele/representation mismatch near a truth variant \| HP>=7 | 943 |
| no truth variant nearby (pure FP) \| HP>=7 | 472 |
| allele/representation mismatch near a truth variant \| TR<=50 | 300 |
| no truth variant nearby (pure FP) \| TR<=50 | 208 |
| allele/representation mismatch near a truth variant \| TR51-200 | 121 |
| allele/representation mismatch near a truth variant \| other | 84 |
| no truth variant nearby (pure FP) \| TR51-200 | 83 |
| no truth variant nearby (pure FP) \| other | 53 |
| genotype error (allele right, GT wrong) \| HP>=7 | 52 |
| genotype error (allele right, GT wrong) \| TR<=50 | 15 |
| genotype error (allele right, GT wrong) \| TR51-200 | 9 |
| genotype error (allele right, GT wrong) \| other | 4 |

**FP called genotype**

- allele/representation mismatch near a truth variant, called 0/1: 848
- allele/representation mismatch near a truth variant, called 1/1: 462
- no truth variant nearby, called 1/1: 413
- no truth variant nearby, called 0/1: 390
- allele/representation mismatch near a truth variant, called 2/1: 138
- genotype error, called 0/1: 48
- genotype error, called 1/1: 32
- no truth variant nearby, called 2/1: 13

**Recall: clustered (another truth variant within 10 bp) vs isolated**

- isolated: 9571/10352 = 0.9246
- clustered: 531/648 = 0.8194

**Median features: TP calls vs pure-FP calls**

| Feature | TP | pure FP |
|---|---|---|
| NALT | 28.00 | 8.00 |
| NREF | 20.00 | 0.00 |
| F | 0.55 | 1.00 |
| HQ | 254.64 | 137.58 |
| QMIN | 0.92 | 0.78 |

Pure FPs with NALT <= 2: 229/816

**Examples**

- *genotype wrong (allele found, GT differs)*
  - `chr20:217716 ATG>A truth GT 0/1 [TR51-200] -> ATG>A PASS F=0.4444444444444444;NALT=16;NREF=20;NOTH=0;CR=0.31128492995909535;RPR`
  - `chr20:240075 C>CAAA truth GT 1/1 [HP>=7] -> C>CAA LowProb F=1;NALT=4;NREF=0;NOTH=20;CR=0.49816965904673927;RPRI=0.24926902618476; C>CAAA PASS F=1;NALT=19;NREF=0;NOTH=5;CR=0.49816965904673927;RPRI=0.24908482952336`
  - `chr20:887636 CT>C truth GT 1/1 [HP>=7] -> CT>C LowProb F=0.75;NALT=6;NREF=2;NOTH=8;CR=0.46569580501495844;RPRI=0.232850381680`
  - `chr20:1138518 C>CA truth GT 1/1 [HP>=7] -> C>CA PASS F=0.7894736842105263;NALT=15;NREF=4;NOTH=3;CR=0.5223084263363668;RPRI=`
- *multi-allelic site (truth 1/2), alleles genotyped separately*
  - `chr20:232033 A>AAAATAAAT,AAAATAAATAAAT truth GT 1/2 [TR<=50] -> A>AAAATAAAT PASS F=1;NALT=13;NREF=0;NOTH=21;CR=0.2954481396332515;RPRI=0.14772406981662; A>AAAATAAATAAAT PASS F=1;NALT=18;NREF=0;NOTH=16;CR=0.2954481396332515;RPRI=0.14772406981662`
  - `chr20:348484 CGTGTGT>CGTGT,C truth GT 2/1 [TR<=50] -> CGT>C PASS F=1;NALT=11;NREF=0;NOTH=12;CR=0.33361649513381403;RPRI=0.1668082475669; CGTGTGT>C PASS F=1;NALT=12;NREF=0;NOTH=11;CR=0.33361649513381403;RPRI=0.1668082475669`
  - `chr20:365798 CTT>CT,C truth GT 2/1 [HP>=7] -> CT>C LowProb F=1;NALT=5;NREF=0;NOTH=6;CR=0.516046853035656;RPRI=0.2580253998929523;; CTT>C PASS F=1;NALT=5;NREF=0;NOTH=6;CR=0.5387741257629287;RPRI=0.2693890362565886`
  - `chr20:511703 TACAC>CACAC,T truth GT 2/1 [TR<=50] -> TACAC>T PASS F=0.55;NALT=33;NREF=27;NOTH=1;CR=0.19626835531556644;RPRI=0.0981341776`
- *candidate found but filtered (LowProb/UNCERTAIN)*
  - `chr20:290494 G>GC truth GT 1/1 [other] -> G>GC LowProb F=1;NALT=23;NREF=0;NOTH=0;CR=0.1418872205440661;RPRI=0.070943610272033`
  - `chr20:787512 CT>C truth GT 1/1 [HP>=7] -> CT>C LowProb F=1;NALT=14;NREF=0;NOTH=6;CR=0.4838248403179828;RPRI=0.241912420158991`
  - `chr20:1311900 C>CTT,CTTTTTT truth GT 1/2 [HP>=7] -> C>CT LowProb F=1;NALT=2;NREF=0;NOTH=20;CR=0.4650257196527999;RPRI=0.657816730404831; C>CTT LowProb F=1;NALT=12;NREF=0;NOTH=10;CR=0.4650257196527999;RPRI=0.23251285982639`
  - `chr20:1882863 CTTTTTTTTT>C truth GT 0/1 [HP>=7] -> CTTTTTTTTT>C LowProb F=0.76;NALT=19;NREF=6;NOTH=4;CR=0.4692768152441251;RPRI=0.234638407622; CTTTTTTTTTT>C LowProb F=0.25;NALT=2;NREF=6;NOTH=21;CR=0.46627051276420733;RPRI=0.23726716477`
- *different allele/representation called nearby*
  - `chr20:316724 GTT>G truth GT 0/1 [HP>=7] -> GT>G PASS F=0.10526315789473684;NALT=2;NREF=17;NOTH=14;CR=0.46977916329734837;RP; GTT>G PASS F=0.43333333333333335;NALT=13;NREF=17;NOTH=3;CR=0.4636269907917969;RPR`
  - `chr20:367246 AC>A truth GT 0/1 [TR<=50] -> AAAC>A PASS F=0.5;NALT=20;NREF=20;NOTH=3;CR=0.11420525244304001;RPRI=0.05710262622; AAAC>A PASS F=0.13043478260869565;NALT=3;NREF=20;NOTH=20;CR=0.09845523058864517;RP`
  - `chr20:367250 AAC>A truth GT 0/1 [TR<=50] -> AAAC>A PASS F=0.5;NALT=20;NREF=20;NOTH=3;CR=0.11420525244304001;RPRI=0.05710262622; AAAC>A PASS F=0.13043478260869565;NALT=3;NREF=20;NOTH=20;CR=0.09845523058864517;RP`
  - `chr20:499958 CA>C truth GT 0/1 [HP>=7] -> CA>C PASS F=0.5294117647058824;NALT=27;NREF=24;NOTH=6;CR=0.4413646313516866;RPRI; CAA>C PASS F=0.1111111111111111;NALT=3;NREF=24;NOTH=30;CR=0.43444308948477983;RPR`
- *never a candidate (discovery miss)*
  - `chr20:1404738 CTTT>C truth GT 0/1 [HP>=7] -> no candidate`
  - `chr20:1914466 AG>A truth GT 1/1 [other] -> no candidate`
  - `chr20:1914707 CT>C truth GT 1/1 [other] -> no candidate`
  - `chr20:1915067 T>TAA truth GT 1/1 [other] -> no candidate`
- *no truth variant nearby (pure FP) | HP>=7*
  - `chr20:92968 GCT>G called 0/1 truth ./. NALT=3 NREF=7 F=0.3 HQ=67.55891322741388`
  - `chr20:152447 TAGC>T called 0/1 truth ./. NALT=3 NREF=9 F=0.25 HQ=36.44825378077667`
  - `chr20:152450 CAAAA>C called 0/1 truth ./. NALT=2 NREF=0 F=1 HQ=110.6358534446531`
- *genotype error (allele right, GT wrong) | TR51-200*
  - `chr20:217716 ATG>A called 1/1 truth 0/1 NALT=16 NREF=20 F=0.4444444444444444 HQ=136.69607024406173`
  - `chr20:10507904 T>TTA called 1/1 truth 0/1 NALT=26 NREF=1 F=0.9629629629629629 HQ=302.226094162942`
  - `chr20:18488953 A>ATG called 1/1 truth 0/1 NALT=27 NREF=2 F=0.9310344827586207 HQ=376.34176178300936`
- *allele/representation mismatch near a truth variant | TR<=50*
  - `chr20:232033 A>AAAATAAAT called 1/1 truth ./. NALT=13 NREF=0 F=1 HQ=143.37122807362243`
  - `chr20:232033 A>AAAATAAATAAAT called 1/1 truth ./. NALT=18 NREF=0 F=1 HQ=107.15009321530462`
  - `chr20:319668 CAA>C called 0/1 truth ./. NALT=2 NREF=40 F=0.047619047619047616 HQ=48.4273014676337`
- *genotype error (allele right, GT wrong) | HP>=7*
  - `chr20:240075 C>CAAA called 0/1 truth 1/1 NALT=19 NREF=0 F=1 HQ=12.605009281055423`
  - `chr20:1138518 C>CA called 0/1 truth 1/1 NALT=15 NREF=4 F=0.7894736842105263 HQ=102.5496330123533`
  - `chr20:2542038 G>GA called 0/1 truth 1/1 NALT=47 NREF=1 F=0.9791666666666666 HQ=254.05347156122116`
- *allele/representation mismatch near a truth variant | HP>=7*
  - `chr20:316724 GTT>GT,G called 2/1 truth ./. NALT=2 NREF=17 F=0.10526315789473684 HQ=24.702795347051733`
  - `chr20:368302 TTG>T called 0/1 truth ./. NALT=3 NREF=27 F=0.1 HQ=45.40994806501885`
  - `chr20:425625 TC>T called 0/1 truth ./. NALT=2 NREF=17 F=0.10526315789473684 HQ=11.925688415713921`
- *no truth variant nearby (pure FP) | TR<=50*
  - `chr20:348484 CGT>C called 1/1 truth ./. NALT=11 NREF=0 F=1 HQ=72.48907745625581`
  - `chr20:348484 CGTGTGT>C called 0/1 truth ./. NALT=12 NREF=0 F=1 HQ=177.8318861559584`
  - `chr20:685696 AAATAATAAT>A called 1/1 truth ./. NALT=23 NREF=0 F=1 HQ=418.9509191239294`
- *no truth variant nearby (pure FP) | TR51-200*
  - `chr20:2011231 CCTCTCTCT>C called 1/1 truth ./. NALT=23 NREF=0 F=1 HQ=755.3464775552713`
  - `chr20:2011231 CCTCTCTCTCTCTCTCT>C called 1/1 truth ./. NALT=25 NREF=0 F=1 HQ=386.0817850424023`
  - `chr20:4779167 CGTGTGT>C called 1/1 truth ./. NALT=12 NREF=0 F=1 HQ=379.23817329261783`
- *genotype error (allele right, GT wrong) | other*
  - `chr20:2499288 AC>A called 1/1 truth 0/1 NALT=18 NREF=13 F=0.5806451612903226 HQ=14.247430174854799`
  - `chr20:22513143 TCTCTCTCC>T called 0/1 truth 1/1 NALT=11 NREF=22 F=0.3333333333333333 HQ=239.88705038687056`
  - `chr20:62491856 C>CCCCCACCCCAGAGTGCCCTCTCTGCTCTCCCAGCT called 0/1 truth 1/1 NALT=24 NREF=0 F=1 HQ=16.52356079239982`
- *no truth variant nearby (pure FP) | other*
  - `chr20:2783142 G>GTTTTTCA called 0/1 truth ./. NALT=2 NREF=27 F=0.06896551724137931 HQ=12.826983344258224`
  - `chr20:2783152 CGATTAGG>C called 0/1 truth ./. NALT=2 NREF=29 F=0.06451612903225806 HQ=12.131570061498792`
  - `chr20:4374492 CA>C called 0/1 truth ./. NALT=7 NREF=36 F=0.16279069767441862 HQ=50.086727723503145`
- *allele/representation mismatch near a truth variant | other*
  - `chr20:3188828 C>CTTT called 1/1 truth ./. NALT=4 NREF=0 F=1 HQ=241.56108957520337`
  - `chr20:4125629 A>AC called 1/1 truth ./. NALT=2 NREF=0 F=1 HQ=57.434779892759934`
  - `chr20:4619359 TC>T called 1/1 truth ./. NALT=2 NREF=57 F=0.03389830508474576 HQ=525.0988928856241`
- *allele/representation mismatch near a truth variant | TR51-200*
  - `chr20:3310333 A>AAG called 0/1 truth ./. NALT=3 NREF=10 F=0.23076923076923078 HQ=108.09735334890388`
  - `chr20:3310333 A>AAAAGAAAG called 1/1 truth ./. NALT=12 NREF=0 F=1 HQ=525.7460709845686`
  - `chr20:3310333 A>AAAAGAAAGAAAGAAAG called 1/1 truth ./. NALT=19 NREF=0 F=1 HQ=267.7263269141321`
- *genotype error (allele right, GT wrong) | TR<=50*
  - `chr20:5454164 AATATAT>A called 1/1 truth 0/1 NALT=21 NREF=1 F=0.9545454545454546 HQ=397.78097319802646`
  - `chr20:9553557 A>AGT called 1/1 truth 0/1 NALT=28 NREF=23 F=0.5490196078431373 HQ=209.91206963947798`
  - `chr20:10476142 C>CTG called 1/1 truth 0/1 NALT=25 NREF=0 F=1 HQ=297.3237529881678`

## Error taxonomy: expAB_67x


**False negatives by cause**

| Cause | Count |
|---|---|
| different allele/representation called nearby | 135 |
| multi-allelic site (truth 1/2), alleles genotyped separately | 126 |
| genotype wrong (allele found, GT differs) | 84 |
| never a candidate (discovery miss) | 51 |
| candidate found but filtered (LowProb/UNCERTAIN) | 3 |

**False positives by cause and repeat context**

| Cause \| context | Count |
|---|---|
| allele/representation mismatch near a truth variant \| HP>=7 | 267 |
| allele/representation mismatch near a truth variant \| TR<=50 | 96 |
| genotype error (allele right, GT wrong) \| HP>=7 | 68 |
| allele/representation mismatch near a truth variant \| other | 41 |
| allele/representation mismatch near a truth variant \| TR51-200 | 39 |
| no truth variant nearby (pure FP) \| other | 34 |
| genotype error (allele right, GT wrong) \| other | 8 |
| genotype error (allele right, GT wrong) \| TR<=50 | 8 |
| genotype error (allele right, GT wrong) \| TR51-200 | 5 |
| no truth variant nearby (pure FP) \| HP>=7 | 4 |
| no truth variant nearby (pure FP) \| TR51-200 | 1 |
| no truth variant nearby (pure FP) \| TR<=50 | 1 |

**FP called genotype**

- allele/representation mismatch near a truth variant, called 0/1: 338
- allele/representation mismatch near a truth variant, called 1/1: 77
- genotype error, called 0/1: 67
- no truth variant nearby, called 0/1: 40
- allele/representation mismatch near a truth variant, called 2/1: 28
- genotype error, called 1/1: 22

**Recall: clustered (another truth variant within 10 bp) vs isolated**

- isolated: 10105/10352 = 0.9761
- clustered: 496/648 = 0.7654

**Median features: TP calls vs pure-FP calls**

| Feature | TP | pure FP |
|---|---|---|
| NALT | n/a | n/a |
| NREF | n/a | n/a |
| F | n/a | n/a |
| HQ | n/a | n/a |
| QMIN | n/a | n/a |

**Examples**

- *genotype wrong (allele found, GT differs)*
  - `chr20:217716 ATG>A truth GT 0/1 [TR51-200] -> ATG>A PASS K=1;DP=42`
  - `chr20:649701 CA>C truth GT 0/1 [HP>=7] -> CA>C PASS K=2;DP=47`
  - `chr20:1122295 CTTTT>C truth GT 1/1 [HP>=7] -> CTTTT>C PASS K=2;DP=45`
  - `chr20:2542038 G>GA truth GT 1/1 [HP>=7] -> G>GA PASS K=3;DP=63`
- *different allele/representation called nearby*
  - `chr20:367246 AC>A truth GT 0/1 [TR<=50] -> AAAC>A PASS K=2;DP=52`
  - `chr20:367250 AAC>A truth GT 0/1 [TR<=50] -> AAAC>A PASS K=2;DP=52`
  - `chr20:499958 CA>C truth GT 0/1 [HP>=7] -> CAA>C PASS K=4;DP=73`
  - `chr20:942273 A>AT truth GT 0/1 [TR<=50] -> A>AT PASS K=3;DP=58`
- *multi-allelic site (truth 1/2), alleles genotyped separately*
  - `chr20:511703 TACAC>CACAC,T truth GT 2/1 [TR<=50] -> TACAC>T PASS K=1;DP=68`
  - `chr20:942277 A>AT,T truth GT 1/2 [TR<=50] -> A>AT PASS K=3;DP=58`
  - `chr20:1903011 CAAAAAAAAA>CAAAAAAAAAAAAA,C truth GT 2/1 [HP>=7] -> CAAAAAAAAA>C PASS K=1;DP=40`
  - `chr20:2647120 TTTTATTTA>TTTTATTTATTTATTTA,T truth GT 2/1 [TR<=50] -> TTTTATTTA>T,TTTTATTTATTTATTTA PASS K=2;DP=40`
- *never a candidate (discovery miss)*
  - `chr20:887636 CT>C truth GT 1/1 [HP>=7] -> no candidate`
  - `chr20:1404738 CTTT>C truth GT 0/1 [HP>=7] -> no candidate`
  - `chr20:1914466 AG>A truth GT 1/1 [other] -> no candidate`
  - `chr20:1914707 CT>C truth GT 1/1 [other] -> no candidate`
- *candidate found but filtered (LowProb/UNCERTAIN)*
  - `chr20:32208648 CAAA>CAAAA,C truth GT 2/1 [HP>=7] -> C>CA LowProb K=1;DP=7`
  - `chr20:53745828 C>CA truth GT 0/1 [HP>=7] -> AC>A LowProb K=1;DP=64`
  - `chr20:59524369 G>GA truth GT 0/1 [HP>=7] -> A>AAG LowProb K=2;DP=58`
- *genotype error (allele right, GT wrong) | TR51-200*
  - `chr20:217716 ATG>A called 1/1 truth 0/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:16801236 T>TTG called 0/1 truth 1/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:22697102 C>CCTTCT called 0/1 truth 1/1 NALT=None NREF=None F=None HQ=None`
- *allele/representation mismatch near a truth variant | TR<=50*
  - `chr20:367249 AAAC>A called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:511703 TACAC>T called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:942277 A>AT called 1/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *allele/representation mismatch near a truth variant | HP>=7*
  - `chr20:458330 TA>T called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:499958 CAA>C called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:560062 AC>A called 1/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *genotype error (allele right, GT wrong) | HP>=7*
  - `chr20:649701 CA>C called 1/1 truth 0/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:1122295 CTTTT>C called 0/1 truth 1/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:2542038 G>GA called 0/1 truth 1/1 NALT=None NREF=None F=None HQ=None`
- *allele/representation mismatch near a truth variant | TR51-200*
  - `chr20:3982103 GTC>G called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:3982105 C>CTCTCTCTCTC called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:6124250 G>GCA,GCACA called 2/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *no truth variant nearby (pure FP) | HP>=7*
  - `chr20:4208933 TG>T called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:25926126 CA>C called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:25928260 AT>A called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *no truth variant nearby (pure FP) | other*
  - `chr20:4374492 CA>C called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:4374503 CTACCTCCAACCTT>C called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:7508395 G>GA called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *allele/representation mismatch near a truth variant | other*
  - `chr20:4513315 AG>A called 1/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:4619359 TC>T called 1/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:9051911 T>TA called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *genotype error (allele right, GT wrong) | other*
  - `chr20:9051906 G>GA called 0/1 truth 1/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:12759430 CG>C called 0/1 truth 1/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:19711132 A>AAC called 0/1 truth 1/1 NALT=None NREF=None F=None HQ=None`
- *genotype error (allele right, GT wrong) | TR<=50*
  - `chr20:9553557 A>AGT called 1/1 truth 0/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:40040821 G>GTCTC called 1/1 truth 0/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:43998719 A>ACATC called 1/1 truth 0/1 NALT=None NREF=None F=None HQ=None`
- *no truth variant nearby (pure FP) | TR51-200*
  - `chr20:25895715 AT>A called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *no truth variant nearby (pure FP) | TR<=50*
  - `chr20:26079727 ACTAATCACTGGCTGCAACT>A called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`

## Error taxonomy: expABC2_67x


**False negatives by cause**

| Cause | Count |
|---|---|
| genotype wrong (allele found, GT differs) | 194 |
| different allele/representation called nearby | 120 |
| multi-allelic site (truth 1/2), alleles genotyped separately | 114 |
| candidate found but filtered (LowProb/UNCERTAIN) | 89 |
| never a candidate (discovery miss) | 51 |

**False positives by cause and repeat context**

| Cause \| context | Count |
|---|---|
| allele/representation mismatch near a truth variant \| HP>=7 | 229 |
| allele/representation mismatch near a truth variant \| TR<=50 | 90 |
| genotype error (allele right, GT wrong) \| HP>=7 | 52 |
| allele/representation mismatch near a truth variant \| TR51-200 | 40 |
| no truth variant nearby (pure FP) \| other | 32 |
| allele/representation mismatch near a truth variant \| other | 30 |
| genotype error (allele right, GT wrong) \| TR<=50 | 10 |
| genotype error (allele right, GT wrong) \| other | 7 |
| genotype error (allele right, GT wrong) \| TR51-200 | 3 |
| no truth variant nearby (pure FP) \| HP>=7 | 3 |
| no truth variant nearby (pure FP) \| TR51-200 | 1 |
| no truth variant nearby (pure FP) \| TR<=50 | 1 |

**FP called genotype**

- allele/representation mismatch near a truth variant, called 0/1: 289
- allele/representation mismatch near a truth variant, called 1/1: 72
- genotype error, called 0/1: 50
- no truth variant nearby, called 0/1: 37
- allele/representation mismatch near a truth variant, called 2/1: 28
- genotype error, called 1/1: 22

**Recall: clustered (another truth variant within 10 bp) vs isolated**

- isolated: 9940/10352 = 0.9602
- clustered: 492/648 = 0.7593

**Median features: TP calls vs pure-FP calls**

| Feature | TP | pure FP |
|---|---|---|
| NALT | n/a | n/a |
| NREF | n/a | n/a |
| F | n/a | n/a |
| HQ | n/a | n/a |
| QMIN | n/a | n/a |

**Examples**

- *genotype wrong (allele found, GT differs)*
  - `chr20:217716 ATG>A truth GT 0/1 [TR51-200] -> ATG>A PASS K=1;DP=36`
  - `chr20:240075 C>CAAA truth GT 1/1 [HP>=7] -> C>CAAA LowProb K=2;DP=24`
  - `chr20:577957 CA>C truth GT 1/1 [HP>=7] -> TCA>TC,TA PASS K=2;DP=23`
  - `chr20:649701 CA>C truth GT 0/1 [HP>=7] -> CA>C PASS K=2;DP=44`
- *candidate found but filtered (LowProb/UNCERTAIN)*
  - `chr20:290494 G>GC truth GT 1/1 [other] -> G>GC LowProb K=1;DP=23`
  - `chr20:787512 CT>C truth GT 1/1 [HP>=7] -> CT>C LowProb K=2;DP=30`
  - `chr20:828060 G>GAA,GAAAA truth GT 1/2 [HP>=7] -> G>GAA LowProb K=2;DP=22`
  - `chr20:1311900 C>CTT,CTTTTTT truth GT 1/2 [HP>=7] -> C>CTT LowProb K=3;DP=22`
- *different allele/representation called nearby*
  - `chr20:367246 AC>A truth GT 0/1 [TR<=50] -> AAAC>A PASS K=2;DP=45`
  - `chr20:367250 AAC>A truth GT 0/1 [TR<=50] -> AAAC>A PASS K=2;DP=45`
  - `chr20:499958 CA>C truth GT 0/1 [HP>=7] -> CAA>C PASS K=2;DP=58`
  - `chr20:942273 A>AT truth GT 0/1 [TR<=50] -> A>AT,ATTTTT PASS K=3;DP=51`
- *multi-allelic site (truth 1/2), alleles genotyped separately*
  - `chr20:511703 TACAC>CACAC,T truth GT 2/1 [TR<=50] -> TACAC>T PASS K=1;DP=61`
  - `chr20:942277 A>AT,T truth GT 1/2 [TR<=50] -> A>AT,ATTTTT PASS K=3;DP=51`
  - `chr20:1903011 CAAAAAAAAA>CAAAAAAAAAAAAA,C truth GT 2/1 [HP>=7] -> CAAAAAAAAA>C PASS K=1;DP=36`
  - `chr20:2314530 TACACAC>TACACACACAC,T truth GT 2/1 [TR<=50] -> TACACAC>T PASS K=2;DP=34`
- *never a candidate (discovery miss)*
  - `chr20:1404738 CTTT>C truth GT 0/1 [HP>=7] -> no candidate`
  - `chr20:1914466 AG>A truth GT 1/1 [other] -> no candidate`
  - `chr20:1914707 CT>C truth GT 1/1 [other] -> no candidate`
  - `chr20:1915067 T>TAA truth GT 1/1 [other] -> no candidate`
- *genotype error (allele right, GT wrong) | TR51-200*
  - `chr20:217716 ATG>A called 1/1 truth 0/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:35321252 C>CTG called 1/1 truth 0/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:61952103 A>AAC called 1/1 truth 0/1 NALT=None NREF=None F=None HQ=None`
- *allele/representation mismatch near a truth variant | TR<=50*
  - `chr20:367249 AAAC>A called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:511703 TACAC>T called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:942277 A>AT,ATTTTT called 2/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *allele/representation mismatch near a truth variant | HP>=7*
  - `chr20:499958 CAA>C called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:560057 C>CA called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:577956 TC>T called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *genotype error (allele right, GT wrong) | HP>=7*
  - `chr20:577957 CA>C called 0/1 truth 1/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:649701 CA>C called 1/1 truth 0/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:2542038 G>GA called 0/1 truth 1/1 NALT=None NREF=None F=None HQ=None`
- *allele/representation mismatch near a truth variant | TR51-200*
  - `chr20:672388 CGT>C called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:3982103 GTC>G called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:3982105 C>CTCTCTCTCTC called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *no truth variant nearby (pure FP) | other*
  - `chr20:4374492 CA>C called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:4374503 CTACCTCCAACCTT>C called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:7508395 G>GA called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *allele/representation mismatch near a truth variant | other*
  - `chr20:4619359 TC>T called 1/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:9051911 T>TA called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:11369715 C>CA called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *genotype error (allele right, GT wrong) | other*
  - `chr20:9051906 G>GA called 0/1 truth 1/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:12759430 CG>C called 0/1 truth 1/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:19711132 A>AAC called 0/1 truth 1/1 NALT=None NREF=None F=None HQ=None`
- *genotype error (allele right, GT wrong) | TR<=50*
  - `chr20:9553557 A>AGT called 1/1 truth 0/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:40040821 G>GTCTC called 1/1 truth 0/1 NALT=None NREF=None F=None HQ=None`
  - `chr20:43998719 A>ACATC called 1/1 truth 0/1 NALT=None NREF=None F=None HQ=None`
- *no truth variant nearby (pure FP) | TR51-200*
  - `chr20:25895715 AT>A called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *no truth variant nearby (pure FP) | HP>=7*
  - `chr20:25926126 CA>C called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:25928260 AT>A called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
  - `chr20:26077020 GA>G called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
- *no truth variant nearby (pure FP) | TR<=50*
  - `chr20:26079727 ACTAATCACTGGCTGCAACT>A called 0/1 truth ./. NALT=None NREF=None F=None HQ=None`
