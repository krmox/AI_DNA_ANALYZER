#!/usr/bin/env bash
# Regression: invalid --region values must fail with a non-zero exit status (they used to exit 0 with an empty VCF).
# usage: test_cli_region.sh NINDELS_CALL BAM FASTA CONTIG
set -u
N=$1; BAM=$2; REF=$3; C=$4; fail=0
expect_fail() { if "$N" call --bam "$BAM" --ref "$REF" --region "$1" --out-vcf /dev/null >/dev/null 2>&1; then echo "FAIL: region '$1' accepted"; fail=1; else echo "ok:   region '$1' rejected"; fi; }
expect_ok()   { if "$N" call --bam "$BAM" --ref "$REF" --region "$1" --out-vcf /dev/null >/dev/null 2>&1; then echo "ok:   region '$1' accepted"; else echo "FAIL: region '$1' rejected"; fail=1; fi; }
expect_fail "$C:200000-100000"
expect_fail "$C:0-1000"
expect_fail "$C:900000000-1000000000"
expect_ok   "$C:100000-120000"
exit $fail
