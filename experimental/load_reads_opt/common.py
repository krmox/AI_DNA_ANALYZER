import sys, os
import pathlib as _pl; _ROOT = str(_pl.Path(__file__).resolve().parents[2])  # project root (path-independent)
W = _ROOT
M = _ROOT
sys.path.insert(0, W); sys.path.insert(0, os.path.join(W, "experimental", "load_reads_opt"))
os.chdir(W)
HG005 = dict(fasta=f"{M}/experimental/stress_test/shared_ref/chr1_chrname.fa", bam=f"{M}/data/giab_hg005_stress/chr1_24Mb_30x.bam",
             vcf=f"{M}/data/giab_hg005_stress/chr1_24Mb.vcf.gz", bed=f"{M}/data/giab_hg005_stress/chr1_24Mb_highconf.bed", contig="chr1")
HG002_15X = dict(fasta=f"{M}/data/reference/chr21_full.fa", bam=f"{M}/data/giab_hg002_chr21_12Mb/hg002_chr21_32_44M_15x.bam",
                 vcf=f"{M}/data/giab_hg002_chr21_12Mb/hg002_chr21_32_44M.vcf.gz", bed=f"{M}/data/giab_hg002_chr21_12Mb/hg002_chr21_32_44M_highconf.bed", contig="chr21")
CHR20 = dict(fasta=f"{M}/data/reference/chr20_full.fa", bam=f"{M}/data/giab_hg002_chr20/chr20_15x.bam",
             vcf=f"{M}/data/giab_hg002_chr20/chr20.vcf.gz", bed=f"{M}/data/giab_hg002_chr20/chr20_highconf.bed", contig="chr20")
DATASETS = {"hg005": HG005, "hg002_15x": HG002_15X, "chr20": CHR20}

def use_original_window_scan():
    """Benchmark/equivalence hook: restore the ORIGINAL O(#intervals) BED scan of _is_confident (the OLD code)."""
    from providers import GiabAlignmentProvider, load_bed_regions
    def _orig(self, start, end):
        if self.high_confidence_bed is None:
            return True
        if self._confident_regions is None:
            self._confident_regions = load_bed_regions(self.high_confidence_bed, self.contig)
        return any(rs <= start and end <= re_ for rs, re_ in self._confident_regions)
    GiabAlignmentProvider._is_confident = _orig
