"""Generate a synthetic VCF whose REF alleles come from the reference FASTA.

REF used to be drawn at random from ACGT instead of read from the FASTA,
so a SNP's ALT landed on the true reference base ~25% of the time. Those
records are *phantom variants*: ``generate_mock_bam.py`` splices ALT into
the reads correctly, but when ALT equals the real reference base the read
is byte-identical to a non-variant read, leaving a label with no evidence
behind it. Reading REF from the FASTA makes ALT != reference base by
construction, so every emitted variant is actually observable in the BAM.

Coordinates follow the same convention as the BAM generator: positions are
0-based half-open internally and written to the VCF as 1-based POS
(``pos_0based + 1``), the inverse of ``load_variants``' ``record.pos - 1``.

Positions whose reference window is not pure ACGT (the slice contains long
N runs) are skipped rather than emitted, because a non-ACGT REF cannot be
represented faithfully in the reads — the BAM generator sanitises anything
outside ACGT to 'A', which would recreate the same evidence-free label.
"""

import random

import pysam

FASTA_PATH = "data/giab_chr21/chr21_slice.fa"
vcf_path = "data/giab_chr21/chr21_slice.vcf"
vcf_gz_path = "data/giab_chr21/chr21_slice.vcf.gz"

#: Variant spacing in bp — ~1 variant per 1000 bp, as in a real genome.
VARIANT_STRIDE = 1000
REGION_START = 1000
REGION_END = 5000000

SNP_PROB = 0.8
INSERTION_PROB = 0.1  # remaining 0.1 is deletions

bases = ['A', 'C', 'G', 'T']


def make_record(seq: str, pos0: int, rand: float, rng: random.Random) -> tuple[str, str] | None:
    """Build one (REF, ALT) pair anchored at 0-based position ``pos0``.

    Args:
        seq: Full contig sequence.
        pos0: 0-based reference position of the variant anchor.
        rand: Uniform draw selecting the variant class.
        rng: Random source for the ALT allele.

    Returns:
        ``(ref, alt)``, or ``None`` if the reference window at this
        position is not pure ACGT and the site must be skipped.
    """
    if rand < SNP_PROB:
        ref = seq[pos0].upper()
        if ref not in bases:
            return None
        return ref, rng.choice([b for b in bases if b != ref])

    if rand < SNP_PROB + INSERTION_PROB:  # insertion: anchor base + one new base
        ref = seq[pos0].upper()
        if ref not in bases:
            return None
        return ref, ref + 'T'

    # deletion: anchor base kept, following reference base removed
    ref = seq[pos0:pos0 + 2].upper()
    if len(ref) < 2 or any(b not in bases for b in ref):
        return None
    return ref, ref[0]


def main() -> None:
    fasta = pysam.FastaFile(FASTA_PATH)
    contig = fasta.references[0]
    contig_len = fasta.lengths[0]
    seq = fasta.fetch(contig)
    fasta.close()

    header = f"""##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID={contig},length={contig_len}>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002
"""

    lines = [header]
    rng = random
    # ``pos`` is the 1-based VCF POS, unchanged from the original generator so
    # the variant loci land on exactly the same coordinates as before; the
    # FASTA is indexed 0-based at ``pos - 1``, matching ``load_variants``.
    for pos in range(REGION_START, REGION_END, VARIANT_STRIDE):
        record = make_record(seq, pos - 1, rng.random(), rng)
        if record is None:
            continue
        ref, alt = record
        lines.append(f"{contig}\t{pos}\t.\t{ref}\t{alt}\t99\tPASS\t.\tGT\t0/1\n")

    with open(vcf_path, "w") as f:
        f.writelines(lines)

    pysam.tabix_compress(vcf_path, vcf_gz_path, force=True)
    pysam.tabix_index(vcf_gz_path, preset="vcf", force=True)
    print(f"✅ Синтетический VCF создан и заиндексирован! Записей: {len(lines) - 1}")


if __name__ == "__main__":
    main()
