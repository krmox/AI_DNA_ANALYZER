"""Generate a synthetic BAM with VCF variants spliced into reads.

Splices each VCF record's ref/alt pair into ~50% of overlapping reads
(heterozygous, VAF ~ 0.5 — the same per-read coin flip already used for
SNPs), producing correct CIGAR operations for SNPs, insertions and
deletions alike. Previously only the SNP branch actually edited a read;
Insertion/Deletion VCF records were read but never spliced in, so every
read's CIGAR was a flat all-match block and the pileup carried zero
evidence for those labels no matter what the VCF said.

Extension point: to support a new variant shape (MNP, symbolic ``<DEL>``,
...), add one classification branch to :func:`classify_variant` and one
splicer function matching the ``VariantSplicer`` signature, then register
it in :data:`VARIANT_SPLICERS`. Nothing else in this module needs to change.
"""

from __future__ import annotations

import argparse
import os
import random
from typing import Callable, Dict, List, Tuple

import pysam

FASTA_PATH = "data/giab_chr21/chr21_slice.fa"
VCF_PATH = "data/giab_chr21/chr21_slice.vcf.gz"
BAM_PATH = "data/giab_chr21/chr21_slice.bam"

READ_LEN = 150
COVERAGE = 20
BASE_QUALITY = 30
MAPPING_QUALITY = 60
#: Per-read probability of carrying the ALT allele at every variant it
#: overlaps. A single coin flip per read (not per variant) is what makes
#: this "heterozygous, VAF ~ 0.5" rather than independently noisy per site.
HET_ALT_READ_FRACTION = 0.5

CIGAR_MATCH = 0
CIGAR_INSERTION = 1
CIGAR_DELETION = 2

Variant = Tuple[str, str]  # (ref, alt)
CigarOp = Tuple[int, int]

#: A splicer turns one variant anchored at absolute reference position
#: ``pos`` into the query bases it contributes, the CIGAR ops covering
#: them, and how many reference positions it consumes — e.g. a deletion
#: consumes more reference than it contributes query bases, so the
#: caller's window walk knows how far to skip forward.
VariantSplicer = Callable[[str, str, List[str], int], Tuple[List[str], List[CigarOp], int]]


def load_variants(vcf_path: str) -> Dict[int, Variant]:
    """Read (ref, alt) pairs from a VCF, keyed by 0-based position.

    Args:
        vcf_path: Path to the (optionally bgzipped/indexed) VCF.

    Returns:
        Mapping of 0-based reference position to ``(ref, alt)``.
    """
    variants: Dict[int, Variant] = {}
    with pysam.VariantFile(vcf_path) as vcf:
        for record in vcf.fetch():
            variants[record.pos - 1] = (record.ref, record.alts[0])
    return variants


def classify_variant(ref: str, alt: str) -> str | None:
    """Classify a VCF ref/alt pair the same way the read splicer dispatches on it.

    Only true single-base SNPs, pure insertions and pure deletions are
    recognised — the only shapes this project's mock VCF generator
    (``generate_mock_vcf.py``) ever emits. Anything else (e.g. an MNP)
    returns ``None`` and is left un-spliced, exactly as the original code
    silently skipped everything but ``len(ref) == len(alt) == 1``.

    Args:
        ref: VCF REF allele.
        alt: VCF ALT allele.

    Returns:
        ``"snp"``, ``"insertion"``, ``"deletion"``, or ``None``.
    """
    if len(ref) == 1 and len(alt) == 1:
        return "snp"
    if len(alt) > len(ref):
        return "insertion"
    if len(ref) > len(alt):
        return "deletion"
    return None


def _splice_snp(ref: str, alt: str, seq: List[str], pos: int) -> Tuple[List[str], List[CigarOp], int]:
    """Substitute the reference base with ALT. Unchanged from the original behaviour."""
    del seq, pos
    return [alt], [(CIGAR_MATCH, 1)], len(ref)


def _splice_insertion(
    ref: str, alt: str, seq: List[str], pos: int
) -> Tuple[List[str], List[CigarOp], int]:
    """Keep the anchor base, insert ALT's extra bases as a real 'I' op.

    VCF spells an insertion as e.g. ``A -> AT``: the leading base is an
    anchor present in both alleles, and ``alt[len(ref):]`` is what is
    actually new and has no reference coordinate of its own.
    """
    anchor = seq[pos]
    inserted = list(alt[len(ref):])
    return [anchor, *inserted], [(CIGAR_MATCH, 1), (CIGAR_INSERTION, len(inserted))], len(ref)


def _splice_deletion(
    ref: str, alt: str, seq: List[str], pos: int
) -> Tuple[List[str], List[CigarOp], int]:
    """Keep the anchor base, drop the deleted reference bases with a real 'D' op."""
    anchor = seq[pos]
    deleted_len = len(ref) - len(alt)
    return [anchor], [(CIGAR_MATCH, 1), (CIGAR_DELETION, deleted_len)], len(ref)


VARIANT_SPLICERS: Dict[str, VariantSplicer] = {
    "snp": _splice_snp,
    "insertion": _splice_insertion,
    "deletion": _splice_deletion,
}


def _merge_cigar_ops(ops: List[CigarOp]) -> List[CigarOp]:
    """Collapse consecutive same-op CIGAR entries, as htslib expects a minimal CIGAR."""
    merged: List[CigarOp] = []
    for op, length in ops:
        if merged and merged[-1][0] == op:
            merged[-1] = (op, merged[-1][1] + length)
        else:
            merged.append((op, length))
    return merged


def build_read_window(
    start: int,
    read_len: int,
    seq: List[str],
    variants: Dict[int, Variant],
    carries_alt: bool,
) -> Tuple[str, List[CigarOp]]:
    """Build one read's query sequence and CIGAR over a fixed reference window.

    The reference footprint is always exactly ``read_len`` positions
    (``start`` .. ``start + read_len``), which is what keeps coverage/depth
    unaffected by indel injection: an insertion adds query bases beyond
    that footprint ('I' consumes no reference), and a deletion consumes
    reference positions inside the footprint without adding query bases
    ('D' consumes no query). Reads that don't carry the alt allele are
    copied verbatim with a single match block, identical to the original
    generator's output.

    A variant that would extend past the window's end is left un-spliced
    rather than truncated, so it can never produce a malformed CIGAR for a
    read near a region boundary.

    Args:
        start: 0-based reference start of this read's window.
        read_len: Fixed reference footprint width.
        seq: Full contig reference sequence, one base per list element.
        variants: All variants keyed by 0-based reference position.
        carries_alt: Whether this read is the "alt" haplotype for every
            variant it overlaps.

    Returns:
        Tuple of ``(query_sequence, cigar_ops)``. The reference-consuming
        ops ('M' + 'D') always sum to ``min(read_len, len(seq) - start)``.
    """
    end = min(start + read_len, len(seq))

    if not carries_alt:
        bases = seq[start:end]
        cigar_ops = [(CIGAR_MATCH, end - start)] if end > start else []
        return "".join(bases), cigar_ops

    query_chars: List[str] = []
    cigar_ops: List[CigarOp] = []
    pos = start
    while pos < end:
        variant = variants.get(pos)
        variant_type = classify_variant(*variant) if variant is not None else None
        splice = VARIANT_SPLICERS.get(variant_type) if variant_type is not None else None

        if splice is not None:
            ref, alt = variant
            if pos + len(ref) <= end:
                chars, ops, consumed = splice(ref, alt, seq, pos)
                query_chars.extend(chars)
                cigar_ops.extend(ops)
                pos += consumed
                continue

        query_chars.append(seq[pos])
        cigar_ops.append((CIGAR_MATCH, 1))
        pos += 1

    return "".join(query_chars), _merge_cigar_ops(cigar_ops)


def _sanitize_bases(sequence: str) -> str:
    """Replace anything outside ACGT (case-insensitive) with 'A', as the original generator did."""
    return "".join(char if char in "ACGTacgt" else "A" for char in sequence)


def generate_reads(
    seq: List[str],
    variants: Dict[int, Variant],
    read_len: int = READ_LEN,
    coverage: int = COVERAGE,
    het_alt_read_fraction: float = HET_ALT_READ_FRACTION,
    rng: random.Random | None = None,
) -> List[Tuple[int, str, List[CigarOp]]]:
    """Simulate reads tiling ``seq`` with variants spliced into ~half of them.

    Args:
        seq: Reference sequence, one base per element.
        variants: Variants keyed by 0-based position.
        read_len: Read length.
        coverage: Target average depth.
        het_alt_read_fraction: Per-read probability of carrying the alt
            allele — see :data:`HET_ALT_READ_FRACTION`.
        rng: Random source; defaults to the module-level ``random`` state.
            Pass a seeded ``random.Random`` for deterministic tests.

    Returns:
        List of ``(reference_start, query_sequence, cigar_ops)`` tuples.
    """
    rng = rng if rng is not None else random
    seq_len = len(seq)
    num_reads = int((seq_len * coverage) / read_len)

    reads: List[Tuple[int, str, List[CigarOp]]] = []
    for _ in range(num_reads):
        start = rng.randint(0, max(0, seq_len - read_len))
        carries_alt = rng.random() < het_alt_read_fraction
        query_sequence, cigar_ops = build_read_window(start, read_len, seq, variants, carries_alt)
        reads.append((start, _sanitize_bases(query_sequence), cigar_ops))
    return reads


def parse_args() -> argparse.Namespace:
    """Path overrides, all defaulting to this module's constants.

    Paths only, deliberately: ``READ_LEN``, ``COVERAGE``, ``BASE_QUALITY``,
    ``MAPPING_QUALITY`` and ``HET_ALT_READ_FRACTION`` are not exposed, so a
    second region is generated under exactly the process that produced the
    first. ``main()`` already derives contig name and length from the FASTA
    itself, so a different-width region needs no other change.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", default=FASTA_PATH, help="Indexed reference FASTA.")
    parser.add_argument("--vcf", default=VCF_PATH, help="Variants to splice into reads.")
    parser.add_argument("--bam", default=BAM_PATH, help="Output BAM path.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    print(f"🧬 Генерация BAM с внедрением VCF-мутаций (VAF ~ {HET_ALT_READ_FRACTION:.2f}) для...")

    variants = load_variants(args.vcf)

    fasta = pysam.FastaFile(args.fasta)
    contig = fasta.references[0]
    seq = list(fasta.fetch(contig))
    fasta.close()

    header = {"HD": {"VN": "1.0"}, "SQ": [{"LN": len(seq), "SN": contig}]}

    reads = generate_reads(seq, variants)

    with pysam.AlignmentFile(args.bam, "wb", header=header) as out_bam:
        for index, (start, query_sequence, cigar_ops) in enumerate(reads):
            alignment = pysam.AlignedSegment()
            alignment.query_name = f"read_{index}"
            alignment.query_sequence = query_sequence
            alignment.flag = 0
            alignment.reference_id = 0
            alignment.reference_start = start
            alignment.mapping_quality = MAPPING_QUALITY
            alignment.cigar = cigar_ops
            alignment.query_qualities = [BASE_QUALITY] * len(query_sequence)
            out_bam.write(alignment)

    print("⚙️ Сортировка и индексирование BAM...")
    # Derived from the output path rather than hardcoded, so generating into
    # a different directory never drops a temp file into (or moves one out
    # of) the original dataset's directory.
    sorted_bam = os.path.splitext(args.bam)[0] + "_sorted.bam"
    pysam.sort("-o", sorted_bam, args.bam)
    pysam.index(sorted_bam)

    os.replace(sorted_bam, args.bam)
    os.replace(sorted_bam + ".bai", args.bam + ".bai")

    print("✅ Валидный BAM файл готов: SNP, Insertion и Deletion теперь спроецированы в риды!")


if __name__ == "__main__":
    main()
