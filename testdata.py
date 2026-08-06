"""Builds genuine indexed FASTA/BAM/VCF fixtures for provider testing.

These are *real* files written and read back through pysam — real BAM
records with real CIGAR strings, a real ``.fai``, a real tabix index, and a
deliberately mismatched ``chr`` prefix between files. That exercises every
code path in :class:`~.providers.GiabAlignmentProvider`: contig
normalisation, the 0-based/1-based VCF conversion, the CIGAR walk, pileup
consensus and VAF.

What this does **not** substitute for. A fixture cannot reproduce real
chr21: uniform base composition instead of isochores and repeats, planted
variants at a density chosen for testing rather than the true ~1:1000,
uniform coverage rather than the depth variation that drives real false
positives, and no mapping ambiguity at all because reads are placed at
known positions rather than aligned. Passing against this fixture means the
parser is *correct*, not that the model will work on GIAB.
"""

from __future__ import annotations

import random
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Tuple

import pysam

NUCLEOTIDES: Tuple[str, ...] = ("A", "C", "G", "T")


@dataclass
class PlantedVariant:
    """A variant deliberately inserted into the fixture.

    Attributes:
        position: 0-based reference position.
        kind: One of ``"snp"``, ``"insertion"`` or ``"deletion"``.
        ref_allele: Reference allele string in VCF convention.
        alt_allele: Alternate allele string in VCF convention.
        genotype: VCF genotype tuple.
    """

    position: int
    kind: str
    ref_allele: str
    alt_allele: str
    genotype: Tuple[int, int] = (0, 1)


@dataclass
class Fixture:
    """Paths and truth data for a generated fixture.

    Attributes:
        fasta_path: Path to the indexed FASTA.
        bam_path: Path to the indexed BAM.
        vcf_path: Path to the bgzipped, tabix-indexed VCF.
        bed_path: Path to the high-confidence BED.
        reference: The reference sequence written to the FASTA.
        variants: The planted variants, in position order.
        fasta_contig: Contig spelling used in the FASTA.
        bam_contig: Contig spelling used in the BAM.
        vcf_contig: Contig spelling used in the VCF.
    """

    fasta_path: Path
    bam_path: Path
    vcf_path: Path
    bed_path: Path
    reference: str
    variants: List[PlantedVariant] = field(default_factory=list)
    fasta_contig: str = "chr21"
    bam_contig: str = "chr21"
    vcf_contig: str = "21"


def _apply_variants_to_read(
    reference: str,
    start: int,
    length: int,
    variants: Dict[int, PlantedVariant],
    carry: Dict[int, bool],
    rng: random.Random,
    error_rate: float,
) -> Tuple[str, List[Tuple[int, int]]]:
    """Build one read sequence and its CIGAR against the reference.

    Args:
        reference: Full reference sequence.
        start: 0-based read start on the reference.
        length: Number of reference bases the read spans.
        variants: Planted variants keyed by position.
        carry: Whether this read carries each variant, keyed by position.
            Heterozygous sites are carried by roughly half the reads, which
            is what produces a realistic VAF near 0.5.
        rng: Random source.
        error_rate: Per-base sequencing error rate.

    Returns:
        Tuple of ``(sequence, cigar)`` where cigar is a list of
        ``(operation, length)`` pairs in pysam's numeric encoding.
    """
    sequence: List[str] = []
    cigar: List[Tuple[int, int]] = []

    def push(operation: int, count: int) -> None:
        """Append to the CIGAR, merging with the previous run if possible."""
        if cigar and cigar[-1][0] == operation:
            cigar[-1] = (operation, cigar[-1][1] + count)
        else:
            cigar.append((operation, count))

    position = start
    end = start + length

    while position < end:
        variant = variants.get(position)
        carries = variant is not None and carry.get(position, False)

        if not carries:
            base = reference[position]
            if rng.random() < error_rate:
                base = rng.choice([n for n in NUCLEOTIDES if n != base])
            sequence.append(base)
            push(0, 1)  # M
            position += 1
            continue

        assert variant is not None

        if variant.kind == "snp":
            sequence.append(variant.alt_allele)
            push(0, 1)
            position += 1

        elif variant.kind == "insertion":
            # VCF anchors the insertion on the preceding base, so the anchor
            # is an M and the novel bases are an I.
            sequence.append(variant.ref_allele)
            push(0, 1)
            inserted = variant.alt_allele[1:]
            sequence.append(inserted)
            push(1, len(inserted))  # I
            position += 1

        else:  # deletion
            sequence.append(variant.ref_allele[0])
            push(0, 1)
            deleted = len(variant.ref_allele) - len(variant.alt_allele)
            push(2, deleted)  # D
            position += 1 + deleted

    return "".join(sequence), cigar


def build_fixture(
    directory: Path,
    contig_length: int = 4000,
    depth: int = 20,
    read_length: int = 120,
    variant_spacing: int = 180,
    error_rate: float = 0.01,
    seed: int = 7,
    fasta_contig: str = "chr21",
    bam_contig: str = "chr21",
    vcf_contig: str = "21",
) -> Fixture:
    """Generate a complete indexed FASTA/BAM/VCF/BED fixture.

    The VCF contig spelling defaults to ``21`` while the FASTA and BAM use
    ``chr21``, so the fixture exercises contig normalisation by default
    rather than only on request.

    Args:
        directory: Output directory; created if absent.
        contig_length: Reference length in bases.
        depth: Approximate read depth.
        read_length: Reference span of each read.
        variant_spacing: Distance between planted variants.
        error_rate: Per-base sequencing error rate in the reads.
        seed: Random seed.
        fasta_contig: Contig spelling for the FASTA.
        bam_contig: Contig spelling for the BAM.
        vcf_contig: Contig spelling for the VCF.

    Returns:
        A :class:`Fixture` describing the generated files.
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    rng = random.Random(seed)

    reference = "".join(rng.choice(NUCLEOTIDES) for _ in range(contig_length))

    # --- plant variants -------------------------------------------------
    variants: List[PlantedVariant] = []
    kinds = ("snp", "insertion", "deletion")
    for index, position in enumerate(range(200, contig_length - 200, variant_spacing)):
        kind = kinds[index % 3]
        reference_base = reference[position]

        if kind == "snp":
            alternate = rng.choice([n for n in NUCLEOTIDES if n != reference_base])
            variants.append(PlantedVariant(position, kind, reference_base, alternate))

        elif kind == "insertion":
            inserted = "".join(rng.choice(NUCLEOTIDES) for _ in range(rng.randint(1, 4)))
            variants.append(
                PlantedVariant(position, kind, reference_base, reference_base + inserted)
            )

        else:  # deletion
            span = rng.randint(1, 3)
            deleted = reference[position : position + 1 + span]
            variants.append(PlantedVariant(position, kind, deleted, reference_base))

    variant_map = {variant.position: variant for variant in variants}

    # --- FASTA ----------------------------------------------------------
    fasta_path = directory / "reference.fa"
    with open(fasta_path, "w") as handle:
        handle.write(f">{fasta_contig}\n")
        for offset in range(0, len(reference), 60):
            handle.write(reference[offset : offset + 60] + "\n")
    pysam.faidx(str(fasta_path))

    # --- BAM ------------------------------------------------------------
    bam_path = directory / "reads.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": bam_contig, "LN": contig_length}],
    }

    records: List[pysam.AlignedSegment] = []
    read_index = 0
    step = max(1, read_length // max(depth // 4, 1))

    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
        for start in range(0, contig_length - read_length, step):
            for replicate in range(max(1, depth // (read_length // step))):
                # Heterozygous sites are carried by a coin flip per read,
                # which is what makes the observed VAF land near 0.5.
                carry = {
                    position: (rng.random() < 0.5 if variant.genotype != (1, 1) else True)
                    for position, variant in variant_map.items()
                }
                sequence, cigar = _apply_variants_to_read(
                    reference, start, read_length, variant_map, carry, rng, error_rate
                )

                segment = pysam.AlignedSegment()
                segment.query_name = f"read_{read_index}_{replicate}"
                segment.query_sequence = sequence
                segment.flag = 0
                segment.reference_id = 0
                segment.reference_start = start
                segment.mapping_quality = 60
                segment.cigartuples = cigar
                segment.query_qualities = pysam.qualitystring_to_array("I" * len(sequence))
                records.append(segment)
                read_index += 1

        for segment in sorted(records, key=lambda record: record.reference_start):
            bam.write(segment)

    pysam.index(str(bam_path))

    # --- VCF ------------------------------------------------------------
    vcf_path = directory / "truth.vcf.gz"
    vcf_header = pysam.VariantHeader()
    vcf_header.add_line(f"##contig=<ID={vcf_contig},length={contig_length}>")
    vcf_header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    vcf_header.add_sample("HG002")

    with pysam.VariantFile(str(vcf_path), "w", header=vcf_header) as vcf:
        for variant in variants:
            record = vcf.new_record(
                contig=vcf_contig,
                # pysam's start is 0-based; it writes the 1-based POS itself.
                start=variant.position,
                stop=variant.position + len(variant.ref_allele),
                alleles=(variant.ref_allele, variant.alt_allele),
            )
            record.samples["HG002"]["GT"] = variant.genotype
            vcf.write(record)

    pysam.tabix_index(str(vcf_path), preset="vcf", force=True)

    # --- BED ------------------------------------------------------------
    bed_path = directory / "highconf.bed"
    with open(bed_path, "w") as handle:
        handle.write(f"{fasta_contig}\t0\t{contig_length}\n")

    return Fixture(
        fasta_path=fasta_path,
        bam_path=bam_path,
        vcf_path=vcf_path,
        bed_path=bed_path,
        reference=reference,
        variants=variants,
        fasta_contig=fasta_contig,
        bam_contig=bam_contig,
        vcf_contig=vcf_contig,
    )


__all__ = ["PlantedVariant", "Fixture", "build_fixture"]
