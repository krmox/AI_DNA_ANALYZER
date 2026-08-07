import os
import random
import pysam

fasta_path = "data/giab_chr21/chr21_slice.fa"
vcf_path = "data/giab_chr21/chr21_slice.vcf.gz"
bam_path = "data/giab_chr21/chr21_slice.bam"

# 1. Считываем мутации из VCF
vcf = pysam.VariantFile(vcf_path)
variants = {}
for rec in vcf.fetch():
    pos_0 = rec.pos - 1  # Приводим к 0-based индексу
    variants[pos_0] = (rec.ref, rec.alts[0])

# 2. Читаем FASTA
fasta = pysam.FastaFile(fasta_path)
contig = fasta.references[0]
seq = list(fasta.fetch(contig))
seq_len = len(seq)

header = {
    'HD': {'VN': '1.0'},
    'SQ': [{'LN': seq_len, 'SN': contig}]
}

read_len = 150
coverage = 20
num_reads = int((seq_len * coverage) / read_len)

print(f"🧬 Генерация BAM с внедрением VCF-мутаций (VAF ~ 0.50) для {contig}...")

with pysam.AlignmentFile(bam_path, "wb", header=header) as out_bam:
    for i in range(num_reads):
        start = random.randint(0, max(0, seq_len - read_len))
        read_bases = list(seq[start:start+read_len])

        # 50% ридов несут альтернативный аллель (гетерозигота, VAF = 0.5)
        is_alt_read = random.random() < 0.5

        if is_alt_read:
            for rel_pos in range(read_len):
                abs_pos = start + rel_pos
                if abs_pos in variants:
                    ref, alt = variants[abs_pos]
                    if len(ref) == 1 and len(alt) == 1:  # SNP
                        read_bases[rel_pos] = alt

        read_seq = ''.join([c if c in 'ACGTacgt' else 'A' for c in read_bases])

        a = pysam.AlignedSegment()
        a.query_name = f"read_{i}"
        a.query_sequence = read_seq
        a.flag = 0
        a.reference_id = 0
        a.reference_start = start
        a.mapping_quality = 60
        a.cigar = ((0, len(read_seq)),)
        a.query_qualities = [30] * len(read_seq)
        out_bam.write(a)

print("⚙️ Сортировка и индексирование BAM...")
sorted_bam = "data/giab_chr21/chr21_slice_sorted.bam"
pysam.sort("-o", sorted_bam, bam_path)
pysam.index(sorted_bam)

os.replace(sorted_bam, bam_path)
os.replace(sorted_bam + ".bai", bam_path + ".bai")

print("✅ Валидный BAM файл с VAF ~ 0.50 готов!")
