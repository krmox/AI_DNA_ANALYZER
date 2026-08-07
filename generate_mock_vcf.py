import random
import pysam

vcf_path = "data/giab_chr21/chr21_slice.vcf"
vcf_gz_path = "data/giab_chr21/chr21_slice.vcf.gz"

header = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr21,length=5000001>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002
"""

bases = ['A', 'C', 'G', 'T']
lines = [header]

# Генерируем ~5000 вариантов (плотность ~1 вариант на 1000 bp, как в реальном геноме)
for pos in range(1000, 5000000, 1000):
    rand = random.random()
    if rand < 0.8: # 80% SNP
        ref = random.choice(bases)
        alt = random.choice([b for b in bases if b != ref])
    elif rand < 0.9: # 10% Insertion
        ref = 'A'
        alt = 'AT'
    else: # 10% Deletion
        ref = 'AT'
        alt = 'A'
    lines.append(f"chr21\t{pos}\t.\t{ref}\t{alt}\t99\tPASS\t.\tGT\t0/1\n")

with open(vcf_path, "w") as f:
    f.writelines(lines)

pysam.tabix_compress(vcf_path, vcf_gz_path, force=True)
pysam.tabix_index(vcf_gz_path, preset="vcf", force=True)
print("✅ Синтетический VCF создан и заиндексирован!")
