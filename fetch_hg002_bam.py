import os
import pysam

# Официальное открытое зеркало GIAB на Google Cloud Storage (без авторизации)
BAM_URL = "https://storage.googleapis.com/genomics-public-data/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/latest/GRCh38/HG002_GRCh38_chr21.bam"

OUTPUT_BAM = "data/giab_chr21/chr21_slice.bam"
OUTPUT_BAI = "data/giab_chr21/chr21_slice.bam.bai"

# Вырезаем окно в 500,000 нуклеотидов (~15 МБ данных)
REGION_CONTIG = "chr21"
START_POS = 10_000_000
END_POS = 10_500_000

print(f"📥 Подключение к публичному хранилищу GIAB (Google Cloud)...")
print(f"🧬 Извлечение реальных прочтений HG002 для {REGION_CONTIG}:{START_POS}-{END_POS}...")

try:
    # htslib стримит по сети только фрагмент байтов для нужного региона
    remote_bam = pysam.AlignmentFile(BAM_URL, "rb")

    count = 0
    temp_bam = "data/giab_chr21/temp_slice.bam"

    with pysam.AlignmentFile(temp_bam, "wb", header=remote_bam.header) as out_bam:
        for read in remote_bam.fetch(REGION_CONTIG, START_POS, END_POS):
            # Сдвигаем координаты на -10,000,000, чтобы привязать к началу локального FASTA (позиция 1)
            read.reference_start = read.reference_start - START_POS
            if read.next_reference_start > 0:
                read.next_reference_start = max(0, read.next_reference_start - START_POS)

            if read.reference_start >= 0:
                out_bam.write(read)
                count += 1

    print(f"⚙️ Извлечено {count} реальных Illumina ридов. Сортировка и индексирование...")

    # Сортируем и индексируем вырезанный BAM
    pysam.sort("-o", OUTPUT_BAM, temp_bam)
    pysam.index(OUTPUT_BAM)

    if os.path.exists(temp_bam):
        os.remove(temp_bam)

    size_mb = os.path.getsize(OUTPUT_BAM) / (1024 * 1024)
    print(f"✅ Успешно! Файл {OUTPUT_BAM} готов (Размер: {size_mb:.2f} МБ).")

except Exception as e:
    print(f"❌ Ошибка выкачки: {e}")
