#!/usr/bin/env bash
set -euo pipefail

echo "🧬 [GIAB Setup] Подготовка окружения для хромосомы 21 (HG002 GRCh38)..."

# 1. Проверка необходимых утилит (Fedora)
REQUIRED_TOOLS=("samtools" "bcftools" "curl")
FOR_INSTALL=()

for tool in "${REQUIRED_TOOLS[@]}"; do
    if ! command -v "$tool" &> /dev/null; then
        FOR_INSTALL+=("$tool")
    fi
done

if [ ${#FOR_INSTALL[@]} -ne 0 ]; then
    echo "⚠️ Не найдены утилиты: ${FOR_INSTALL[*]}"
    echo "Установи их через dnf:"
    echo "  sudo dnf install -y samtools bcftools htslib"
    exit 1
fi

# 2. Создание структуры папок
DATA_DIR="data/giab_chr21"
mkdir -p "$DATA_DIR"
cd "$DATA_DIR"

echo "📁 Данные будут сохранены в: $(pwd)"

# Зеркало Ensembl / EBI (работает стабильнее UCSC/NCBI)
ENSEMBL_CHR21_URL="https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz"
VCF_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"

REGION="21:10000000-15000000" # У Ensembl префикс без 'chr', parser его сам нормализует!

echo "⏳ 1. Извлечение референса FASTA для $REGION из зеркала Ensembl..."
if [ ! -f "chr21_slice.fa" ]; then
    # Качаем с User-Agent и автоповторами
    curl -A "Mozilla/5.0" -C - --retry 5 --retry-connrefused -sSL "$ENSEMBL_CHR21_URL" -o chr21_full.fa.gz
    gunzip -f chr21_full.fa.gz

    samtools faidx chr21_full.fa

    # В Ensembl хромосома называется "21", а не "chr21"
    samtools faidx chr21_full.fa "$REGION" > chr21_slice.fa
    samtools faidx chr21_slice.fa

    # Чистим временные файлы
    rm -f chr21_full.fa chr21_full.fa.fai
    echo "✅ Референс готов."
else
    echo "ℹ️ Референс уже скачан."
fi

echo "⏳ 2. Извлечение истинных мутаций VCF (GIAB Benchmark)..."
if [ ! -f "chr21_slice.vcf.gz" ]; then
    # Пробуем вытащить регион, если упадет - вытащим весь VCF (он всего ~95K в сжатом виде для региона)
    bcftools view -r "chr21:10000000-15000000" "$VCF_URL" -O z -o chr21_slice.vcf.gz || \
    bcftools view -r "21:10000000-15000000" "$VCF_URL" -O z -o chr21_slice.vcf.gz

    bcftools index -t chr21_slice.vcf.gz
    echo "✅ VCF готов."
else
    echo "ℹ️ VCF уже скачан."
fi

echo "🎉 Все мини-файлы готовы для тестирования!"
echo "Файлы в $(pwd):"
ls -lh chr21_slice.*
