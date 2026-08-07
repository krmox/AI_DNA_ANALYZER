import inspect
import sys
import torch
from torch.utils.data import DataLoader

try:
    from providers import GiabAlignmentProvider, ProviderDataset
except ImportError as e:
    print(f"❌ Ошибка импорта: {e}")
    sys.exit(1)

FASTA_PATH = "data/giab_chr21/chr21_slice.fa"
VCF_PATH = "data/giab_chr21/chr21_slice.vcf.gz"
BAM_PATH = "data/giab_chr21/chr21_slice.bam"
CONTIG = "chr21"

def build_provider():
    """Адаптивно собирает провайдер под сигнатуру __init__"""
    sig = inspect.signature(GiabAlignmentProvider.__init__)
    params = sig.parameters
    kwargs = {}

    # Сопоставляем пути к файлам
    for key in ['fasta_path', 'fasta']:
        if key in params: kwargs[key] = FASTA_PATH
    for key in ['vcf_path', 'vcf']:
        if key in params: kwargs[key] = VCF_PATH
    for key in ['bam_path', 'bam']:
        if key in params: kwargs[key] = BAM_PATH

    if 'contig' in params: kwargs['contig'] = CONTIG

    for key in ['start', 'region_start', 'start_pos']:
        if key in params: kwargs[key] = 1
    for key in ['end', 'region_end', 'end_pos']:
        if key in params: kwargs[key] = 5000000

    return GiabAlignmentProvider(**kwargs)

def test_pipeline():
    print("🚀 Инициализация GiabAlignmentProvider & ProviderDataset...")

    provider = build_provider()

    # Открываем контекст провайдера для чтения htslib/pysam
    with provider:
        dataset = ProviderDataset(provider)
        total_samples = len(dataset)
        print(f"📊 Всего сэмплов (окон) в датасете: {total_samples}")

        if total_samples == 0:
            print("❌ Ошибка: Датасет пуст! Проверь границы региона или нарезку окон.")
            return

        BATCH_SIZE = 8
        loader = DataLoader(dataset, batch_size=BATCH_SIZE, shuffle=True, num_workers=0)

        print(f"⚙️ Загрузка первого батча (batch_size={BATCH_SIZE})...\n")

        batch = next(iter(loader))

        print("==================================================")
        print("ПРОВЕРКА ТЕНЗОРОВ И СТРУКТУРЫ БАТЧА")
        print("==================================================")

        if isinstance(batch, dict):
            for k, v in batch.items():
                if isinstance(v, torch.Tensor):
                    has_nan = torch.isnan(v).any().item() or torch.isinf(v).any().item()
                    print(f"  • Key: {k:<15} | Shape: {str(list(v.shape)):<18} | Dtype: {v.dtype} | Has NaN/Inf: {has_nan}")
                else:
                    print(f"  • Key: {k:<15} | Non-tensor value: {type(v)}")
        else:
            print(f"  • Raw batch type: {type(batch)}")

        # Проверим доступность CUDA для Mamba
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        print("\n==================================================")
        print(f"ПРОВЕРКА ТРАНСФЕРА НА DEVICE ({device})")
        print("==================================================")
        if isinstance(batch, dict):
            for k, v in batch.items():
                if isinstance(v, torch.Tensor):
                    v_dev = v.to(device)
                    print(f"  ✅ {k} -> {v_dev.device}")

        print("\n✅ DataLoader работает корректно! Тензоры готовы к передаче в модель.")

if __name__ == "__main__":
    test_pipeline()
