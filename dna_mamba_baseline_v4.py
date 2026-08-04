"""
DNA Mutation Detection PoC — Mamba (SSM) Baseline v4: Reference-Aware Input
=============================================================================

КОРНЕВОЙ ФИКС по итогам диагностики v3: до этого момента модель училась
находить мутации, глядя ТОЛЬКО на итоговую (уже мутировавшую) наблюдаемую
последовательность — без референса. Поскольку каждая буква (что на месте
SNP, что на нормальной позиции) генерируется как равновероятный нуклеотид,
а соседние позиции независимы (i.i.d.), задача была статистически
неразрешима: SNP и Normal неотличимы без референса для сравнения. Это и
объясняло `f1_mutations_macro=0.0000`, стабильно воспроизводившееся при
любых гиперпараметрах.

ЧТО ИЗМЕНИЛОСЬ:

  1. Датасет (`SyntheticDNADataset.__getitem__`) теперь возвращает ДВА
     выровненных по позициям канала одинаковой длины seq_len:
       - 'reference_ids' — эталонная последовательность (без мутаций).
       - 'input_ids'     — наблюдаемая (потенциально мутировавшая)
                            последовательность.
     Плюс, как и раньше, 'labels' — класс события в каждой позиции.

     Чтобы сохранить строгое 1:1 позиционное выравнивание (без сдвигов
     длины, которые в v1-v3 разрушали смысл per-position классификации),
     Insertion и Deletion переопределены как события "в рамках одной и той
     же координаты референса" (по аналогии с пайлап-представлением в
     variant calling / SAM-CIGAR, где делеция в риде относительно
     референса представляется как gap на выровненной позиции):
       - Normal:    observed[i] = reference[i]                    label=0
       - SNP:       observed[i] = случайная буква != reference[i]  label=1
                     (гарантированное несовпадение с референсом)
       - Insertion: observed[i] = случайная буква (референс НЕ
                     исключается — может случайно совпасть, что
                     соответствует, например, дупликации основания)
                                                                     label=2
       - Deletion:  observed[i] = 'N' (пропуск/gap — базовый вызов
                     отсутствует на этой выровненной позиции)        label=3
     Референс на этой позиции всегда содержит исходную "правильную" букву
     независимо от типа события — модель должна научиться сравнивать
     observed[i] с reference[i], чтобы обнаружить несовпадение и
     классифицировать его тип.

     ЧЕСТНАЯ ОГОВОРКА: это упрощение жертвует биологической точностью
     истинных inDel-событий (которые в реальности сдвигают длину
     последовательности) ради того, чтобы задача была вообще разрешима в
     формате строгой per-position классификации без сложных механизмов
     выравнивания (attention/alignment-слоёв, которые явно запрещены
     текущим заданием). Это стандартное упрощение для pileup-style PoC.

  2. Модель (`DNAMambaBaseline`) получила ДВА отдельных Embedding-слоя
     (для reference и для observed) и слой слияния (fusion):
       - 'concat' (по умолчанию): конкатенация эмбеддингов по последней
         размерности -> Linear(2*d_model -> d_model). Даёт модели
         возможность выучить произвольную (в т.ч. похожую на разностную)
         функцию сравнения двух каналов.
       - 'sum': поэлементное сложение эмбеддингов (дешевле, менее
         выразительно).
     Никаких новых блоков/attention не добавлено — расширен только входной
     слой перед стеком MambaBlock, сам MambaBlock не изменён.

Зависимости: только torch (+ стандартная библиотека Python).
"""

from __future__ import annotations

import copy
import random
from typing import Dict, List, Optional, Tuple

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader


# ============================================================================
# 1. Токенизатор ДНК
# ============================================================================

class DNATokenizer:
    """
    Токенизатор нуклеотидов.

    Алфавит:
        'N' -> 0   (padding / unknown / gap-маркер делеции в observed-канале)
        'A' -> 1
        'C' -> 2
        'G' -> 3
        'T' -> 4
    """

    VOCAB: Dict[str, int] = {"N": 0, "A": 1, "C": 2, "G": 3, "T": 4}
    INV_VOCAB: Dict[int, str] = {v: k for k, v in VOCAB.items()}
    NUCLEOTIDES: List[str] = ["A", "C", "G", "T"]

    def __init__(self) -> None:
        self.vocab_size = len(self.VOCAB)

    def encode(self, sequence_str: str) -> torch.Tensor:
        ids = [self.VOCAB.get(ch, self.VOCAB["N"]) for ch in sequence_str]
        return torch.tensor(ids, dtype=torch.long)

    def decode(self, tensor: torch.Tensor) -> str:
        ids = tensor.tolist() if isinstance(tensor, torch.Tensor) else list(tensor)
        return "".join(self.INV_VOCAB.get(int(i), "N") for i in ids)


# ============================================================================
# 2. Синтетический датасет: reference_ids + input_ids (ИСПРАВЛЕНИЕ: reference-aware)
# ============================================================================

class SyntheticDNADataset(Dataset):
    """
    Генерирует ПАРУ выровненных по позициям последовательностей длины
    seq_len — референс и наблюдаемую (потенциально мутировавшую) — плюс
    метку класса события в каждой позиции.

    Метки:
        0 — Normal:    observed[i] == reference[i]
        1 — SNP:       observed[i] = случайная буква != reference[i]
        2 — Insertion: observed[i] = случайная буква (reference не исключается)
        3 — Deletion:  observed[i] = 'N' (gap на выровненной позиции)
    """

    LABEL_NORMAL = 0
    LABEL_SNP = 1
    LABEL_INSERTION = 2
    LABEL_DELETION = 3

    def __init__(
        self,
        num_samples: int,
        seq_len: int,
        mutation_rate: float = 0.05,
        seed: Optional[int] = None,
    ) -> None:
        if not (0.0 <= mutation_rate <= 1.0):
            raise ValueError("mutation_rate должен быть в диапазоне [0, 1]")

        self.num_samples = num_samples
        self.seq_len = seq_len
        self.mutation_rate = mutation_rate
        self.tokenizer = DNATokenizer()
        self._rng = random.Random(seed)

    def __len__(self) -> int:
        return self.num_samples

    def _random_nucleotide(self, exclude: Optional[str] = None) -> str:
        choices = self.tokenizer.NUCLEOTIDES
        if exclude is not None:
            choices = [n for n in choices if n != exclude]
        return self._rng.choice(choices)

    def _generate_reference(self) -> str:
        return "".join(self._rng.choice(self.tokenizer.NUCLEOTIDES) for _ in range(self.seq_len))

    def _generate_observed(self, reference: str) -> Tuple[str, List[int]]:
        """
        Строит observed-последовательность и метки строго позиция-в-позицию
        относительно reference (без изменения длины — см. docstring модуля).
        """
        observed_chars: List[str] = []
        labels: List[int] = []

        for ref_base in reference:
            if self._rng.random() < self.mutation_rate:
                mutation_type = self._rng.choice(["snp", "insertion", "deletion"])

                if mutation_type == "snp":
                    observed_chars.append(self._random_nucleotide(exclude=ref_base))
                    labels.append(self.LABEL_SNP)

                elif mutation_type == "insertion":
                    observed_chars.append(self._random_nucleotide())  # exclude=None
                    labels.append(self.LABEL_INSERTION)

                else:  # deletion
                    observed_chars.append("N")
                    labels.append(self.LABEL_DELETION)

            else:
                observed_chars.append(ref_base)
                labels.append(self.LABEL_NORMAL)

        return "".join(observed_chars), labels

    def __getitem__(self, idx: int) -> Dict[str, torch.Tensor]:
        reference = self._generate_reference()
        observed, labels = self._generate_observed(reference)

        reference_ids = self.tokenizer.encode(reference)
        input_ids = self.tokenizer.encode(observed)
        labels_tensor = torch.tensor(labels, dtype=torch.long)

        return {
            "reference_ids": reference_ids,
            "input_ids": input_ids,
            "labels": labels_tensor,
        }


# ============================================================================
# 3. Bidirectional Causal Gated-Convolutional PoC блок ("Bi-Mamba") — без изменений
# ============================================================================

class CausalGatedConvCore(nn.Module):
    """
    Однонаправленное каузальное "ядро" (строительный блок для
    BidirectionalMambaBlock): LayerNorm -> Depthwise Causal Conv1d ->
    LayerNorm -> SiLU -> Selective Gate. Residual/комбинация направлений —
    на уровне BidirectionalMambaBlock.
    """

    def __init__(self, d_model: int, conv_kernel_size: int = 3, dropout: float = 0.15) -> None:
        super().__init__()
        self.norm_pre = nn.LayerNorm(d_model)

        self.conv_padding = conv_kernel_size - 1
        self.depthwise_conv = nn.Conv1d(
            in_channels=d_model,
            out_channels=d_model,
            kernel_size=conv_kernel_size,
            padding=self.conv_padding,
            groups=d_model,
        )

        self.norm_post_conv = nn.LayerNorm(d_model)
        self.activation = nn.SiLU()
        self.selective_gate = nn.Linear(d_model, d_model)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x_norm = self.norm_pre(x)

        x_conv_in = x_norm.transpose(1, 2)
        x_conv_out = self.depthwise_conv(x_conv_in)
        x_conv_out = x_conv_out[:, :, : x_norm.size(1)]
        x_conv_out = x_conv_out.transpose(1, 2)

        x_conv_out = self.norm_post_conv(x_conv_out)
        x_conv_out = self.activation(x_conv_out)
        x_conv_out = self.dropout(x_conv_out)

        gate = torch.sigmoid(self.selective_gate(x_norm))
        h = gate * x_conv_out + (1.0 - gate) * x_norm
        return h


class MambaBlock(nn.Module):
    """
    Bidirectional Causal Gated-Convolutional PoC блок ("Bi-Mamba").
    Не изменён относительно v3 — см. предыдущую версию за подробным
    докстрингом. forward_core (левый контекст) + backward_core (правый
    контекст через flip) -> concat -> Linear(2D->D) -> residual.
    """

    def __init__(self, d_model: int, conv_kernel_size: int = 3, dropout: float = 0.15) -> None:
        super().__init__()
        self.forward_core = CausalGatedConvCore(d_model, conv_kernel_size, dropout)
        self.backward_core = CausalGatedConvCore(d_model, conv_kernel_size, dropout)

        self.combine_proj = nn.Linear(2 * d_model, d_model)
        self.combine_dropout = nn.Dropout(dropout)
        self.out_norm = nn.LayerNorm(d_model)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        residual = x

        h_forward = self.forward_core(x)

        x_reversed = torch.flip(x, dims=[1])
        h_backward_rev = self.backward_core(x_reversed)
        h_backward = torch.flip(h_backward_rev, dims=[1])

        combined = torch.cat([h_forward, h_backward], dim=-1)
        combined = self.combine_proj(combined)
        combined = self.combine_dropout(combined)

        out = residual + combined
        out = self.out_norm(out)
        return out


# ============================================================================
# 4. Полная модель — расширенный входной слой (reference + observed)
# ============================================================================

class DNAMambaBaseline(nn.Module):
    """
    Baseline-модель для token-level классификации мутаций в ДНК,
    ОСОЗНАЮЩАЯ РЕФЕРЕНС (reference-aware).

    Pipeline:
        input_ids [B, L], reference_ids [B, L]  (long)
          -> Embedding(observed) [B, L, D]  +  Embedding(reference) [B, L, D]
          -> Fusion ('concat'+Linear ИЛИ 'sum')  -> [B, L, D]
          -> N x MambaBlock (Bidirectional Causal Gated-Conv PoC) [B, L, D]
          -> LayerNorm -> Classification Head (Linear D -> num_classes)

    Два отдельных Embedding-слоя (не общие веса) для reference/observed,
    т.к. одна и та же буква играет разную "роль" в зависимости от канала —
    так модели проще выучить сравнение, чем если бы оба канала шарили
    один и тот же Embedding.
    """

    def __init__(
        self,
        vocab_size: int = 5,
        d_model: int = 128,
        num_layers: int = 6,
        num_classes: int = 4,
        conv_kernel_size: int = 3,
        dropout: float = 0.15,
        fusion: str = "concat",
    ) -> None:
        super().__init__()

        if fusion not in ("concat", "sum"):
            raise ValueError("fusion должен быть 'concat' или 'sum'")
        self.fusion = fusion

        self.observed_embedding = nn.Embedding(num_embeddings=vocab_size, embedding_dim=d_model, padding_idx=0)
        self.reference_embedding = nn.Embedding(num_embeddings=vocab_size, embedding_dim=d_model, padding_idx=0)

        if fusion == "concat":
            self.fusion_proj: Optional[nn.Linear] = nn.Linear(2 * d_model, d_model)
        else:
            self.fusion_proj = None

        self.embed_dropout = nn.Dropout(dropout)

        self.mamba_blocks = nn.ModuleList(
            [
                MambaBlock(d_model=d_model, conv_kernel_size=conv_kernel_size, dropout=dropout)
                for _ in range(num_layers)
            ]
        )

        self.final_norm = nn.LayerNorm(d_model)
        self.classifier = nn.Linear(d_model, num_classes)

    def forward(self, input_ids: torch.Tensor, reference_ids: torch.Tensor) -> torch.Tensor:
        """
        Args:
            input_ids:     [B, L] long tensor наблюдаемой последовательности
            reference_ids: [B, L] long tensor эталонной последовательности
        Returns:
            logits: [B, L, num_classes]
        """
        obs_emb = self.observed_embedding(input_ids)  # [B, L, D]
        ref_emb = self.reference_embedding(reference_ids)  # [B, L, D]

        if self.fusion == "concat":
            x = torch.cat([obs_emb, ref_emb], dim=-1)  # [B, L, 2D]
            x = self.fusion_proj(x)  # [B, L, D]
        else:  # 'sum'
            x = obs_emb + ref_emb  # [B, L, D]

        x = self.embed_dropout(x)

        for block in self.mamba_blocks:
            x = block(x)

        x = self.final_norm(x)
        logits = self.classifier(x)  # [B, L, num_classes]
        return logits


# ============================================================================
# 5. Focal Loss с калиброванными (sqrt + capped) весами классов
# ============================================================================

class FocalLoss(nn.Module):
    """Focal Loss (Lin et al., 2017) с per-class alpha-весами и label smoothing."""

    def __init__(
        self,
        alpha: Optional[torch.Tensor] = None,
        gamma: float = 2.0,
        label_smoothing: float = 0.0,
        reduction: str = "mean",
    ) -> None:
        super().__init__()
        self.alpha = alpha
        self.gamma = gamma
        self.label_smoothing = label_smoothing
        self.reduction = reduction

    def forward(self, logits: torch.Tensor, targets: torch.Tensor) -> torch.Tensor:
        num_classes = logits.size(-1)
        log_probs = F.log_softmax(logits, dim=-1)
        probs = log_probs.exp()

        if self.label_smoothing > 0.0:
            with torch.no_grad():
                true_dist = torch.full_like(log_probs, fill_value=self.label_smoothing / (num_classes - 1))
                true_dist.scatter_(1, targets.unsqueeze(1), 1.0 - self.label_smoothing)
            ce = -(true_dist * log_probs).sum(dim=1)
        else:
            ce = F.nll_loss(log_probs, targets, reduction="none")

        pt = probs.gather(1, targets.unsqueeze(1)).squeeze(1).clamp(min=1e-6, max=1.0)
        focal_factor = (1.0 - pt) ** self.gamma

        if self.alpha is not None:
            alpha_factor = self.alpha.to(logits.device)[targets]
        else:
            alpha_factor = torch.ones_like(pt)

        loss = alpha_factor * focal_factor * ce

        if self.reduction == "mean":
            return loss.mean()
        elif self.reduction == "sum":
            return loss.sum()
        return loss


def compute_calibrated_class_weights(
    dataset: SyntheticDNADataset,
    num_classes: int,
    sample_size: int = 200,
    max_weight_cap: float = 6.0,
) -> torch.Tensor:
    """sqrt-сглаженные, нормированные и capped веса классов."""
    counts = torch.zeros(num_classes, dtype=torch.float)

    for i in range(min(sample_size, len(dataset))):
        labels = dataset[i]["labels"]
        for c in range(num_classes):
            counts[c] += (labels == c).sum().item()

    counts = torch.clamp(counts, min=1.0)
    total = counts.sum()

    weights_raw = torch.sqrt(total / (num_classes * counts))
    weights_norm = weights_raw / weights_raw.mean()
    weights_final = torch.clamp(weights_norm, max=max_weight_cap)

    return weights_final


# ============================================================================
# 6. Метрики
# ============================================================================

def compute_f1_per_class(
    preds: torch.Tensor, targets: torch.Tensor, num_classes: int
) -> Dict[str, float]:
    metrics: Dict[str, float] = {}
    f1_scores: List[float] = []

    for c in range(num_classes):
        tp = ((preds == c) & (targets == c)).sum().item()
        fp = ((preds == c) & (targets != c)).sum().item()
        fn = ((preds != c) & (targets == c)).sum().item()

        precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

        metrics[f"f1_class_{c}"] = f1
        f1_scores.append(f1)

    metrics["macro_f1"] = sum(f1_scores) / len(f1_scores)
    mutation_f1s = [metrics[f"f1_class_{c}"] for c in range(1, num_classes)]
    metrics["f1_mutations_macro"] = sum(mutation_f1s) / len(mutation_f1s) if mutation_f1s else 0.0

    return metrics


# ============================================================================
# 7. EarlyStopping (мониторит val_loss, mode='min')
# ============================================================================

class EarlyStopping:
    def __init__(self, patience: int = 5, min_delta: float = 1e-4, mode: str = "max") -> None:
        if mode not in ("min", "max"):
            raise ValueError("mode должен быть 'min' или 'max'")
        self.patience = patience
        self.min_delta = min_delta
        self.mode = mode
        self.best_score: float = float("inf") if mode == "min" else -float("inf")
        self.best_state_dict: Optional[dict] = None
        self.counter: int = 0
        self.should_stop: bool = False

    def step(self, score: float, model: nn.Module) -> None:
        if self.mode == "max":
            improved = score > self.best_score + self.min_delta
        else:
            improved = score < self.best_score - self.min_delta

        if improved:
            self.best_score = score
            self.best_state_dict = copy.deepcopy(model.state_dict())
            self.counter = 0
        else:
            self.counter += 1
            if self.counter >= self.patience:
                self.should_stop = True

    def restore_best(self, model: nn.Module) -> None:
        if self.best_state_dict is not None:
            model.load_state_dict(self.best_state_dict)


# ============================================================================
# 8. Обучающий и валидационный циклы (теперь передают reference_ids)
# ============================================================================

def train_epoch(
    model: nn.Module,
    dataloader: DataLoader,
    optimizer: torch.optim.Optimizer,
    criterion: nn.Module,
    device: torch.device,
) -> float:
    model.train()
    total_loss = 0.0
    total_tokens = 0

    for batch in dataloader:
        input_ids = batch["input_ids"].to(device)
        reference_ids = batch["reference_ids"].to(device)
        labels = batch["labels"].to(device)

        optimizer.zero_grad()

        logits = model(input_ids, reference_ids)
        loss = criterion(logits.reshape(-1, logits.size(-1)), labels.reshape(-1))

        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
        optimizer.step()

        batch_tokens = labels.numel()
        total_loss += loss.item() * batch_tokens
        total_tokens += batch_tokens

    return total_loss / total_tokens


@torch.no_grad()
def validate(
    model: nn.Module,
    dataloader: DataLoader,
    criterion: nn.Module,
    device: torch.device,
    num_classes: int,
) -> Tuple[float, Dict[str, float]]:
    model.eval()
    total_loss = 0.0
    total_tokens = 0

    all_preds: List[torch.Tensor] = []
    all_targets: List[torch.Tensor] = []

    for batch in dataloader:
        input_ids = batch["input_ids"].to(device)
        reference_ids = batch["reference_ids"].to(device)
        labels = batch["labels"].to(device)

        logits = model(input_ids, reference_ids)
        loss = criterion(logits.reshape(-1, logits.size(-1)), labels.reshape(-1))

        batch_tokens = labels.numel()
        total_loss += loss.item() * batch_tokens
        total_tokens += batch_tokens

        preds = torch.argmax(logits, dim=-1)
        all_preds.append(preds.reshape(-1).cpu())
        all_targets.append(labels.reshape(-1).cpu())

    val_loss = total_loss / total_tokens

    preds_flat = torch.cat(all_preds)
    targets_flat = torch.cat(all_targets)
    metrics = compute_f1_per_class(preds_flat, targets_flat, num_classes=num_classes)

    return val_loss, metrics


def build_warmup_cosine_scheduler(
    optimizer: torch.optim.Optimizer, warmup_epochs: int, total_epochs: int
) -> torch.optim.lr_scheduler.SequentialLR:
    warmup_epochs = max(1, min(warmup_epochs, total_epochs - 1))

    warmup_scheduler = torch.optim.lr_scheduler.LinearLR(
        optimizer, start_factor=0.1, end_factor=1.0, total_iters=warmup_epochs
    )
    cosine_scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimizer, T_max=max(1, total_epochs - warmup_epochs)
    )
    return torch.optim.lr_scheduler.SequentialLR(
        optimizer, schedulers=[warmup_scheduler, cosine_scheduler], milestones=[warmup_epochs]
    )


# ============================================================================
# 9. Блок запуска
# ============================================================================

if __name__ == "__main__":
    SEQ_LEN = 64
    NUM_CLASSES = 4
    D_MODEL = 128
    NUM_LAYERS = 6
    DROPOUT = 0.15
    BATCH_SIZE = 16
    NUM_EPOCHS = 30
    WARMUP_EPOCHS = 3
    LR = 1e-3
    # Возвращаем реалистичный mutation_rate (был временно поднят до 0.18
    # только для отладки в v3) — теперь, когда добавлен референс, сигнал
    # должен быть детектируемым и на исходном разреженном режиме.
    MUTATION_RATE = 0.05
    FOCAL_GAMMA = 2.0
    LABEL_SMOOTHING = 0.0
    MAX_WEIGHT_CAP = 6.0
    EARLY_STOPPING_PATIENCE = 5
    FUSION = "concat"  # 'concat' (Linear(2D->D)) или 'sum'

    LABEL_NAMES = {0: "Normal", 1: "SNP", 2: "Insertion", 3: "Deletion"}

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Устройство: {device}")

    torch.manual_seed(42)

    tokenizer = DNATokenizer()

    train_dataset = SyntheticDNADataset(
        num_samples=1000, seq_len=SEQ_LEN, mutation_rate=MUTATION_RATE, seed=None
    )
    val_dataset = SyntheticDNADataset(
        num_samples=200, seq_len=SEQ_LEN, mutation_rate=MUTATION_RATE, seed=123
    )

    train_loader = DataLoader(train_dataset, batch_size=BATCH_SIZE, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=BATCH_SIZE, shuffle=False)

    model = DNAMambaBaseline(
        vocab_size=tokenizer.vocab_size,
        d_model=D_MODEL,
        num_layers=NUM_LAYERS,
        num_classes=NUM_CLASSES,
        dropout=DROPOUT,
        fusion=FUSION,
    ).to(device)

    num_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"Модель DNAMambaBaseline (v4, Reference-Aware, fusion={FUSION}) создана. Параметров: {num_params:,}")

    class_weights = compute_calibrated_class_weights(
        train_dataset, num_classes=NUM_CLASSES, sample_size=200, max_weight_cap=MAX_WEIGHT_CAP
    ).to(device)
    print(f"Калиброванные веса классов (sqrt + cap={MAX_WEIGHT_CAP}): {class_weights.tolist()}")

    criterion = FocalLoss(
        alpha=class_weights, gamma=FOCAL_GAMMA, label_smoothing=LABEL_SMOOTHING, reduction="mean"
    )

    optimizer = torch.optim.AdamW(model.parameters(), lr=LR)
    scheduler = build_warmup_cosine_scheduler(optimizer, warmup_epochs=WARMUP_EPOCHS, total_epochs=NUM_EPOCHS)
    early_stopping = EarlyStopping(patience=EARLY_STOPPING_PATIENCE, mode="min")

    print("\n" + "=" * 78)
    print("ОБУЧЕНИЕ (v4, Reference-Aware)")
    print("=" * 78)

    for epoch in range(1, NUM_EPOCHS + 1):
        train_loss = train_epoch(model, train_loader, optimizer, criterion, device)
        val_loss, val_metrics = validate(model, val_loader, criterion, device, num_classes=NUM_CLASSES)

        current_lr = optimizer.param_groups[0]["lr"]
        scheduler.step()

        per_class_str = ", ".join(
            f"F1[{LABEL_NAMES[c]}]={val_metrics[f'f1_class_{c}']:.3f}" for c in range(NUM_CLASSES)
        )

        print(
            f"Epoch {epoch:2d}/{NUM_EPOCHS} | LR: {current_lr:.6f} | "
            f"Train Loss: {train_loss:.4f} | Val Loss: {val_loss:.4f} | "
            f"Macro F1: {val_metrics['macro_f1']:.4f} | "
            f"Mutations Macro F1: {val_metrics['f1_mutations_macro']:.4f}"
        )
        print(f"           {per_class_str}")

        early_stopping.step(val_loss, model)
        if early_stopping.should_stop:
            print(
                f"\n[EarlyStopping] Нет улучшения val_loss "
                f"{EARLY_STOPPING_PATIENCE} эпох подряд. Остановка на эпохе {epoch}. "
                f"Лучший val_loss: {early_stopping.best_score:.4f}"
            )
            break

    early_stopping.restore_best(model)
    print(f"\nВосстановлены веса лучшей модели (val_loss={early_stopping.best_score:.4f}).")

    print("\n" + "=" * 78)
    print("ПРИМЕР ПРЕДСКАЗАНИЯ (случайный валидационный сэмпл, лучшая модель)")
    print("=" * 78)

    model.eval()
    sample_idx = random.randint(0, len(val_dataset) - 1)
    sample = val_dataset[sample_idx]

    input_ids = sample["input_ids"].unsqueeze(0).to(device)
    reference_ids = sample["reference_ids"].unsqueeze(0).to(device)
    ground_truth = sample["labels"]

    with torch.no_grad():
        logits = model(input_ids, reference_ids)
        prediction = torch.argmax(logits, dim=-1).squeeze(0).cpu()

    decoded_ref = tokenizer.decode(sample["reference_ids"])
    decoded_obs = tokenizer.decode(sample["input_ids"])

    print(f"Reference:   {decoded_ref}")
    print(f"Observed:    {decoded_obs}")
    print(f"Ground Truth:{' '.join(str(l.item()) for l in ground_truth)}")
    print(f"Prediction:  {' '.join(str(p.item()) for p in prediction)}")

    print("\nПозиции с мутациями (Ground Truth vs Prediction):")
    any_flagged = False
    for pos in range(SEQ_LEN):
        gt = ground_truth[pos].item()
        pr = prediction[pos].item()
        if gt != 0 or pr != 0:
            any_flagged = True
            match_marker = "OK" if gt == pr else "MISMATCH"
            print(
                f"  pos={pos:3d} | ref={decoded_ref[pos]} obs={decoded_obs[pos]} | "
                f"GT={LABEL_NAMES[gt]:<10} | Pred={LABEL_NAMES[pr]:<10} | {match_marker}"
            )
    if not any_flagged:
        print("  (в этом сэмпле нет ни истинных мутаций, ни ложных срабатываний)")
