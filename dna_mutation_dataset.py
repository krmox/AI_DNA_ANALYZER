"""
DNA Mutation Detection PoC — Data Generation Module
=====================================================

Модуль содержит:
  1. DNATokenizer          — простой токенизатор нуклеотидов A/C/G/T (+ N как pad/unknown)
  2. SyntheticDNADataset   — генератор синтетических ДНК-последовательностей
                              с внедрёнными мутациями (SNP / Insertion / Deletion)
                              для задачи посимвольной (token-level) классификации.

Зависимости: только torch (+ стандартная библиотека Python).
"""

from __future__ import annotations

import random
from typing import Dict, List

import torch
from torch.utils.data import Dataset, DataLoader


# ============================================================================
# 1. Токенизатор
# ============================================================================

class DNATokenizer:
    """
    Простой токенизатор для последовательностей ДНК.

    Алфавит:
        'N' -> 0   (padding / unknown — используется для заполнения после делеций)
        'A' -> 1
        'C' -> 2
        'G' -> 3
        'T' -> 4
    """

    # Токен -> id
    VOCAB: Dict[str, int] = {"N": 0, "A": 1, "C": 2, "G": 3, "T": 4}
    # id -> Токен (обратный словарь)
    INV_VOCAB: Dict[int, str] = {v: k for k, v in VOCAB.items()}

    # Нуклеотиды, из которых генерируется "настоящая" последовательность
    # (N используется только как служебный/паддинг символ, случайно не генерируется)
    NUCLEOTIDES: List[str] = ["A", "C", "G", "T"]

    def __init__(self) -> None:
        self.vocab_size = len(self.VOCAB)

    def encode(self, sequence_str: str) -> torch.Tensor:
        """
        Преобразует строку нуклеотидов в тензор индексов (LongTensor).

        Неизвестные символы (не входящие в алфавит) кодируются как 'N' (0).

        Args:
            sequence_str: строка вида "ACGTAG..."

        Returns:
            torch.LongTensor формы [len(sequence_str)]
        """
        ids = [self.VOCAB.get(ch, self.VOCAB["N"]) for ch in sequence_str]
        return torch.tensor(ids, dtype=torch.long)

    def decode(self, tensor: torch.Tensor) -> str:
        """
        Преобразует тензор индексов обратно в строку нуклеотидов.

        Args:
            tensor: 1D тензор (или список) индексов токенов

        Returns:
            Строка нуклеотидов, например "ACGTNN..."
        """
        if isinstance(tensor, torch.Tensor):
            ids = tensor.tolist()
        else:
            ids = list(tensor)
        return "".join(self.INV_VOCAB.get(int(i), "N") for i in ids)


# ============================================================================
# 2. Синтетический датасет с мутациями
# ============================================================================

class SyntheticDNADataset(Dataset):
    """
    Генерирует случайные ДНК-последовательности фиксированной длины и
    внедряет в них мутации, размечая каждую позицию классом мутации.

    Классы меток (per-token classification):
        0 — Normal      (нормальный нуклеотид, мутаций нет)
        1 — SNP         (Single Nucleotide Polymorphism — точечная замена)
        2 — Insertion   (вставка "лишнего" нуклеотида)
        3 — Deletion    (удаление нуклеотида — позиция помечается как делеция,
                          а сама последовательность "схлопывается" и справа
                          добивается символом 'N', чтобы сохранить длину seq_len)

    Логика генерации одного сэмпла:
        1. Генерируется случайная "эталонная" последовательность длины seq_len.
        2. Для каждой позиции с вероятностью `mutation_rate` выбирается один
           из трёх типов мутации (равновероятно между SNP/Insertion/Deletion).
        3. Мутации применяются последовательно (слева направо) к рабочему
           списку нуклеотидов и соответствующему списку меток:
             - SNP: нуклеотид в текущей позиции заменяется на другой
               (гарантированно отличный от исходного) -> label = 1.
             - Insertion: перед текущим нуклеотидом вставляется новый
               случайный нуклеотид -> вставленная позиция получает label = 2.
               Оригинальный нуклеотид, стоявший на этом месте, сдвигается
               вправо и сохраняет свой label (обычно 0).
             - Deletion: текущий нуклеотид удаляется из последовательности
               -> label = 3 ставится на позицию, где произошло удаление
               (т.е. на нуклеотид, "занявший" освободившееся место).
        4. После применения мутаций последовательность обрезается или
           дополняется символом 'N' (label = 0, т.к. паддинг — не мутация)
           до фиксированной длины seq_len.

    Итоговый семпл — это словарь:
        {
            "input_ids": LongTensor [seq_len],
            "labels":    LongTensor [seq_len],
        }
    """

    # Метки классов мутаций
    LABEL_NORMAL = 0
    LABEL_SNP = 1
    LABEL_INSERTION = 2
    LABEL_DELETION = 3

    def __init__(
        self,
        num_samples: int,
        seq_len: int,
        mutation_rate: float = 0.05,
        seed: int | None = None,
    ) -> None:
        """
        Args:
            num_samples:   количество сэмплов в датасете (виртуальный размер).
            seq_len:       длина каждой генерируемой последовательности (в токенах).
            mutation_rate: вероятность того, что в конкретной позиции референсной
                            последовательности возникнет мутация (0.0 - 1.0).
            seed:          опциональный seed для воспроизводимости.
        """
        if not (0.0 <= mutation_rate <= 1.0):
            raise ValueError("mutation_rate должен быть в диапазоне [0, 1]")

        self.num_samples = num_samples
        self.seq_len = seq_len
        self.mutation_rate = mutation_rate
        self.tokenizer = DNATokenizer()

        self._rng = random.Random(seed)

    def __len__(self) -> int:
        return self.num_samples

    def _random_nucleotide(self, exclude: str | None = None) -> str:
        """Возвращает случайный нуклеотид, опционально исключая заданный."""
        choices = self.tokenizer.NUCLEOTIDES
        if exclude is not None:
            choices = [n for n in choices if n != exclude]
        return self._rng.choice(choices)

    def _generate_reference(self) -> str:
        """Генерирует случайную эталонную последовательность длины seq_len."""
        return "".join(self._rng.choice(self.tokenizer.NUCLEOTIDES) for _ in range(self.seq_len))

    def _apply_mutations(self, reference: str) -> tuple[List[str], List[int]]:
        """
        Проходит по эталонной последовательности слева направо и с вероятностью
        mutation_rate внедряет в текущую позицию одну из мутаций (SNP / Insertion
        / Deletion), формируя итоговую (рабочую) последовательность и метки.

        Returns:
            (nucleotides, labels) — списки одинаковой длины (до финального
            выравнивания по seq_len).
        """
        nucleotides: List[str] = []
        labels: List[int] = []

        for ref_base in reference:
            if self._rng.random() < self.mutation_rate:
                mutation_type = self._rng.choice(["snp", "insertion", "deletion"])

                if mutation_type == "snp":
                    # Замена нуклеотида на другой (гарантированно отличный)
                    nucleotides.append(self._random_nucleotide(exclude=ref_base))
                    labels.append(self.LABEL_SNP)

                elif mutation_type == "insertion":
                    # Вставляем лишний случайный нуклеотид перед оригинальным
                    nucleotides.append(self._random_nucleotide())
                    labels.append(self.LABEL_INSERTION)
                    # Затем сохраняем оригинальный нуклеотид как нормальный
                    nucleotides.append(ref_base)
                    labels.append(self.LABEL_NORMAL)

                else:  # deletion
                    # Оригинальный нуклеотид "удалён" — просто не добавляем его
                    # в последовательность. Метку деляции ставим на позицию
                    # непосредственно перед пропуском (сигнал модели о разрыве).
                    if nucleotides:
                        labels[-1] = self.LABEL_DELETION
                    else:
                        # Если делеция — самый первый символ, помечаем следующий
                        # добавленный нуклеотид как содержащий делецию слева.
                        # Практически: вставим служебный маркер отдельно.
                        nucleotides.append(ref_base if False else "")  # no-op placeholder
                        # Убираем placeholder — делецию в самой первой позиции
                        # просто пропускаем без метки (крайне редкий edge-case).
                        nucleotides.pop()
                    # Нуклеотид НЕ добавляется — это и есть "удаление".

            else:
                nucleotides.append(ref_base)
                labels.append(self.LABEL_NORMAL)

        return nucleotides, labels

    def _pad_or_truncate(
        self, nucleotides: List[str], labels: List[int]
    ) -> tuple[List[str], List[int]]:
        """Приводит последовательность и метки к фиксированной длине seq_len."""
        if len(nucleotides) > self.seq_len:
            nucleotides = nucleotides[: self.seq_len]
            labels = labels[: self.seq_len]
        elif len(nucleotides) < self.seq_len:
            pad_len = self.seq_len - len(nucleotides)
            nucleotides = nucleotides + ["N"] * pad_len
            labels = labels + [self.LABEL_NORMAL] * pad_len
        return nucleotides, labels

    def __getitem__(self, idx: int) -> Dict[str, torch.Tensor]:
        """
        Генерирует один сэмпл "на лету" (данные не кешируются, поэтому при
        каждом обращении к одному и тому же idx может вернуться разная
        последовательность — это нормально для синтетического генератора).

        Returns:
            dict с ключами:
              'input_ids': LongTensor [seq_len]
              'labels':    LongTensor [seq_len]
        """
        reference = self._generate_reference()
        nucleotides, labels = self._apply_mutations(reference)
        nucleotides, labels = self._pad_or_truncate(nucleotides, labels)

        sequence_str = "".join(nucleotides)
        input_ids = self.tokenizer.encode(sequence_str)
        labels_tensor = torch.tensor(labels, dtype=torch.long)

        return {"input_ids": input_ids, "labels": labels_tensor}


# ============================================================================
# 3. Блок проверки
# ============================================================================

if __name__ == "__main__":
    LABEL_NAMES = {0: "Normal", 1: "SNP", 2: "Insertion", 3: "Deletion"}

    tokenizer = DNATokenizer()
    dataset = SyntheticDNADataset(
        num_samples=100,
        seq_len=32,
        mutation_rate=0.05,
        seed=42,
    )

    dataloader = DataLoader(dataset, batch_size=4, shuffle=True)

    batch = next(iter(dataloader))
    input_ids = batch["input_ids"]
    labels = batch["labels"]

    print("=" * 70)
    print("ПРОВЕРКА DataLoader")
    print("=" * 70)
    print(f"input_ids.shape: {tuple(input_ids.shape)}")
    print(f"labels.shape:    {tuple(labels.shape)}")
    print()

    print("Пример декодированных последовательностей батча:")
    print("-" * 70)
    for i in range(input_ids.shape[0]):
        seq_str = tokenizer.decode(input_ids[i])
        label_row = labels[i].tolist()
        label_str = " ".join(str(l) for l in label_row)

        print(f"[{i}] Sequence : {seq_str}")
        print(f"    Labels   : {label_str}")

        mutated_positions = [
            (pos, LABEL_NAMES[l]) for pos, l in enumerate(label_row) if l != 0
        ]
        if mutated_positions:
            details = ", ".join(f"pos={p}:{name}" for p, name in mutated_positions)
            print(f"    Mutations: {details}")
        else:
            print(f"    Mutations: (нет мутаций в этом сэмпле)")
        print()
