# Final Forensic Audit

Дата аудита: 2026-09-20. Режим: **audit-only** (код, пороги, verdict'ы, рукопись, данные не менялись; создан только этот файл).
Объект проверки: `research/PAPER_MANUSCRIPT.md` (+ `TABLES_FINAL.md`, `FIGURE_PLAN.md`, `RESEARCH_TIMELINE.md`, `LIMITATIONS_AND_OPEN_QUESTIONS.md`, `CLAIM_EVIDENCE_MAP.csv`, вторичные документы `research/*.md`).
Источники истины: сырые `results/bench_v13…v20/*.json`, `experimental/*/happy/*/happy.summary.csv`, `experimental/head_to_head/HEAD_TO_HEAD_RESULTS.csv`, `results/bench_v20/chr20_results.json`, код в рабочем дереве, тесты, `History/*_DEVLOG.md`.

**Что я реально пересчитал/перепроверил сам в этой сессии** (остальное сверено между документами и помечено как «по записи»):
- SHA-256 `cascade.py`, `method_c_regression.py`, `read_level_pileup.py`, `providers.py`, `pileup_counts.py`, `extract_bench_v12.py`, `native/*.c`, `bench_v20_chr20.py` — совпадают с `ANALYSIS_FREEZE.json`/`FROZEN_MANIFEST.json`; файлы основного дерева побайтно совпадают с `.claude/worktrees/opt-chr20` (`cmp`).
- chr20: `results/bench_v20/chr20_results.json` и оба `happy.summary.csv` (все числа §D ниже).
- hap.py HG002 chr21 (AI ×2) из `experimental/unified_happy/happy/*`; DeepVariant/Clair3/GATK из `HEAD_TO_HEAD_RESULTS.csv`.
- Тесты: `test_native_reads_equivalence.py` → **17 passed**; `test_bench_v12/test_cascade/test_cheap_router/test_robustness_benchmark/test_bench_v13` → **79 passed**.
- Наличие всех файлов, на которые ссылается `CLAIM_EVIDENCE_MAP.csv` (найдено 13 отсутствующих).
- Бутстрэпы v13–v19 я **не пересчитывал**; их числа сверены между рукописью, `TABLES_FINAL.md`, `CLAIM_EVIDENCE_MAP.csv` и (для v14/v15) `results/*/summary.md`, `History/14_DEVLOG.md §21`.

---

## 1. Executive verdict

Рукопись в её нынешнем виде **внутренне согласована в пределах своего периметра и честна к отрицательным результатам** (v14 = NEGATIVE, отзыв «never loses a SNP», withdrawn 65×/«17–80×», история coordinate bug и провала первой native-версии, M-1 с исправленным «~35 %» → 3.5 %, явное предупреждение о перекрытии U-H2 с областью подбора констант). Все hap.py-числа AI и внешних callers из вашего списка подтверждены сырыми файлами. **Однако рукопись написана до двух последних экспериментов и не содержит ни одного из них**: chr20 (целая хромосома, strict verdict **DEGRADED**, ΔF1 −5.6·10⁻⁵) и оптимизации `load_reads` (26.9×/35.1×, bottleneck → PB 93–96 %). Поэтому ряд утверждений теперь **OUTDATED или CONTRADICTED**: «no chromosome- or genome-scale run», «`load_reads` — no speedup», «1.22× end-to-end / 98.7 % на load_reads+PB», «seven experiments», диапазон routed 0.04–0.11 %, состав failure taxonomy, а во вторичных документах — «cascade preserves PB-only accuracy in every tested experiment». Кроме того, (i) 13 файлов-источников из `CLAIM_EVIDENCE_MAP.csv` физически отсутствуют в чекауте (потеряны вместе с `/tmp`), включая **сырые post-fix hap.py HG005 и cache, на которых стоят M-1 HG005 и состав rescue 13/129** — а в `experimental/stress_test/happy/` лежит только *инвалидированный* pre-fix прогон (F1 0.0085); (ii) два теста (`test_native_pileup_equivalence.py`, `test_hg005_extraction_integrity.py`) не существуют, поэтому «102/102», «12-category suite» и «11 regression tests» сейчас невоспроизводимы; (iii) `experimental/` содержит захардкоженные пути на worktree `opt-chr20`. Блокирующие проблемы — **не научные (данные не опровергнуты), а интеграционные и provenance**: рукопись нельзя замораживать, пока в неё не внесены chr20 и load_reads и пока не восстановлены/пересозданы утраченные источники.

---

## 2. Critical findings

### CF-1 — CRITICAL — chr20 отсутствует в рукописи; ряд утверждений им противоречится
- **Claim:** Abstract/§1/§3.1/§4.1/§5/§6: cascade «tracks/preserves» PB; «seven regional experiments»; «no chromosome- or genome-scale run for any sample»; «largest single span 12 Mb».
- **Evidence:** `research/CHR20_VALIDATION_REPORT.md`, `results/bench_v20/chr20_results.json`: весь chr20 (0–64,444,167; 56,266,816 loci extracted, 70,324 SNP в frame), verdict **DEGRADED** (locus CI [−1.049e−4, −1.41e−5], block CI [−1.08e−4, −7.5e−6], McNemar 9:1 p=0.021), 10 расхождений, hap.py ΔF1 −5.5e−5.
- **Problem:** самый сильный по мощности эксперимент проекта не отражён; при этом он единственный «chromosome-scale» и единственный с полным strict DEGRADED. Утверждение «no chromosome-scale run» — ложно на сегодня.
- **Correction:** добавить v20/chr20 в Table 1/2/3, Abstract, §3.1, §3.3 (Table 4), §4, §5, §6. Сохранить **оба** факта: «statistically classified DEGRADED under the pre-registered rule» и «ΔF1 = −0.000056 (−0.0056 п.п.), 18× меньше порога 0.001, 10 discordant loci из 56.2 M, net −1 TP / +7 FP / +1 FN». Не писать «no statistically significant degradation» и не писать «large degradation».

### CF-2 — CRITICAL — утрачены сырые источники ключевых HG005-чисел; в дереве лежит только инвалидированный прогон
- **Claim:** Table 3/§3.1 (HG005 hap.py 0.9521/0.9519, ΔF1 −0.0002), §3.3 M-1 (4,090/3,948/142, 3.5 %), §3.4 (13 true-SNP rescue / 129 FP avoidance), HG005 strata/enrichment, ERROR_ANALYSIS.
- **Evidence:** `CLAIM_EVIDENCE_MAP.csv` ссылается на 13 файлов, которых нет: `cache/hg005_stress_postfix/chr1_1000000_4000000.npz`, `experimental/stress_test/{HG005_POSTFIX_RESULTS.csv, HG005_POSTFIX_STRESS_TEST_REPORT.md, HG005_COORDINATE_INTEGRITY_REPORT.md, HG005_EXTRACTION_BUG_ROOT_CAUSE.md, HG005_HISTORICAL_AUDIT.md, HG005_ERROR_ANALYSIS.csv, HG005_RESCUE_ANALYSIS.csv}`, `experimental/performance/{CHTSLIB_INTEGRATION_REPORT.md, PERFORMANCE_OPTIMIZATION_REPORT.md}` (второй сохранился только в `.claude/worktrees/agent-a4ced…`), `bench_htslib_io/{END_TO_END_MULTIPROCESSING.md, PROFILING_ANALYSIS.md}`. В `experimental/stress_test/happy/ai_{pb_only,cascade}_C/happy.summary.csv` лежит **pre-fix** результат (TRUTH.TOTAL 28,435, F1 0.0085).
- **Problem:** числа HG005 внутренне согласованы (142/4,090 = 3.47 %; 3,948+142 = 4,090; 3,846/3,948 = 0.9742) и не противоречат ничему, но **в этом чекауте не верифицируемы**; более того, невалидный F1 0.0085 лежит рядом под теми же именами и может быть принят за результат.
- **Correction:** восстановить файлы из бэкапа или пересоздать (post-fix HG005 hap.py + region-scoped truth + cache); пометить/переместить pre-fix `happy/` в `_INVALIDATED_prefix`. До этого HG005-строки помечать «по записи; raw отсутствует».

### CF-3 — HIGH — «Native backend» и тесты: заявленное покрытие не существует
- **Claim:** §3.7 «12-category permanent test suite», «102/102 tests pass (90 + 12)», C1 «11 regression tests», Data availability: файлы `test_hg005_extraction_integrity.py`, `test_native_pileup_equivalence.py` существуют в worktree.
- **Evidence:** обоих файлов нет (`ls`); `LOAD_READS_OPTIMIZATION_REPORT.md` («Not restored»). Фактическое покрытие сейчас: 17 (reads) + 79 (5 suites) = **96 passed**; `test_providers.py` не собирается (нет `bimamba_variant_caller`).
- **Problem:** 12-категорийный тест на **counts**-backend (в т.ч. overlapping mates, orphans) и 11 coordinate-тестов утрачены; «102/102» — историческое измерение, не воспроизводимое. См. также CF-4.
- **Correction:** удалить/пометить «102/102» как «исторически, тесты утрачены»; переписать два теста заново; в тексте утверждать только текущее покрытие (96 passed + не собираемый `test_providers.py`). Не писать «all historical tests preserved».

### CF-4 — MEDIUM/HIGH — «bit-exact для counts» шире, чем текущее доказательство
- **Claim:** §2.7/§4.5: одна streaming `bam_mplp` проходка по всему span'у — «bit-exact replacement for per-window pysam loops».
- **Evidence:** `native/pileup_native.c:128–143` — один `bam_mplp_auto` по всему span'у с `bam_mplp_init_overlaps`. Для `reads_native` тот же приём **провалил** gate (30 клеток расходится на BAM с повторяющимися именами ридов: overlap-adjustment спаривает риды по имени среди живых в pileup), и финальный дизайн — независимый проход на каждое 64-nt окно (`LOAD_READS_OPTIMIZATION_REPORT.md`, A4). Для `pileup_native` эквивалентность подтверждена на реальных данных (0/109,813,120; 0/2,403,008 loci; chr20 gate 980,928 loci), но adversarial-fuzz на repeated names для **counts** утрачен вместе с тестом.
- **Problem:** тот же механизм расхождения теоретически применим к counts. Реальные данные его не проявили; «exact by construction» для counts не доказано.
- **Correction:** сузить формулировку до «bit-identical on all tested real regions (list)»; вернуть fuzz-тест (repeated names, overlapping mates) для counts; **и** обязательно сохранить в тексте историю провала первой reads-native версии (30 клеток) — сейчас в рукописи её нет (есть только провал F10 с ignore_orphans/overlaps).

### CF-5 — HIGH — устаревшая performance-картина (1.22×, «load_reads: none», 98.7 %)
- **Claim:** Abstract, §3.8, §4.6, §5, §6, Fig. 6: end-to-end 1.22×, `load_reads` без ускорения, PB+`load_reads` = 98.7 % wall.
- **Evidence:** `LOAD_READS_OPTIMIZATION_REPORT.md`: function-level 26.9× (HG005 30×, 400 kb) и 35.1× (HG002 15×); pipeline HG005 800 kb: 2.72× serial, native+4 процесса 7.11× (vs pure-Python serial); HG005 3 Mb: 309 s (1 proc) / 146 s (4 proc); PB = 287/309 s = 93 % (chr20: 95.9 % stage compute).
- **Problem:** 1.22× относится к конфигурации «native counts only» и верно как историческое измерение, но не описывает текущий production. Bottleneck теперь PB (93–96 %), а не `load_reads`; «Amdahl 3.1× для load_reads» — это ограничение уже реализованного шага (наблюдено 2.72× при потолке 3.06×).
- **Correction:** оставить 1.22× как «stage-1 (counts only)»; добавить load_reads-этап отдельной строкой; не смешивать function/extraction/pipeline/projection (см. §6).

### CF-6 — MEDIUM — воспроизводимость: hard-coded пути и tracked-артефакты
- **Evidence:** захардкожены абсолютные пути `/home/mark/…` и **`.claude/worktrees/opt-chr20`** в: `experimental/chr20_validation/{old_path_check.py, freeze_manifest.py, run_chr20_pipeline.py, run_happy_chr20.sh, run_extract_chr20.sh}`, `experimental/load_reads_opt/{common.py, run_e2e_hg005.sh}`, `bench_v20_chr20.py`, `test_native_reads_equivalence.py:221,228`, а также `experimental/{head_to_head/analyze_cascade.py, unified_happy/build_vcfs.py}` (`/tmp/head_to_head_work`). 79 `.pyc` в `__pycache__/` отслежены git'ом (4 изменены); `.gitignore` = `*.npz`, `data/`; `native/*.so` не игнорируется и не отслеживается; вся нативная работа, chr20, `research/` — untracked.
- **Correction:** см. §11 (MUST FIX).

### CF-7 — MEDIUM — рукописное «Data availability» устарело
- Указывает worktree под `/tmp/claude-…/hg005-fix-worktree` (стёрт при reboot 2026-09-20), файлы `experimental/{performance,stress_test}/` (первой нет), два несуществующих теста. Память проекта тоже говорит о `/tmp` — уже неверно; актуально: рабочее дерево `main` (modified/untracked) + `.claude/worktrees/opt-chr20` (ветка `fix/hg005-extraction-integrity`, base `a3d5761`).

### CF-8 — MEDIUM — глобальный «CONFIRMED»: правило неопределимо при текущих данных
- §2.8: «Overall CONFIRMED only if every powered cell and pooled result is PRESERVED». Пять pooled результатов = IMPROVED (v15, v17, v18, v19, v14-pooled), v14 pre-registered = NEGATIVE, v20 = DEGRADED. Правило как написано невыполнимо, при этом рукопись нигде не объявляет итоговый глобальный verdict — читатель может достроить его сам. См. §8 (предлагаемая формулировка).

---

## 3. Numerical consistency audit

| Quantity | Manuscript value | Current evidence | Status |
|---|---|---|---|
| AI PB-only+C hap.py P/R/F1 | 0.9785 / 0.9319 / 0.9547 | `unified_happy/.../ai_pb_only_C`: TP 15,748, FP 346, F1 0.954656 | CONFIRMED |
| AI Cascade+C hap.py P/R/F1 | 0.9783 / 0.9319 / 0.9546 | FP 349, F1 0.954569 | CONFIRMED |
| DeepVariant P/R/F1 | 0.9340 / 0.9654 / 0.9494 | H2H CSV 0.933997/0.965381/0.94943 | CONFIRMED |
| Clair3 P/R/F1 | 0.8965 / 0.9638 / 0.9290 | 0.896516/0.963842/0.928961 | CONFIRMED |
| GATK HC P/R/F1 | 0.9877 / 0.9493 / 0.9682 | 0.987748/0.949343/0.968165 | CONFIRMED |
| Старая пара 0.9893/0.9228 | явно отвергнута (§3.6) | 0.9894/0.9594 — internal evaluator, 16,597 denom.; 0.9228 нигде нет | CONFIRMED (как опровержение). Вторичные док-ты: не найдено использования как hap.py |
| Cascade-specific true-SNP events | 9 (Table 4; 7 v14+v17, +1 HG004, +1 v15) | 5+1+2+1 = 9 **без chr20**; chr20 добавляет 10-е (`chr20:46,110,517`) | AMBIGUOUS — определить, включён ли chr20 |
| «7» в тексте | Table 4 «7 in HG002 discovery/validation as counted in earlier reports» | происхождение объяснено | CONFIRMED (объяснено); прочих неспецифицированных 7/8 в рукописи не найдено |
| `Lost SNP` covariates | VAF 0.110–0.130, BQ 36.7–38.8, depth 27–107 | chr20: VAF 0.125, **BQ 39.2**, depth 28, `alldifficult` но **не** segdup | OUTDATED (диапазон BQ, признак segdup) |
| v14 verdict | NEGATIVE; pooled IMPROVED | `History/14_DEVLOG.md §21.1` | CONFIRMED |
| v15 «WEAK POSITIVE/NEGATIVE» | относится к safety layer | `results/bench_v15/summary.md:197,208` | CONFIRMED |
| Safety layer: fixed 5 / broke 69, captured 74/80 (92.5 %) | так | claim map (History/15) | CONFIRMED (по записи) |
| chr20 truth SNP in BED | (нет в рукописи) | 71,387 | CONFIRMED (json), отсутствует в рукописи |
| chr20 dropped | (нет) | 1,030 = 1.44 % (0.014428) | CONFIRMED |
| chr20 «retained/evaluable» | — | **три разных знаменателя**: retained 70,357; в SNP-frame 70,324 (−33 перекрыты indel-label); hap.py TRUTH.TOTAL 71,333 | CORRECT BUT NEEDS QUALIFICATION — «70,324» не равно 71,387−1,030 |
| chr20 locus F1 PB / Cascade | — | 0.97809 / 0.97803 | CONFIRMED |
| chr20 hap.py F1 PB / Cascade | — | 0.962851 / 0.962796 | CONFIRMED |
| chr20 ΔF1 | — | −5.64e-5, CI [−1.049e-4, −1.41e-5]; block [−1.08e-4, −7.5e-6] | CONFIRMED (CI — locus-level; для hap.py CI нет) |
| chr20 disagreements | — | 10 из 56,242,693 scored | CONFIRMED |
| chr20 routed | — | 0.1373 % (77,248) | CONFIRMED |
| chr20 rescuable / captured | — | 2,826; 2,817 (99.68 %); true-SNP 1,328 (1,327 захвачено), FP-avoid 1,498 (1,490) | CONFIRMED |
| Диапазон routed | 0.003–0.22 % (cell), 0.04–0.11 % (pooled) | chr20 0.137 %; HG002 chr21 0.155 % | OUTDATED |
| HG005 M-1 | 4,090 / 3,948 / 142 = 3.5 % | арифметика верна; raw утрачен | CONFIRMED по записи† |
| HG002 chr21 M-1 | 16,597 vs 16,898 → 301 (1.8 %) | claim map | CONFIRMED по записи |
| «~35 %» | опровергнуто | 142/4,090 = 3.47 % | CONFIRMED; остаётся в `LIMITATIONS_AND_OPEN_QUESTIONS.md` как пометка «wrong» — приемлемо |
| HG005 rescue 142 = 13 + 129 | так | raw утрачен | CONFIRMED по записи† |
| Rescue recall HG005 0.986 | в таблице §3.4 с оговоркой | `RESULTS.md:65`, `RESEARCH_SUMMARY_1PAGE.md:40` — без оговорки | OVERSTATED во вторичных док-тах |
| load_counts speedup | 21.0× (HG005 3 Mb), 65.1× (HG002 12 Mb), 29.9× (0.2 Mb) | `PERFORMANCE_RESULTS.csv`; 65.1× — переиспользованное измерение | CONFIRMED (stage-level, 65.1× оговорено) |
| load_reads function-level | «none» | 26.9× (HG005 30×), 35.1× (HG002 15×) | OUTDATED |
| E2E HG005 3 Mb | 1,392 → 1,140 s = 1.22× | 309 s (1 proc), 146 s (4 proc) новый код | OUTDATED (как описание текущего кода) |
| Pipeline speedups | не упомянуты | 2.72× serial; 7.11× native+4 proc (vs pure-Py serial; вклад только процессов ≈ 99.9/38.2 = 2.6×) | отсутствует |
| Threads | htslib BGZF ≈ 0×; multiproc 3.1× | Python threads на load_reads 0.73× (`LOAD_READS…`) | отсутствует в рукописи |
| Async | не упомянут (grep: 0) | 1.21× (Python reads), 1.05× (native) — только pipeline-level | OK: не переоценён; при добавлении не называть load_reads-ускорением |
| PB share | 98.7 % (load_reads+PB) | PB 93 % (287/309 s HG005), 95.9 % (chr20 stage compute) | OUTDATED |
| Amdahl | 1.33× (counts), 3.1× (load_reads), 1.44× (PB) | 3.06× потолок vs 2.72× наблюдено; PB-free ≈ 14× | CORRECT BUT NEEDS QUALIFICATION (исторические bounds vs новые) |
| Equivalence load_reads | — | 0/6.7·10⁸ cells; 8 fuzz BAM; 5 real chunks; 253,707 + 320,000 = 573,707 ≈ 574k BED-запросов; 0 mismatches | CONFIRMED (report); тест 17 passed воспроизведён |
| Tests | 102/102 | сейчас 96 passed (79+17); 2 файла утрачены | UNSUPPORTED (102) |
| Whole-genome estimates | помечены extrapolation | chr20: linear extrapolation «~1.5 days» помечена | CONFIRMED (как extrapolation) |

† raw-источник отсутствует в чекауте (CF-2).

---

## 4. Scientific claim audit

Категории: CONFIRMED · OUTDATED · UNSUPPORTED · OVERSTATED · CONTRADICTED · AMBIGUOUS · CORRECT BUT NEEDS QUALIFICATION. «†» — raw утрачен.

| ID | Claim в manuscript | Где | Источник доказательства | Подтверждено? | Число актуально? | Scope корректен? | Действие |
|---|---|---|---|---|---|---|---|
| C01 | Константы 7.0 / 10.5 / 5.4119 подобраны на chr21:32–40 Mb validation, заморожены | §2.4, Abstract | `cascade.py` sha b0ee9f4b… (пересчитан); FROZEN_MANIFEST | CONFIRMED | да | да | — |
| C02 | `cascade.py`, genotype-модуль byte-identical до/после | §2.4 | sha256 пересчитан = manifest chr20 | CONFIRMED | да | да | — |
| C03 | Конфликт devlog 11 (F-beta) vs worktree (0.1 % quantile) раскрыт | §2.4 | worktree `agent-aaff…` существует | CONFIRMED | да | да | — |
| C04 | U-H2 (chr21:32–44) содержит 8 Мб подбора, не blind holdout | Abstract, §2.4, §3.6, §5 | claim map; §2.4 | CONFIRMED | да | да | сохранить; добавить ту же оговорку к chr20 (тот же образец) |
| C05 | Внешнее сравнение «does not support a ranking» | Abstract, §3.6, §4.7 | H2H CSV | CONFIRMED | да | да | — |
| C06 | Pre-registered rule: PRESERVED/IMPROVED/DEGRADED/UNDERPOWERED, defect клаузы раскрыт | §2.8 | `bench_v14_crosschrom.classify`; devlog 13 | CONFIRMED | да | да | — |
| C07 | «Overall CONFIRMED only if every powered cell is PRESERVED» | §2.8 | Table 2 + v20 | AMBIGUOUS | — | — | CF-8; определить итоговую формулировку (§8) |
| C08 | v13: PB 0.98162, C 0.98172, PRESERVED (12/13) | Table 2 | claim map; TABLES_FINAL | CONFIRMED | да | да | — |
| C09 | v14: +0.000581 [+0.000075,+0.001089], **NEGATIVE**, pooled IMPROVED | Table 2, §3.1 | `History/14 §21.1` | CONFIRMED | да | да | — |
| C10 | v14 NEGATIVE вызван одной ячейкой chr20_neutral (2 FP) | §3.1 | `results/bench_v14/summary.md:5` | CONFIRMED | да | ⚠ имя | уточнить: это `chr20:33–34 Mb`, не v20 |
| C11 | v15 IMPROVED cascade-vs-PB; «WEAK POSITIVE/NEGATIVE» — safety layer | Table 2, §3.1 | `bench_v15/summary.md` | CONFIRMED | да | да | — |
| C12 | v16 PRESERVED; v17/v18/v19 IMPROVED (числа) | Table 2 | claim map / TABLES_FINAL | CONFIRMED | да | да | — |
| C13 | «IMPROVED» = артефакт PB FP в segdup, «not a better caller» | §3.1, §4.1 | 26/41, 36/38, 83/87 | CONFIRMED (как интерпретация) | да | да | — |
| C14 | «Seven regional experiments and two hap.py evaluations» | §1, Abstract | v20 существует | OUTDATED | нет | — | 8 internal + 3 hap.py |
| C15 | «|ΔF1| ≤ 0.0027 in seven experiments» | Abstract | chr20 −5.6e-5 внутри диапазона | CORRECT BUT NEEDS QUALIFICATION | числово да | — | добавить chr20 + DEGRADED |
| C16 | routed 0.003–0.22 % / 0.04–0.11 % | Abstract, §3.1 | chr20 0.137 %, chr21 0.155 % | OUTDATED | нет | — | обновить диапазоны |
| C17 | HG002 hap.py ΔF1 −0.0001; 3 из 10,975,654 loci | §3.1 | happy.summary (пересчитан) | CONFIRMED | да | да | — |
| C18 | HG005 hap.py ΔF1 −0.0002; 2 loci, оба cascade FP | §3.1 | по записи† | CONFIRMED† | да | да (1 регион) | CF-2 |
| C19 | «It is not lossless: nine cell-level events» | Abstract | 5+1+2+1 | CORRECT BUT NEEDS QUALIFICATION | 9 без chr20; 10 с ним | да | зафиксировать определение (C-1) |
| C20 | «Cascade never loses a true SNP» — false | §3.3 | devlog 14 | CONFIRMED | да | да | — Формулировки «never» в рукописи не найдено |
| C21 | True-SNP losses редки, концентрируются в difficult/low-VAF | §3.3, §6 | Table 4 + chr20 (1 из 10 не segdup) | CORRECT BUT NEEDS QUALIFICATION | — | segdup-часть сигнатуры не универсальна | «low-VAF, high-BQ, unrouted; чаще segdup/alldifficult» |
| C22 | Cascade-introduced FP: 15 events (low BQ) | Table 4 | v14 8 + v16 2 + v19 3 + HG005 2 | CONFIRMED (без v20) | 15 → 23 с chr20 | да | добавить 8 из chr20 |
| C23 | PB FP avoided: 153 | Table 4 | 26+36+7+83+1 | CONFIRMED (без v20) | 154 с chr20 | да | добавить 1 |
| C24 | «Дом. класс IMPROVED — segdup PB FP; cascade нечего выигрывать» | §4.1 | chr20: только 1 такое событие; 8/10 = router_fp | OUTDATED | — | — | добавить chr20-контрпример; знак ΔF1 зависит от состава участка |
| C25 | Failure mode низкий BQ характеризует FP-класс, не lost-SNP | §3.3 | chr20 signature (BQ 11–22; lost BQ 39.2) | CONFIRMED | да | да | — |
| C26 | Router band [1.588, 12.412]; router ≠ error detector | §2.4 | arithmetic; chr20 margins 5.50–6.65 / 6.33 | CONFIRMED | да | да | — |
| C27 | Router rescue HG005: 142 = 13 + 129; 0.986 не «true-SNP recall» | §3.4 | по записи† | CONFIRMED† | да | ⚠ | сохранить оговорку |
| C28 | «band captured 98.6–99.7 % of PB-rescuable loci» | §4.2 | HG005 91 % FP-avoid; chr20 true-SNP 0.9992 / FP 0.9947 | CORRECT BUT NEEDS QUALIFICATION | да | смешивает виды | развести на две популяции |
| C29 | Composition 13/129 не переносится на другие наборы | §3.4 | chr20 1,328/1,498 | CONFIRMED | да | да | в рукописи не переносится (OK); chr20 — отдельная строка |
| C30 | Метод C: GT acc 0.9888, F1 0.9538 (16,124-call scope) | Table 5a | claim map | CONFIRMED | да | да | — |
| C31 | GT accuracy 0.646 → 0.989 «and F1 to 0.955» | Abstract | 0.9888 — scope 16,124; 0.9547 — scope 16,094 (сама рукопись: «не объединять») | CORRECT BUT NEEDS QUALIFICATION | да | смешение scope | привести числа одного scope |
| C32 | Genotype layer не меняет allele call set | §3.5 | symmetric diff 0 | CONFIRMED | да | да | — |
| C33 | Genotype effect = evaluation-representation, not detection | §4.4 | root cause 5,630/5,777 | CONFIRMED | да | да | — |
| C34 | F1 0.6238 = writer artefact | §3.5 | unified_happy summary | CONFIRMED | да | да | — |
| C35 | External callers: region chr21:32–44, thresholds fitted 32–40 | §2.9, §3.6 | §2.4 | CONFIRMED | да | да | не называть blind holdout — не называется |
| C36 | «AI arms between DeepVariant and GATK in F1» | Abstract | 0.9547 vs 0.9494 / 0.9682 | CONFIRMED | да | да (с оговоркой) | — |
| C37 | M-1 HG005 4,090/3,948/142 = 3.5 % | §3.3 | †; арифметика | CONFIRMED† | да | да | CF-2 |
| C38 | M-1 HG002 chr21: 301/16,898 = 1.8 % | §3.3 | claim map | CONFIRMED | да | да | — |
| C39 | M-1 chr20: 1,030/71,387 = 1.44 % | (нет) | json | CONFIRMED (отсутствует) | — | — | добавить |
| C40 | Каждый F1 имеет явный evaluator и denominator | §2.8, таблицы | таблицы имеют evaluator; denominators M-1 частично | CORRECT BUT NEEDS QUALIFICATION | — | chr20: 3 знаменателя | явная таблица denominators |
| C41 | «Every AI recall is depressed by M-1; Table 3 conservative for AI recall» | §3.3 | HG005/chr21/chr20 (1.4–3.5 %) | CONFIRMED | да | да | — |
| C42 | Native: first attempt failed (≤282,335/1,684,800; ignore_orphans/overlaps) | §3.7, F10 | claim map | CONFIRMED | да | да | — |
| C43 | Native reads: first design (one pass over run of windows) → 30 mismatching cells; final = per-64-nt window | (нет) | LOAD_READS report A4 | UNSUPPORTED (в рукописи отсутствует) | — | — | добавить, не прятать |
| C44 | Native counts «bit-exact replacement» | §2.7, §4.5 | CF-4 | OVERSTATED | — | да, для перечисленных данных | сузить до «identical on tested regions» |
| C45 | 0/109,813,120 cells (12 Mb); 0/2,403,008 (HG005) | §3.7 | по записи (`experimental/performance/…` отсутствует) | CONFIRMED† | да | да | CF-2 |
| C46 | «12-category permanent test suite pass» | §3.7 | файл утрачен | UNSUPPORTED | — | — | CF-3 |
| C47 | «102 / 102 tests pass» | §3.7 | сейчас 96 | UNSUPPORTED | нет | — | CF-3 |
| C48 | C1: «11 regression tests» для coordinate fix | §3.9 C1, disclosure | файл утрачен | UNSUPPORTED | — | — | CF-3 |
| C49 | Coordinate bug: не инвалидировал HG002/3/4; v16–v19 не перезапускались с исправленным extractor | §3.9 C1 | claim map «blast radius one script» (HG005_HISTORICAL_AUDIT утрачен) | CONFIRMED† + честно раскрыто | да | да | CF-2 |
| C50 | Старые fabricated positions/`chunk_start+arange` не использованы в текущих числах | §3.9 C1 | grep: v13–v19 — row-index; chr20 — реальные positions + 0 ref-mismatches, `positions` strictly increasing | CONFIRMED | да | да | — |
| C51 | load_counts 21.0× / 65.1× — extraction stage only | Abstract, §3.7 | PERFORMANCE_RESULTS.csv | CONFIRMED | да | да | 65.1× — переиспользованное |
| C52 | «Extraction speedup 65×, but pipeline 1.22×» | Abstract, §6 | — | OUTDATED | 1.22× исторически | — | добавить load_reads-этап |
| C53 | `load_reads` без ускорения (696→775 s) | §3.8 | LOAD_READS 26.9×/35.1× | OUTDATED | нет | — | обновить |
| C54 | PB + load_reads = 98.7 % | Abstract, §3.8 | PB 93–96 % | OUTDATED | нет | — | обновить |
| C55 | Amdahl 3.1× для load_reads | §3.8 | 3.06× потолок, 2.72× наблюдено | CORRECT BUT NEEDS QUALIFICATION | ≈ | да | указать наблюдённое |
| C56 | Withdrawn: «65× end-to-end», «17–80× faster than callers» | §3.8 | — | CONFIRMED (honest) | да | да | — |
| C57 | Routed fraction ≠ CPU savings; 1.65× измерено (HG004) | §3.8 | timing.json | CONFIRMED | да | да | chr20 подтверждает: PB считается для всех loci |
| C58 | Threads: htslib BGZF ≈ 0×; Python threads 0.73× | §3.8 (BGZF only) | LOAD_READS | CORRECT BUT NEEDS QUALIFICATION | — | — | добавить 0.73× |
| C59 | Async — не заявлен как ускорение load_reads | (нет упоминаний) | LOAD_READS: 1.21×/1.05× pipeline | CONFIRMED (нет overclaim) | — | — | при добавлении — pipeline-level только |
| C60 | «Chromosome- or genome-scale run: none» | §5 | chr20 whole | CONTRADICTED | — | — | обновить: один хромосомный run (chr20, HG002, 15×) |
| C61 | «Largest single span 12 Mb» | §5 | chr20 64.4 Mb | CONTRADICTED | — | — | обновить |
| C62 | «All timings single-machine, mostly single-run» | §5 | chr20 6 workers, contention-inflated | CORRECT BUT NEEDS QUALIFICATION | — | — | добавить оговорку chr20 |
| C63 | Related-trio caveat HG003/HG004 | §3.2, §5 | — | CONFIRMED | да | да | — |
| C64 | HG005 «supports neither generalises nor population claim» | §3.2 | — | CONFIRMED | да | да | вторичные: `LIMITATIONS.md:25` «demonstrate … generalizes across trio members» — OVERSTATED |
| C65 | SNP-only, GRCh38, GIAB; indels/MNP/SV не исследованы | §2.1, §5 | — | CONFIRMED | да | да | — |
| C66 | Конституты глубина-специфичны (11.5/14.0 @30×) | §2.4, §5 | worktree devlog 10 | CONFIRMED (по записи) | да | да | chr20 также 15× |
| C67 | chr20: hap.py CI отсутствует, locus CI есть | (нет) | report §6 | CONFIRMED | — | — | явно разделить в тексте |
| C68 | chr20 частично не unseen (1 Mb в v14) | (нет) | sensitivity excl. 33–34 Mb: −5.70e-5 | CONFIRMED | да | да | добавить |
| C69 | chr20 — same sample, chromosome-level separation | (нет) | report §2 | CONFIRMED | да | да | добавить; не называть independent |
| C70 | F5: «No cutoff between 4.5 and 8.0 removed the DEGRADED chr20 cell» | §3.9 F5 | `History/14 …:641` — про **v14 chr20_neutral**; v20: sweep −30 %…+20 % (F1 0.97722–0.97808), потеря не устраняется | AMBIGUOUS | — | коллизия имён | указать «v14 chr20:33–34 Mb» vs «v20 whole chr20» |
| C71 | Вторичные: «Cascade preserves PB-only accuracy in every tested experiment» | `RESEARCH_SUMMARY_1PAGE.md:37` | v14 NEGATIVE, v20 DEGRADED | CONTRADICTED | — | — | пометить superseded |
| C72 | Вторичные: Rescue recall 0.986 без оговорки | `RESULTS.md:65`, `SUMMARY_1PAGE:40` | 13/129 | OVERSTATED | — | — | superseded-метка |
| C73 | Вторичные: `RESULTS.md` v14 «IMPROVED» | `RESULTS.md:16` | pre-registered NEGATIVE | OVERSTATED (уже записано D1 в TABLES_FINAL) | — | — | superseded-метка |
| C74 | Стоимость PB/binomial «≈240–490×» | §4.2 | 2.6–3.1e6 vs 6–10e3 loci/s → 260–520×; 0.12–0.15 ms | AMBIGUOUS | ≈ | — | выровнять с §3.8 |
| C75 | «Constants frozen before any later validation experiments» | §4.1 | v13 включает 3 региона внутри 32–40 Mb (раскрыто) | CORRECT BUT NEEDS QUALIFICATION | — | — | уже частично оговорено |
| C76 | Hardware/env versions, hap.py digest | Data availability | ANALYSIS_FREEZE | CONFIRMED | да | да | — |
| C77 | Data availability: worktree `/tmp/…`, 2 тестовых файла, `experimental/performance` | §Data availability | CF-7 | CONTRADICTED | — | — | обновить |
| C78 | chr20 hap.py conclusion == locus conclusion | (нет) | −5.5e-5 vs −5.6e-5 | CONFIRMED | да | да | добавить |

**Итог по таблице (78 claims):** CONFIRMED 48 (из них 8 «†/по записи»), CORRECT BUT NEEDS QUALIFICATION 10, OUTDATED 6, UNSUPPORTED 4, OVERSTATED 3, CONTRADICTED 4, AMBIGUOUS 3.

---

## 5. Internal contradiction audit

| # | Места | Противоречие |
|---|---|---|
| IC-1 | Abstract/§1/§3.1 («seven … two hap.py») ↔ существование v20 | Пересчёт: 8 internal + 3 hap.py. |
| IC-2 | §3.8 «`load_reads`: none / 98.7 %» и Fig. 6 ↔ `LOAD_READS_OPTIMIZATION_REPORT.md` (26.9×/35.1×, PB 93–96 %) | Рукопись описывает промежуточный, не текущий production. |
| IC-3 | §5 «no chromosome-scale run; largest span 12 Mb» ↔ chr20 (64.4 Mb) | Прямое противоречие. |
| IC-4 | §2.8 правило «CONFIRMED только если все PRESERVED» ↔ Table 2 (5×IMPROVED, 1×NEGATIVE) | Правило невыполнимо; итоговый verdict нигде не задан (CF-8). |
| IC-5 | Abstract «nine events (eight in consolidated catalogue + one v15)» ↔ Table 4 «7 + 1 + 1» ↔ Fig. 4 «eight … with covariates» | Согласовано пояснениями, но три способа счёта; chr20 сделает 10. |
| IC-6 | §2.7 «одна проходка по span'у» (counts) ↔ LOAD_READS «whole-run pass diverged (30 cells) → per-window» | Не противоречие данных, а несимметрия доказательств (CF-4). |
| IC-7 | §3.7 «102/102 …12-category suite» ↔ §Data availability (те же тесты «в worktree») ↔ LOAD_READS «Not restored» | Утверждение о существовании ↔ факт утраты. |
| IC-8 | Abstract: GT acc 0.989 и F1 0.955 в одной фразе ↔ §3.5.7 «0.9547 и 0.9538 не объединять» | Смешение scope. |
| IC-9 | §3.3 signature lost SNP «segdup» ↔ chr20 lost SNP (только alldifficult) | Сигнатура нуждается в расширении. |
| IC-10 | F5 «chr20 cell» ↔ v20 «chr20» | Коллизия имён (chr20:33–34 Mb vs whole chr20; первый включён во второй). |
| IC-11 | §3.2 «HG005 = only unrelated sample» и чистые HG005-числа ↔ raw утрачен | Внутренне ок, не проверяемо (CF-2). |
| IC-12 | Figure legends Fig. 4/6/7 ↔ новые данные | Fig. 6 (1,392→1,140 s) описывает промежуточную конфигурацию; Fig. 4 не содержит chr20. |

Между Methods/Results/Limitations внутренних числовых расхождений сверх перечисленных не найдено.

---

## 6. Performance audit

Четыре уровня, которые нельзя смешивать (проверено, что в LOAD_READS-отчёте они разделены; в рукописи load_reads-уровня нет):

| Уровень | Что измерено | Значение | Статус |
|---|---|---|---|
| **Function** (`load_reads`, только эта функция) | HG005 30×, 400 kb, native+bisect BED, serial | **26.9×** (79.95 → 2.98 s) | CONFIRMED (по отчёту; 3 повтора; фон ≈1–2 ядра) |
| Function | HG002 15×, 400 kb | **35.1×** (89.63 → 2.56 s; baseline n=2) | CONFIRMED (n=2 у baseline) |
| **Function** (`load_counts`, extraction stage) | HG005 3 Mb / HG002 12 Mb / 0.2 Mb | 21.0× / 65.1× / 29.9× | CONFIRMED; 65.1× — переиспользованное, single run |
| Python threads (load_reads) | 4 threads | **0.73×** | CONFIRMED (отчёт); в рукописи отсутствует |
| Multiprocessing (Python) | 4 / 8 workers, load_reads | 2.86× / 4.39× | CONFIRMED |
| **Pipeline E2E** (HG005 800 kb, serial) | native reads + counts | **2.72×** vs pure-Python serial | CONFIRMED |
| Pipeline E2E, native + 4 процесса | 800 kb | **7.11×** vs pure-Python serial (**вклад процессов относительно native serial ≈ 2.6×**) | CORRECT BUT NEEDS QUALIFICATION — не называть «multiprocessing 7.1×» без базы |
| Pipeline HG005 3 Mb | новый код | 309 s (1 proc), 146 s (4 proc) → 2.1× от процессов | CONFIRMED |
| Историч. pipeline stage-1 | native counts only | 1,392 → 1,140 s = 1.22× | CONFIRMED как историческое, OUTDATED как описание кода |
| Async (pipeline-level) | asyncio prefetch | 1.21× (Python reads), 1.05× (native) | не измерено на `load_reads` отдельно; не заявлять как ускорение load_reads |
| Bottleneck | PB | 287/309 s = 93 % (HG005); 95.9 % stage compute (chr20, при 6 workers → contention-inflated) | CONFIRMED; рукопись 98.7 % устарела |
| Amdahl | eliminating `load_reads` | потолок 3.06× (p = 183/272 = 0.67), наблюдено 2.72× (89 % потолка); PB-free ≈ 14× | CONFIRMED (отчёт); историч. 1.33/3.1/1.44 в рукописи — маркировать |
| Проекции | 171–396× caller-stage; genome ≈ 6.1 ч / 399 ч / «~1.5 суток» | линейные экстраполяции | не результаты (в рукописи помечены) |
| Equivalence | 0/6.7·10⁸ cells; 8 fuzz BAM; 5 real chunks; ≈574k BED-запросов | 17 passed воспроизведено | CONFIRMED |

Замечания: (i) время HG005 3 Mb 1,140 s vs новый PB 287 s не сопоставимы напрямую (другая нагрузка машины; отчёт это признаёт); (ii) RSS и «ядра» — с фоновой нагрузкой; (iii) chr20 PB-время 15,036 s — сумма по 6 workers на 4 физ. ядрах, inflated.

---

## 7. Validation audit

| Эксперимент | Что доказано | Ограничения (подтверждены) |
|---|---|---|
| HG002 v13 (13 регионов chr21) | PRESERVED 12/13; 0 lost | 3 региона внутри 32–40 Mb (подбор констант) |
| HG002 v14 (4 контига) | pooled IMPROVED; **pre-reg NEGATIVE**; 5 lost | chr20_neutral (2 FP) вызвал DEGRADED |
| HG002 v15 (held-out) | cascade vs PB IMPROVED; safety-layer rejected | verdict «WEAK POSITIVE» ≠ cascade |
| HG002 v16 / v17 | PRESERVED / IMPROVED (segdup) | 2 lost в v17 |
| HG003 v18, HG004 v19 | IMPROVED; 1 lost (HG004) | одна трио, один батч — слабая независимость |
| HG005 chr1:1–4 Mb (hap.py) | ΔF1 −0.0002; 2 loci | 1 регион; raw утрачен (CF-2); нет CI |
| HG002 chr21:32–44 (hap.py + внешние) | AI 0.9547/0.9546 между DV и GATK | перекрытие с подбором констант |
| **HG002 chr20 целиком** | 70,324 SNP, 56.2 M loci, DEGRADED, ΔF1 −5.6e-5, 10 discordant | тот же образец; chromosome-level, не sample-level; 1 Mb ранее в v14; только 15× |

Дополнительно: независимость v14+ регионов (`independence_audit.json` all_ok) не перепроверялась; определение M-1 (window-level BED) применяется единообразно, denominator drop: 3.5 % (HG005), 1.8 % (chr21), 1.44 % (chr20).

---

## 8. Statistical audit

- **Preregistered rule.** Реализация `classify` = как описано в §2.8; дефект клаузы `|ΔF1| ≥ 0.001` раскрыт. chr20: клауза 1 (CI целиком < 0) → **DEGRADED**. Верно также для block-bootstrap и sensitivity. Verdict не менять.
- **Magnitude отдельно от significance.** chr20: ΔF1 = −5.64·10⁻⁵ = −0.0056 п.п.; 18× внутри порога 0.001; net −1 TP / +7 FP / +1 FN на 70,324 SNP; 10 discordant loci; McNemar exact p = 0.021 (9:1). Допустимая формулировка: «*classified DEGRADED by the pre-registered rule; the effect is extremely small*». Недопустимо: «no significant degradation», «PRESERVED», «substantial degradation».
- **Мощность.** UNDERPOWERED-ячеек нет (CI half-width chr20 ≈ 5·10⁻⁵). Замечание отчёта: p = 0.021 опирается на 10 событий.
- **Bootstrap.** Locus-bootstrap считает соседние loci независимыми; block-CI шире (2.3× в v19; в chr20 block-CI верхняя граница −7.5e-6 — всё ещё < 0). hap.py — без CI (в обоих местах честно указано).
- **Множественные сравнения:** нет коррекции (раскрыто §5); 8 experiments + под-ячейки → один DEGRADED из большого числа проверок ожидаем; но правило pre-registered, поэтому verdict не отменяется.
- **Слова.** Проверены «significant/insignificant/negligible/meaningful/preserves/improves/degrades/confirmed» в манускрипте: «significant» встречается только в контексте прежних статистик; «preserves» в Abstract/Conclusions — **нуждается в квалификации после chr20**; слова «negligible» и «CONFIRMED» как итоговых оценок в рукописи нет (CONFIRMED — только определение §2.8). `IMPROVED` корректно ограничен интерпретацией.
- **Предлагаемая формулировка глобального итога (вместо «CONFIRMED»; данные не менять):**
  > «Under the pre-registered per-cell rule, cascade-vs-PB results were PRESERVED or IMPROVED in the large majority of evaluated cells/experiments; one pre-registered cell (v14, chr20:33–34 Mb, full depth) and the whole-chromosome chr20 15× run (v20) were classified DEGRADED. The v20 effect is very small (ΔF1 = −5.6×10⁻⁵, 10 discordant loci of 56.2 M) but statistically detectable with the available power. The broad statement ‘accuracy is preserved’ is therefore supported only as ‘preserved to within small margins in the evaluated GIAB/GRCh38 SNP-only setting’, not as an unconditional confirmation.»

---

## 9. Scope audit

| Что доказано | Что НЕ доказано |
|---|---|
| Внутри GIAB v4.2.1 high-confidence regions, GRCh38, **SNP-only** | Indels/MNP/SV не оценивались |
| Illumina 2×250 novoalign (HG002–4); HG005 300× HiSeq novoalign ≤250 bp, downsampled | Другие платформы, aligners (bwa-mem), read lengths |
| HG002/HG003/HG004 (одна трио, один батч) + HG005 (1 регион, 3 Mb) | Популяционное разнообразие, «generalises across samples» |
| Пороги заморожены на chr21:32–40 Mb @≈15×; не переподбирались | Оптимальность на других глубинах (на 30× оптимум 11.5/14.0) |
| chr20 целиком, HG002, 15× | Native-depth/30× chr20; другие хромосомы; геном целиком |
| U-H2 сравнение с DV/Clair3/GATK на одном регионе | Ранжирование callers; независимый blind-holdout |
| Exact equivalence native на перечисленных реальных данных | Общая «bit-exact» эквивалентность для произвольных BAM (counts) |

Запрещённых формулировок вроде «works on genomes in general» в рукописи не найдено; но вторичные документы содержат более сильные (`LIMITATIONS.md:25`, `RESEARCH_SUMMARY_1PAGE.md:37`, C64/C71).

---

## 10. Reproducibility audit

| Пункт | Находка |
|---|---|
| Hard-coded абсолютные пути | Да: `/home/mark/Documents/Projects/AI_DNA_ANALYZER/…` в 12+ файлах; **`.claude/worktrees/opt-chr20`** в `experimental/chr20_validation/{old_path_check.py, freeze_manifest.py, run_chr20_pipeline.py, run_happy_chr20.sh}`, `experimental/load_reads_opt/{common.py, run_e2e_hg005.sh}`; `/tmp/head_to_head_work` в `experimental/head_to_head/analyze_cascade.py`, `overlap_analysis_DISCARDED.py` → **MUST FIX** |
| Ссылки на удалённый `/tmp` worktree | `research/REPRODUCIBILITY.md:9–10`, рукопись «Data availability», `RESEARCH_TIMELINE.md` («Not resolved») |
| Отсутствующие файлы | 13 источников из `CLAIM_EVIDENCE_MAP.csv` (CF-2); `experimental/performance/` целиком; `test_hg005_extraction_integrity.py`, `test_native_pileup_equivalence.py`; `cache/hg005_stress_postfix/` |
| Инвалидированный прогон под теми же именами | `experimental/stress_test/happy/ai_*_C` (F1 0.0085) |
| Тесты | 17 + 79 passed воспроизведены; `test_providers.py` не собирается (модуль `bimamba_variant_caller` отсутствует, pre-existing) |
| Сгенерированные `.so` | `native/*.so` не отслеживаются и не игнорируются; сборка `native/build.sh` (нужны Python headers + pysam htslib) |
| Кэши в git | 79 `.pyc` отслежены (`__pycache__/`), 4 изменены; `.gitignore` не содержит `__pycache__/` |
| Фиксация | Всё нативное, chr20, `research/`, `results/bench_v20`, HG005-инфраструктура — untracked/modified (нет коммитов, по вашему правилу); `git` HEAD = `a3d5761` |
| Определённость окружения | Плюс: версии Python/numpy/scipy/pysam/hap.py digest зафиксированы; SHA-256 скриптов анализа зафиксированы **до** открытия результатов (`ANALYSIS_FREEZE.json`). Минус: нет lockfile |
| Детерминизм | chr20 re-run воспроизвёл потерянный запуск (56,266,816 loci, 70,324 SNP, 69,291/69,297 calls, 2,826 rescuable) — сильное подтверждение |
| Gate оптимизированного пути | Два полных chunk'а chr20 (980,928 loci) байт-в-байт равны пути без native |

---

## 11. Required fixes before manuscript freeze

### MUST FIX
1. **Включить v20/chr20 в рукопись** (Table 1/2/3, Abstract, §3.1, §3.3, §3.4, §4, §5, §6): DEGRADED + magnitude −0.0056 п.п. одновременно; 8/1/1 taxonomy; router 0.137 %; rescue 1,328/1,498; M-1 1.44 % с тремя знаменателями; same-sample/1 Mb overlap оговорки (CF-1, C14–C16, C39, C60–C61, C67–C70).
2. **Обновить performance-секции** (load_reads 26.9×/35.1×; E2E 309/146 s, 2.72×/7.11×; PB 93–96 %; Amdahl 3.06 vs 2.72; threads 0.73×; async — только pipeline-level 1.21×/1.05×); 1.22× оставить как «stage-1» (CF-5, C52–C55, C58).
3. **Заменить глобальный «CONFIRMED»** (§2.8) формулировкой из §8 и не оставлять невыполнимое правило без пояснения (CF-8, C07).
4. **Восстановить/пересоздать утраченные источники HG005** (post-fix hap.py, region-scoped truth, cache, error/rescue CSV, `experimental/performance` отчёты) и убрать риск путаницы с pre-fix `happy/` (CF-2, C18, C27, C37, C45, C49).
5. **Убрать утверждения о несуществующих тестах** («102/102», «12-category suite», «11 regression tests») либо переписать тесты; в тексте — только текущее покрытие (CF-3, C46–C48). Обновить Data availability (CF-7, C77).
6. **Заменить hard-coded пути в `experimental/` (`opt-chr20` и др.)** на относительные/конфигурируемые (CF-6).
7. **Определить счёт true-SNP событий** (9 без chr20 / 10 с ним) и расширить сигнатуру (BQ до 39.2; segdup **или** alldifficult) (C19, C21, IC-5, IC-9).
8. **Историю первой reads-native версии (30 mismatching cells)** включить в текст; сузить «bit-exact» для counts до «identical on tested regions» (CF-4, C43–C44).

### SHOULD FIX
- Развести true-SNP rescue и FP-avoidance в §4.2 (C28); привести Abstract-фразу про GT-accuracy/F1 к одному scope (C31).
- Явная таблица denominators (BED / retained / scored-frame / hap.py TRUTH.TOTAL) для HG005, chr21, chr20 (C40).
- Указать, что chr20 CI — locus/block-bootstrap, hap.py без CI (C67); переименовать v14 `chr20_neutral` в тексте, чтобы не путать с v20 (C10, C70).
- Пометить как superseded: `RESULTS.md:16,65`, `RESEARCH_SUMMARY_1PAGE.md:37–41`, `LIMITATIONS.md:25`, `SCIENTIFIC_CLAIMS.md` (~65×) (C71–C73, C64).
- Обновить Fig. 4/6/7 и легенды (chr20; новая bottleneck-картина).
- Вернуть fuzz-тест для counts-backend (repeated names, overlapping mates); убрать tracked `__pycache__`, добавить `.gitignore` для `__pycache__/`, `native/*.so`; lockfile.
- Выровнять стоимость PB/binomial «240–490×» с §3.8 (C74).

### OPTIONAL
- Указать, что chr20 PB-время inflated contention'ом; n=2 у HG002 baseline.
- Добавить sensitivity (excl. v14 1 Mb) в основной текст.
- Ссылочные `[REFERENCE NEEDED]` и `[DATA VERIFICATION NEEDED]` закрыть при финальной вычитке.

---

## 12. Final recommendation

**NOT READY**

Обоснование (только факты): научные числа, которые можно перепроверить в чекауте (chr20, hap.py AI и внешние callers, v14/v15 verdict'ы, хэши замороженных файлов, 96 тестов), **подтверждаются**, отрицательные результаты не смягчены. Но рукопись (а) не содержит chr20 и оптимизации `load_reads`, из-за чего как минимум 6 утверждений устарели и 2 (chromosome-scale, span 12 Mb) прямо противоречат существующим данным, (б) ссылается на 13 источников и 2 теста, которых физически нет, включая сырьё HG005-чисел, (в) не имеет определённого глобального verdict'а при невыполнимом правиле CONFIRMED. Это не «минорные правки»: нужна интеграция двух экспериментов и восстановление provenance. Новых экспериментов **не требуется** (кроме, возможно, пересоздания утраченного HG005 post-fix прогона, если файлы не найдутся в бэкапе). После пунктов MUST FIX статус меняется на **READY AFTER MINOR CORRECTIONS**.
