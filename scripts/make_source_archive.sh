#!/usr/bin/env bash
# Build a lightweight source archive: code, tests, models, build scripts, docs and the few files
# verify.sh reads -- WITHOUT research/, checkpoints*/, deprecated/checkpoints/, experimental/, results/, data/, datasets/
# and docs/history.  The full repository (with all evidence) stays on GitHub.
#
#   bash scripts/make_source_archive.sh [OUT.tar.gz]              committed tree (git archive HEAD)
#   bash scripts/make_source_archive.sh --worktree [OUT.tar.gz]   tracked + untracked-not-ignored files of
#                                                                 the working tree (to test before committing)
set -euo pipefail
cd "$(dirname "$0")/.."
MODE=head; [ "${1:-}" = "--worktree" ] && { MODE=worktree; shift; }
OUT=${1:-varion-src-$(git rev-parse --short HEAD).tar.gz}
EXCL='^(checkpoints[^/]*|deprecated/checkpoints|research|data|datasets|experimental|results|docs/history|third_party|[.]venv|[.]worktrees|NEW version)(/|$)'
KEEP='^(research/stage4/frozen_R2[.]model|research/final_independent_eval/run/ai_A_C_chr8_300x[.]vcf[.]gz([.]tbi)?|research/final_independent_eval/checksums[.]txt)$'
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
if [ $MODE = head ]; then
  git ls-tree -r -z --name-only HEAD > "$TMP/all"
else
  { git ls-files -z --cached --others --exclude-standard; } | sort -zu > "$TMP/all"
fi
# keep = not excluded, or explicitly kept; must exist in the working tree (worktree mode)
tr '\0' '\n' < "$TMP/all" | awk -v ex="$EXCL" -v keep="$KEEP" '($0 !~ ex) || ($0 ~ keep)' > "$TMP/list"
if [ $MODE = worktree ]; then while IFS= read -r f; do [ -e "$f" ] && echo "$f"; done < "$TMP/list" > "$TMP/list2"; mv "$TMP/list2" "$TMP/list"; fi
if [ $MODE = head ]; then
  tr '\n' '\0' < "$TMP/list" | xargs -0 git archive --format=tar.gz -o "$OUT" HEAD --
else
  tr '\n' '\0' < "$TMP/list" | tar --null -czf "$OUT" -T -
fi
echo "wrote $OUT ($(du -h "$OUT" | cut -f1), $(wc -l < "$TMP/list") files, mode: $MODE)"
