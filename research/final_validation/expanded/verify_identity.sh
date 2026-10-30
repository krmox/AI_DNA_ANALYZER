#!/usr/bin/env bash
# Frozen-R2 identity check (read-only). Prints PASS/FAIL lines; exit 1 on any mismatch.
R=/home/kramkost/Projects/Varion_DNA/AI_DNA_ANALYZER; rc=0
chk(){ [ "$2" = "$3" ] && echo "OK   $1" || { echo "FAIL $1 got=$2 want=$3"; rc=1; }; }
chk HEAD "$(git -C $R rev-parse --short HEAD)" 3cdfe5a
n=$(cd $R/native_r2 && sha256sum -c /mnt/hdd/AI_DNA_ANALYZER_r2/gates/source_sha256.txt 2>&1 | grep -c ": OK"); chk sources_ok_of_21 "$n" 21
chk r2_binary "$(sha256sum /mnt/hdd/AI_DNA_ANALYZER_r2/build/dnav2 | cut -d' ' -f1)" 6c675174854917ca3bbbfc1192c9205ac0043a696bd37d6ff2b2c08c3a3304be
chk r2_model "$(sha256sum $R/research/stage4/frozen_R2.model | cut -d' ' -f1)" 649215d52e57cf00d9f69b8631d7c0d9a1848db2550fc593b7d464c1b4cab6ac
chk tracked_modified "$(git -C $R status --porcelain --untracked-files=no | wc -l)" 0
r11=$(sha256sum /mnt/hdd/AI_DNA_ANALYZER_stage3/build/dnav2 | cut -d' ' -f1); echo "info r1.1_binary $r11"; echo "info r1.1_model $(sha256sum $R/research/stage3/frozen/R1.model | cut -d' ' -f1)"
exit $rc
