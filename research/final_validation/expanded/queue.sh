#!/usr/bin/env bash
cd /home/kramkost/Projects/Varion_DNA/AI_DNA_ANALYZER/research/final_validation/expanded
for sc in "HG007 chr17" "HG007 chr19" "HG001 chr17" "HG001 chr19"; do ./matrix_expanded.sh $sc; done
for sc in "HG007 chr19" "HG001 chr19" "HG007 chr17" "HG001 chr17"; do ./matrix_expanded.sh $sc dv; done
echo QUEUE_DONE
