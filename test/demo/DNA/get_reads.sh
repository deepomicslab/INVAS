#!/usr/bin/env bash
set -euo pipefail

FA="chr12_concat.fa"
R1="inv_wgs_R1.fq"
R2="inv_wgs_R2.fq"

L=$(grep -v '^>' "$FA" | tr -d '\n' | wc -c)
echo "Region length (bp): $L"

DEPTH=60
RL=150
INS=350
SD=50

N=$(( (DEPTH * L + (2 * RL - 1)) / (2 * RL) ))
echo "Read pairs to simulate: $N (PE${RL}, insert ${INS}±${SD})"

# Run wgsim (no inline comments)
wgsim \
  -N "${N}" \
  -1 "${RL}" -2 "${RL}" \
  -d "${INS}" -s "${SD}" \
  -e 0 \
  -r 0 \
  -R 0 \
  "$FA" "$R1" "$R2"

echo "Done."
echo "Outputs:"
echo "  - $R1"
echo "  - $R2"
