#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GRAVMAG_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

RSPHERE_KM="1737.4"
INPUT_IN="${SCRIPT_DIR}/reiner_gamma_table4_direct.in"
OUTPUT_XYZ="${SCRIPT_DIR}/output/reiner_gamma_table4_direct_xyz.txt"
OUTPUT_BRTP="${SCRIPT_DIR}/output/reiner_gamma_table4_direct_brtp.txt"
mkdir -p -- "${SCRIPT_DIR}/output"

"${GRAVMAG_DIR}/build_gravmag_tools.sh"

echo "+ running point_source solver on ${INPUT_IN}"
"${GRAVMAG_DIR}/gravmag_sphere_point_source" \
  "${RSPHERE_KM}" \
  "${INPUT_IN}" \
  "${OUTPUT_XYZ}" \
  "1" "0" "0" "0"

echo "+ converting XYZ output to Br/Btheta/Bphi"
"${GRAVMAG_DIR}/gravmag_xyz_to_brtp" \
  "${OUTPUT_XYZ}" \
  "${OUTPUT_BRTP}"

echo "Input IN: ${INPUT_IN}"
echo "Output XYZ: ${OUTPUT_XYZ}"
echo "Output BRTP: ${OUTPUT_BRTP}"
