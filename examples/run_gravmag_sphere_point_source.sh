#!/usr/bin/env bash
set -euo pipefail

# Point-source solver runner (XYZ output)
# Usage:
#   ./run_gravmag_sphere_point_source.sh <R_sphere_km> <input_in_file> [xyz_output_file] [refine_factor] [source_nlat] [source_nlon] [source_nr]
#
# Equivalent executable call:
#   ./gravmag_sphere_point_source R input.in output_xyz.txt refine source_nlat source_nlon source_nr

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_SPHERE_KM="${1:-}"
INPUT_FILE="${2:-}"
XYZ_OUT="${3:-}"
REFINE_FACTOR="${4:-2}"
SOURCE_NLAT="${5:-0}"
SOURCE_NLON="${6:-0}"
SOURCE_NR="${7:-0}"

if [[ -z "${R_SPHERE_KM}" || -z "${INPUT_FILE}" ]]; then
  echo "Usage: $0 <R_sphere_km> <input_in_file> [xyz_output_file] [refine_factor] [source_nlat] [source_nlon] [source_nr]"
  exit 2
fi

if [[ ! -f "${INPUT_FILE}" ]]; then
  echo "ERROR: input file not found: ${INPUT_FILE}"
  exit 2
fi

if [[ -z "${XYZ_OUT}" ]]; then
  base="$(basename "${INPUT_FILE}")"
  base="${base%.in}"
  # keep derived inputs in their owning case's output directory
  output_dir="$(cd -- "$(dirname -- "${INPUT_FILE}")" && pwd)"
  parent_dir="${output_dir}"
  while [[ "${parent_dir}" != "/" ]]; do
    if [[ "${parent_dir##*/}" == "output" || "${parent_dir##*/}" == "figs" ]]; then
      output_dir="$(dirname -- "${parent_dir}")"
      break
    fi
    parent_dir="$(dirname -- "${parent_dir}")"
  done
  output_dir="${output_dir%/}/output"
  XYZ_OUT="${output_dir}/${base}_xyz.txt"
fi
mkdir -p -- "$(dirname -- "${XYZ_OUT}")"

"${ROOT_DIR}/build_gravmag_tools.sh"

echo "+ running point_source solver"
echo "+ ./gravmag_sphere_point_source ${R_SPHERE_KM} ${INPUT_FILE} ${XYZ_OUT} ${REFINE_FACTOR} ${SOURCE_NLAT} ${SOURCE_NLON} ${SOURCE_NR}"
"${ROOT_DIR}/gravmag_sphere_point_source" "${R_SPHERE_KM}" "${INPUT_FILE}" "${XYZ_OUT}" "${REFINE_FACTOR}" "${SOURCE_NLAT}" "${SOURCE_NLON}" "${SOURCE_NR}"

echo "Done: ${XYZ_OUT}"
