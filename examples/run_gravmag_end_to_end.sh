#!/usr/bin/env bash
set -euo pipefail

# End-to-end pipeline:
#   input .in -> point_source solver XYZ components -> spherical Br/Btheta/Bphi output
#
# Usage:
#   ./run_gravmag_end_to_end.sh <R_sphere_km> <input_in_file> [xyz_output_file] [spherical_output_file] [refine_factor] [source_nlat] [source_nlon] [source_nr]
#
# This is the new clean "single-case point_source solver" wrapper

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_SPHERE_KM="${1:-}"
INPUT_FILE="${2:-}"
XYZ_OUT="${3:-}"
SPH_OUT="${4:-}"
REFINE_FACTOR="${5:-2}"
SOURCE_NLAT="${6:-0}"
SOURCE_NLON="${7:-0}"
SOURCE_NR="${8:-0}"

if [[ -z "${R_SPHERE_KM}" || -z "${INPUT_FILE}" ]]; then
  echo "Usage: $0 <R_sphere_km> <input_in_file> [xyz_output_file] [spherical_output_file] [refine_factor] [source_nlat] [source_nlon] [source_nr]"
  exit 2
fi

if [[ ! -f "${INPUT_FILE}" ]]; then
  echo "ERROR: input file not found: ${INPUT_FILE}"
  exit 2
fi

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

if [[ -z "${XYZ_OUT}" ]]; then
  XYZ_OUT="${output_dir}/${base}_xyz.txt"
fi

if [[ -z "${SPH_OUT}" ]]; then
  SPH_OUT="${output_dir}/${base}_brtp.txt"
fi
mkdir -p -- "$(dirname -- "${XYZ_OUT}")" "$(dirname -- "${SPH_OUT}")"

"${ROOT_DIR}/build_gravmag_tools.sh"

echo "+ step 1/2: point_source solver (xyz)"
echo "+ ./gravmag_sphere_point_source ${R_SPHERE_KM} ${INPUT_FILE} ${XYZ_OUT} ${REFINE_FACTOR} ${SOURCE_NLAT} ${SOURCE_NLON} ${SOURCE_NR}"
"${ROOT_DIR}/gravmag_sphere_point_source" "${R_SPHERE_KM}" "${INPUT_FILE}" "${XYZ_OUT}" "${REFINE_FACTOR}" "${SOURCE_NLAT}" "${SOURCE_NLON}" "${SOURCE_NR}"

echo "+ step 2/2: xyz -> spherical conversion"
echo "+ ./gravmag_xyz_to_brtp ${XYZ_OUT} ${SPH_OUT}"
"${ROOT_DIR}/gravmag_xyz_to_brtp" "${XYZ_OUT}" "${SPH_OUT}"

echo "Done."
echo "  XYZ output      : ${XYZ_OUT}"
echo "  Spherical output: ${SPH_OUT}"
