#!/usr/bin/env bash
set -euo pipefail

# Run all example inputs end-to-end in bash:
#   input .in -> XYZ -> Br/Btheta/Bphi
#
# Usage:
#   ./run_all_examples_end_to_end.sh [solver] [R_sphere_km] [examples_dir] [output_dir] \
#     [refine_factor] [lmax] [ntheta_fit] [nphi_fit] [reg_lambda] [reg_power] \
#     [source_nlat] [source_nlon] [source_nr]
#
# solver options:
#   - spectral (default)
#   - direct_sum
#
# method:
#   - compiles tools once
#   - runs every example .in (except "*comment*" if they exist)
#   - writes xyz and brtp text outputs

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

SOLVER="${1:-spectral}"

# accept legacy selectors while dispatching canonical method names
case "${SOLVER}" in
  direct) SOLVER="direct_sum" ;;
  gauss) SOLVER="spectral" ;;
esac
R_SPHERE_KM="${2:-1737.4}"
EXAMPLES_DIR="${3:-${ROOT_DIR}/lunar_examples}"
OUTPUT_DIR="${4:-}"
REFINE_FACTOR="${5:-2}"
LMAX="${6:-24}"
NTHETA_FIT="${7:-72}"
NPHI_FIT="${8:-144}"
REG_LAMBDA="${9:-0.2}"
REG_POWER="${10:-4.0}"
SOURCE_NLAT="${11:-0}"
SOURCE_NLON="${12:-0}"
SOURCE_NR="${13:-0}"

if [[ "${SOLVER}" != "spectral" && "${SOLVER}" != "direct_sum" && "${SOLVER}" != "gauss_legendre" ]]; then
  echo "ERROR: solver must be 'spectral', 'direct_sum', or 'gauss_legendre' (got: ${SOLVER})"
  exit 2
fi

if [[ "${EXAMPLES_DIR}" != /* ]]; then EXAMPLES_DIR="${ROOT_DIR}/${EXAMPLES_DIR}"; fi
if [[ -n "${OUTPUT_DIR}" && "${OUTPUT_DIR}" != /* ]]; then OUTPUT_DIR="${ROOT_DIR}/${OUTPUT_DIR}"; fi
if [[ ! -d "${EXAMPLES_DIR}" ]]; then
  echo "ERROR: examples directory not found: ${EXAMPLES_DIR}"
  exit 2
fi

EXAMPLES_DIR="$(cd -- "${EXAMPLES_DIR}" && pwd)"
OUTPUT_DIR="${OUTPUT_DIR:-${EXAMPLES_DIR}/output}"
mkdir -p -- "${OUTPUT_DIR}"
OUTPUT_DIR="$(cd -- "${OUTPUT_DIR}" && pwd)"

CASES="$(find "${EXAMPLES_DIR}" -maxdepth 1 -type f -name "*.in" | sort)"
if [[ -z "${CASES}" ]]; then
  echo "ERROR: no .in example files found in ${EXAMPLES_DIR}"
  exit 2
fi

"${ROOT_DIR}/build_gravmag_tools.sh"

count=0
while IFS= read -r input_file; do
  [[ -z "${input_file}" ]] && continue
  name="$(basename "${input_file}")"
  low_name="$(echo "${name}" | tr '[:upper:]' '[:lower:]')"
  if [[ "${low_name}" == *comment* ]]; then
    continue
  fi

  base="${name%.in}"
  xyz_out="${OUTPUT_DIR}/${base}_xyz.txt"
  brtp_out="${OUTPUT_DIR}/${base}_brtp.txt"

  echo "[RUN] ${name} (${SOLVER})"
  # step 1: input cards -> XYZ field table
  SKIP_BUILD=1 "${ROOT_DIR}/run_input_to_xyz.sh" "${SOLVER}" "${R_SPHERE_KM}" "${input_file}" "${xyz_out}" \
    "${REFINE_FACTOR}" "${LMAX}" "${NTHETA_FIT}" "${NPHI_FIT}" "${REG_LAMBDA}" "${REG_POWER}" \
    "${SOURCE_NLAT}" "${SOURCE_NLON}" "${SOURCE_NR}"

  # step 2: XYZ -> Br/Btheta/Bphi
  SKIP_BUILD=1 "${ROOT_DIR}/run_xyz_to_brtp.sh" "${xyz_out}" "${brtp_out}"
  echo "[OK ] wrote ${xyz_out} and ${brtp_out}"
  count=$((count + 1))
done <<< "${CASES}"

echo "Completed ${count} example runs."
