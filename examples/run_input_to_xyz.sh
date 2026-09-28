#!/usr/bin/env bash
set -euo pipefail

# Step 1: input .in -> solver XYZ output.
#
# Usage:
#   ./run_input_to_xyz.sh <solver> <R_sphere_km> <input_in_file> \
#     [xyz_output_file] [refine_factor] [lmax] [ntheta_fit] [nphi_fit] \
#     [reg_lambda] [reg_power] [source_nlat] [source_nlon] [source_nr] \
#     [auto_mode] [joint_strength] [edge_correction] [hybrid_mode] \
#     [hybrid_band_deg] [complex_vertex_threshold] [hybrid_transition_deg]
#   Optional args:
#     --refine-factor N --lmax N --ntheta-fit N --nphi-fit N
#     --reg-lambda X --reg-power P
#     --source-nlat N --source-nlon N --source-nr N
#     --auto-mode 0|1 --joint-strength X --edge-correction 0|1 --no-local-correction
#     --hybrid-mode 0|1|2 --hybrid-band-deg X --complex-vertex-threshold N --hybrid-transition-deg X
#
# solver:
#   - spectral : gravmag_sphere_spectral
#   - point_source   : gravmag_sphere_point_source
#   - gauss_legendre: gravmag_sphere_gauss_legendre
#     source_nr/source_nlat/source_nlon set node orders; refine_factor sets subdivisions
#
# output:
#   - writes one XYZ table (body_id lon lat Fx Fy Fz Ftot)
#   - does not perform spherical component conversion.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

SOLVER="${1:-spectral}"

# accept legacy selectors while dispatching canonical method names
case "${SOLVER}" in
  direct) SOLVER="point_source" ;;
  gauss) SOLVER="spectral" ;;
esac
R_SPHERE_KM="${2:-}"
INPUT_FILE="${3:-}"
shift $(( $# >= 3 ? 3 : $# ))

XYZ_OUT=""
REFINE_FACTOR="2"
LMAX="24"
NTHETA_FIT="72"
NPHI_FIT="144"
REG_LAMBDA="0.2"
REG_POWER="4.0"
SOURCE_NLAT="0"
SOURCE_NLON="0"
SOURCE_NR="0"
AUTO_MODE="1"
JOINT_STRENGTH="1.0"
EDGE_CORRECTION="1"
HYBRID_MODE="1"
HYBRID_BAND_DEG="1.5"
COMPLEX_VERTEX_THRESHOLD="12"
HYBRID_TRANSITION_DEG="0.75"
SPECTRAL_FORCE_FULL="0"
point_source_force_full="0"

positional=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --no-local-correction)
      EDGE_CORRECTION="0"
      SPECTRAL_FORCE_FULL="1"
      shift
      ;;
    --refine-factor)
      REFINE_FACTOR="${2:-}"
      SPECTRAL_FORCE_FULL="1"
      point_source_force_full="1"
      shift 2
      ;;
    --lmax)
      LMAX="${2:-}"
      SPECTRAL_FORCE_FULL="1"
      shift 2
      ;;
    --ntheta-fit)
      NTHETA_FIT="${2:-}"
      SPECTRAL_FORCE_FULL="1"
      shift 2
      ;;
    --nphi-fit)
      NPHI_FIT="${2:-}"
      SPECTRAL_FORCE_FULL="1"
      shift 2
      ;;
    --reg-lambda)
      REG_LAMBDA="${2:-}"
      SPECTRAL_FORCE_FULL="1"
      shift 2
      ;;
    --reg-power)
      REG_POWER="${2:-}"
      SPECTRAL_FORCE_FULL="1"
      shift 2
      ;;
    --source-nlat)
      SOURCE_NLAT="${2:-}"
      SPECTRAL_FORCE_FULL="1"
      point_source_force_full="1"
      shift 2
      ;;
    --source-nlon)
      SOURCE_NLON="${2:-}"
      SPECTRAL_FORCE_FULL="1"
      point_source_force_full="1"
      shift 2
      ;;
    --source-nr)
      SOURCE_NR="${2:-}"
      SPECTRAL_FORCE_FULL="1"
      point_source_force_full="1"
      shift 2
      ;;
    --auto-mode)
      AUTO_MODE="${2:-}"
      SPECTRAL_FORCE_FULL="1"
      shift 2
      ;;
    --joint-strength)
      JOINT_STRENGTH="${2:-}"
      SPECTRAL_FORCE_FULL="1"
      shift 2
      ;;
    --edge-correction)
      EDGE_CORRECTION="${2:-}"
      SPECTRAL_FORCE_FULL="1"
      shift 2
      ;;
    --hybrid-mode)
      HYBRID_MODE="${2:-}"
      SPECTRAL_FORCE_FULL="1"
      shift 2
      ;;
    --hybrid-band-deg)
      HYBRID_BAND_DEG="${2:-}"
      SPECTRAL_FORCE_FULL="1"
      shift 2
      ;;
    --complex-vertex-threshold)
      COMPLEX_VERTEX_THRESHOLD="${2:-}"
      SPECTRAL_FORCE_FULL="1"
      shift 2
      ;;
    --hybrid-transition-deg)
      HYBRID_TRANSITION_DEG="${2:-}"
      SPECTRAL_FORCE_FULL="1"
      shift 2
      ;;
    *)
      positional+=("$1")
      shift
      ;;
  esac
done

if [[ ${#positional[@]} -ge 1 ]]; then XYZ_OUT="${positional[0]}"; fi
if [[ ${#positional[@]} -ge 2 ]]; then
  REFINE_FACTOR="${positional[1]}"
  SPECTRAL_FORCE_FULL="1"
  point_source_force_full="1"
fi
if [[ ${#positional[@]} -ge 3 ]]; then LMAX="${positional[2]}"; SPECTRAL_FORCE_FULL="1"; fi
if [[ ${#positional[@]} -ge 4 ]]; then NTHETA_FIT="${positional[3]}"; SPECTRAL_FORCE_FULL="1"; fi
if [[ ${#positional[@]} -ge 5 ]]; then NPHI_FIT="${positional[4]}"; SPECTRAL_FORCE_FULL="1"; fi
if [[ ${#positional[@]} -ge 6 ]]; then REG_LAMBDA="${positional[5]}"; SPECTRAL_FORCE_FULL="1"; fi
if [[ ${#positional[@]} -ge 7 ]]; then REG_POWER="${positional[6]}"; SPECTRAL_FORCE_FULL="1"; fi
if [[ ${#positional[@]} -ge 8 ]]; then SOURCE_NLAT="${positional[7]}"; SPECTRAL_FORCE_FULL="1"; point_source_force_full="1"; fi
if [[ ${#positional[@]} -ge 9 ]]; then SOURCE_NLON="${positional[8]}"; SPECTRAL_FORCE_FULL="1"; point_source_force_full="1"; fi
if [[ ${#positional[@]} -ge 10 ]]; then SOURCE_NR="${positional[9]}"; SPECTRAL_FORCE_FULL="1"; point_source_force_full="1"; fi
if [[ ${#positional[@]} -ge 11 ]]; then AUTO_MODE="${positional[10]}"; SPECTRAL_FORCE_FULL="1"; fi
if [[ ${#positional[@]} -ge 12 ]]; then JOINT_STRENGTH="${positional[11]}"; SPECTRAL_FORCE_FULL="1"; fi
if [[ ${#positional[@]} -ge 13 ]]; then EDGE_CORRECTION="${positional[12]}"; SPECTRAL_FORCE_FULL="1"; fi
if [[ ${#positional[@]} -ge 14 ]]; then HYBRID_MODE="${positional[13]}"; SPECTRAL_FORCE_FULL="1"; fi
if [[ ${#positional[@]} -ge 15 ]]; then HYBRID_BAND_DEG="${positional[14]}"; SPECTRAL_FORCE_FULL="1"; fi
if [[ ${#positional[@]} -ge 16 ]]; then COMPLEX_VERTEX_THRESHOLD="${positional[15]}"; SPECTRAL_FORCE_FULL="1"; fi
if [[ ${#positional[@]} -ge 17 ]]; then HYBRID_TRANSITION_DEG="${positional[16]}"; SPECTRAL_FORCE_FULL="1"; fi

if [[ -z "${R_SPHERE_KM}" || -z "${INPUT_FILE}" ]]; then
  echo "Usage: $0 <solver> <R_sphere_km> <input_in_file> [xyz_output_file] [refine_factor] [lmax] [ntheta_fit] [nphi_fit] [reg_lambda] [reg_power] [source_nlat] [source_nlon] [source_nr] [auto_mode] [joint_strength] [edge_correction] [hybrid_mode] [hybrid_band_deg] [complex_vertex_threshold] [hybrid_transition_deg] [--refine-factor N] [--lmax N] [--ntheta-fit N] [--nphi-fit N] [--reg-lambda X] [--reg-power P] [--source-nlat N] [--source-nlon N] [--source-nr N] [--auto-mode 0|1] [--joint-strength X] [--edge-correction 0|1] [--hybrid-mode 0|1|2] [--hybrid-band-deg X] [--complex-vertex-threshold N] [--hybrid-transition-deg X] [--no-local-correction]"
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

if [[ "${SKIP_BUILD:-0}" != "1" ]]; then
  # build once unless caller explicitly skipped it.
  "${ROOT_DIR}/build_gravmag_tools.sh"
fi

case "${SOLVER}" in
  spectral)
    # spectral pipeline: source mesh -> harmonic fit -> XYZ field table.
    if [[ "${SPECTRAL_FORCE_FULL}" == "0" ]]; then
      echo "+ step1 spectral(auto): ./gravmag_sphere_spectral ${R_SPHERE_KM} ${INPUT_FILE} ${XYZ_OUT}"
      "${ROOT_DIR}/gravmag_sphere_spectral" \
        "${R_SPHERE_KM}" "${INPUT_FILE}" "${XYZ_OUT}"
    else
      echo "+ step1 spectral: ./gravmag_sphere_spectral ${R_SPHERE_KM} ${INPUT_FILE} ${XYZ_OUT} ${LMAX} ${REFINE_FACTOR} ${NTHETA_FIT} ${NPHI_FIT} ${REG_LAMBDA} ${REG_POWER} ${SOURCE_NLAT} ${SOURCE_NLON} ${SOURCE_NR} ${AUTO_MODE} ${JOINT_STRENGTH} ${EDGE_CORRECTION} ${HYBRID_MODE} ${HYBRID_BAND_DEG} ${COMPLEX_VERTEX_THRESHOLD} ${HYBRID_TRANSITION_DEG}"
      "${ROOT_DIR}/gravmag_sphere_spectral" \
        "${R_SPHERE_KM}" "${INPUT_FILE}" "${XYZ_OUT}" \
        "${LMAX}" "${REFINE_FACTOR}" "${NTHETA_FIT}" "${NPHI_FIT}" \
        "${REG_LAMBDA}" "${REG_POWER}" "${SOURCE_NLAT}" "${SOURCE_NLON}" "${SOURCE_NR}" \
        "${AUTO_MODE}" "${JOINT_STRENGTH}" "${EDGE_CORRECTION}" \
        "${HYBRID_MODE}" "${HYBRID_BAND_DEG}" "${COMPLEX_VERTEX_THRESHOLD}" "${HYBRID_TRANSITION_DEG}"
    fi
    ;;
  point_source)
    # point_source pipeline: source mesh -> point_source Green's function summation.
    if [[ "${point_source_force_full}" == "0" ]]; then
      echo "+ step1 point_source(auto): ./gravmag_sphere_point_source ${R_SPHERE_KM} ${INPUT_FILE} ${XYZ_OUT}"
      "${ROOT_DIR}/gravmag_sphere_point_source" \
        "${R_SPHERE_KM}" "${INPUT_FILE}" "${XYZ_OUT}"
    else
      echo "+ step1 point_source: ./gravmag_sphere_point_source ${R_SPHERE_KM} ${INPUT_FILE} ${XYZ_OUT} ${REFINE_FACTOR} ${SOURCE_NLAT} ${SOURCE_NLON} ${SOURCE_NR}"
      "${ROOT_DIR}/gravmag_sphere_point_source" \
        "${R_SPHERE_KM}" "${INPUT_FILE}" "${XYZ_OUT}" \
        "${REFINE_FACTOR}" "${SOURCE_NLAT}" "${SOURCE_NLON}" "${SOURCE_NR}"
    fi
    ;;
  gauss_legendre)
    "${ROOT_DIR}/gravmag_sphere_gauss_legendre" \
      "${R_SPHERE_KM}" "${INPUT_FILE}" "${XYZ_OUT}" \
      "${SOURCE_NR}" "${SOURCE_NLAT}" "${SOURCE_NLON}" "${REFINE_FACTOR}"
    ;;
  *)
    echo "ERROR: solver must be 'spectral', 'point_source', or 'gauss_legendre' (got: ${SOLVER})"
    exit 2
    ;;
esac

echo "Done: ${XYZ_OUT}"
