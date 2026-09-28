#!/usr/bin/env bash
set -euo pipefail

# run GravMagSphere comparison pipelines end-to-end
#
# Covered test groups:
#   - external:  Fortran spectral vs SciPy SH LSQ vs SHTOOLS LSQ (shared residual suite)
#   - shtools:   dedicated GravMagSphere vs SHTOOLS residual/side-by-side outputs
#   - harmonica: GravMagSphere spectral vs Harmonica EQS in XYZ
#   - complexlarge: direct_sum vs spectral boundary-visibility analysis
#   - hybrid:    hybrid-vs-spectral vertex sweep
#
# Notes:
#   - uses python scripts in same directory
#   - can auto-clean intermediate text/csv artifacts

diagnostics_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${diagnostics_dir}/.." && pwd)"
example_root="${repo_root}/examples"

PYTHON_BIN="${PYTHON_BIN:-python3}"
TESTS_RAW="external,shtools,harmonica,complexlarge,hybrid,quadrature"
EXAMPLES_RAW=""

SKIP_BUILD=0
CLEAN_INTERMEDIATE=1
STRICT_DEPS=0

RSPHERE_KM="${RSPHERE_KM:-1737.4}"
REFINE_FACTOR="${REFINE_FACTOR:-2}"
FORTRAN_AUTO_MODE="${FORTRAN_AUTO_MODE:-1}"
LMAX="${LMAX:-24}"
NTHETA_FIT="${NTHETA_FIT:-72}"
NPHI_FIT="${NPHI_FIT:-144}"
REG_LAMBDA="${REG_LAMBDA:-0.2}"
REG_POWER="${REG_POWER:-4.0}"

SCIPY_LMAX_GRID="${SCIPY_LMAX_GRID:-12,18,24,30,36}"
SCIPY_RIDGE_GRID="${SCIPY_RIDGE_GRID:-0,1e-10,1e-8,1e-6}"
SHTOOLS_LMAX_GRID="${SHTOOLS_LMAX_GRID:-12,18,24,30,36}"

HARMONICA_DAMPING_GRID="${HARMONICA_DAMPING_GRID:-0,1e-8,1e-6,1e-4}"
HARMONICA_DEPTH_FACTOR_GRID="${HARMONICA_DEPTH_FACTOR_GRID:-1,2,4}"
HARMONICA_BLOCK_FACTOR="${HARMONICA_BLOCK_FACTOR:-30}"

HYBRID_VERTICES="${HYBRID_VERTICES:-4,6,8,10,12,14,16,20,24,30,36,48}"
HYBRID_SHAPES="${HYBRID_SHAPES:-regular,wavy}"

COMPLEXLARGE_DIR="diagnostics/complexlarge_direct_vs_spectral"
FORMATTED_ROOT="diagnostics/external_solver_inputs"

usage() {
  cat <<'EOF'
Usage: bash diagnostics/run_comparison_tests.sh [options]

Options:
  --python <bin>             Python executable (default: python3 or $PYTHON_BIN)
  --tests <csv>              Test groups to run.
                             Available: external,shtools,harmonica,complexlarge,hybrid,quadrature
                             Default: external,shtools,harmonica,complexlarge,hybrid,quadrature
  --examples <csv>           Optional example filters (e.g. "case1,case2" or "case1.in,case2.in")
                             Applied to external/shtools/harmonica workflows.
  --skip-build               Do not rebuild Fortran tools before running tests
  --no-clean                 Keep intermediate txt/csv workspace artifacts
  --strict-deps              Fail immediately when Python dependencies are missing
  --help                     Show this help text

Important env overrides:
  RSPHERE_KM, REFINE_FACTOR, FORTRAN_AUTO_MODE,
  LMAX, NTHETA_FIT, NPHI_FIT, REG_LAMBDA, REG_POWER,
  SCIPY_LMAX_GRID, SCIPY_RIDGE_GRID, SHTOOLS_LMAX_GRID,
  HARMONICA_DAMPING_GRID, HARMONICA_DEPTH_FACTOR_GRID, HARMONICA_BLOCK_FACTOR,
  HYBRID_VERTICES, HYBRID_SHAPES

Examples:
  bash diagnostics/run_comparison_tests.sh
  bash diagnostics/run_comparison_tests.sh --tests external,shtools --examples gravmag_sphere_1body_mag_polygon_inc90_dec0_base
  PYTHON_BIN=/Users/danywaller/code/venvs/gravmagpy/bin/python bash diagnostics/run_comparison_tests.sh --tests harmonica,hybrid
EOF
}

log() {
  printf '[compare-tests] %s\n' "$*"
}

have_test() {
  local q="$1"
  local t
  for t in "${TESTS[@]}"; do
    [[ "$t" == "$q" ]] && return 0
  done
  return 1
}

python_has_module() {
  local mod="$1"
  "$PYTHON_BIN" -c "import importlib.util,sys; sys.exit(0 if importlib.util.find_spec('$mod') else 1)" >/dev/null 2>&1
}

require_modules() {
  local missing=()
  local mod
  for mod in "$@"; do
    if ! python_has_module "$mod"; then
      missing+=("$mod")
    fi
  done
  if ((${#missing[@]} > 0)); then
    if [[ "$STRICT_DEPS" == "1" ]]; then
      printf 'ERROR: missing Python module(s): %s\n' "${missing[*]}" >&2
      exit 2
    fi
    log "Skipping test due to missing Python module(s): ${missing[*]}"
    return 1
  fi
  return 0
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --python)
      PYTHON_BIN="${2:-}"
      shift 2
      ;;
    --tests)
      TESTS_RAW="${2:-}"
      shift 2
      ;;
    --examples|--lunar_examples)
      EXAMPLES_RAW="${2:-}"
      shift 2
      ;;
    --skip-build)
      SKIP_BUILD=1
      shift
      ;;
    --no-clean)
      CLEAN_INTERMEDIATE=0
      shift
      ;;
    --strict-deps)
      STRICT_DEPS=1
      shift
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 2
      ;;
  esac
done

cd "$repo_root"

IFS=',' read -r -a TESTS <<<"$TESTS_RAW"
if ((${#TESTS[@]} == 0)); then
  echo "ERROR: no tests selected." >&2
  exit 2
fi

for t in "${TESTS[@]}"; do
  case "$t" in
    external|shtools|harmonica|complexlarge|hybrid|quadrature) ;;
    *)
      echo "ERROR: unknown test group '$t'." >&2
      echo "Valid values: external,shtools,harmonica,complexlarge,hybrid,quadrature" >&2
      exit 2
      ;;
  esac
done

FORMAT_EXAMPLES=()
COMPARE_CASES=()
if [[ -n "$EXAMPLES_RAW" ]]; then
  IFS=',' read -r -a RAW_EXAMPLES <<<"$EXAMPLES_RAW"
  for ex in "${RAW_EXAMPLES[@]}"; do
    ex="${ex//[[:space:]]/}"
    [[ -z "$ex" ]] && continue

    local_name="$ex"
    if [[ "$local_name" != *.in && -f "$example_root/examples/${local_name}.in" ]]; then
      local_name="${local_name}.in"
    fi

    FORMAT_EXAMPLES+=("$local_name")
    COMPARE_CASES+=("${local_name%.in}")
  done
fi

if [[ "$SKIP_BUILD" == "0" ]]; then
  log "Building Fortran tools"
  "$example_root/build_gravmag_tools.sh"
else
  log "Skipping build (--skip-build)"
fi

run_external_format() {
  local cmd=(
    "$PYTHON_BIN" diagnostics/external_solver_compare.py format
    --lunar_examples-dir examples/lunar_examples
    --out-root "$FORMATTED_ROOT"
    --rsphere-km "$RSPHERE_KM"
    --refine-factor "$REFINE_FACTOR"
  )
  if ((${#FORMAT_EXAMPLES[@]} > 0)); then
    cmd+=(--lunar_examples "${FORMAT_EXAMPLES[@]}")
  fi
  log "Formatting comparison inputs"
  "${cmd[@]}"
}

run_compare_external() {
  local tmp_csv="/tmp/gravmag_external_solver_comparison.csv"
  local cmd=(
    "$PYTHON_BIN" diagnostics/external_solver_compare.py compare
    --formatted-root "$FORMATTED_ROOT"
    --report-md diagnostics/external_solver_comparison.md
    --report-csv "$tmp_csv"
    --residual-dir diagnostics/external_solver_residuals
    --side-by-side-dir diagnostics/gravmag_vs_shtools_side_by_side
    --fortran-auto-mode "$FORTRAN_AUTO_MODE"
    --rsphere-km "$RSPHERE_KM"
    --refine-factor "$REFINE_FACTOR"
    --lmax "$LMAX"
    --ntheta-fit "$NTHETA_FIT"
    --nphi-fit "$NPHI_FIT"
    --reg-lambda "$REG_LAMBDA"
    --reg-power "$REG_POWER"
    --scipy-lmax-grid "$SCIPY_LMAX_GRID"
    --scipy-ridge-grid "$SCIPY_RIDGE_GRID"
    --shtools-lmax-grid "$SHTOOLS_LMAX_GRID"
  )
  if ((${#COMPARE_CASES[@]} > 0)); then
    cmd+=(--lunar_examples "${COMPARE_CASES[@]}")
  fi
  log "Running external solver comparison (Fortran vs SciPy vs SHTOOLS)"
  "${cmd[@]}"
  [[ "$CLEAN_INTERMEDIATE" == "1" ]] && rm -f "$tmp_csv"
}

run_compare_shtools() {
  local tmp_csv="/tmp/gravmag_vs_shtools_residuals.csv"
  local cmd=(
    "$PYTHON_BIN" diagnostics/external_solver_compare.py compare
    --formatted-root "$FORMATTED_ROOT"
    --report-md diagnostics/gravmag_vs_shtools_residuals.md
    --report-csv "$tmp_csv"
    --residual-dir diagnostics/gravmag_vs_shtools_residuals
    --side-by-side-dir diagnostics/gravmag_vs_shtools_side_by_side
    --fortran-auto-mode "$FORTRAN_AUTO_MODE"
    --rsphere-km "$RSPHERE_KM"
    --refine-factor "$REFINE_FACTOR"
    --lmax "$LMAX"
    --ntheta-fit "$NTHETA_FIT"
    --nphi-fit "$NPHI_FIT"
    --reg-lambda "$REG_LAMBDA"
    --reg-power "$REG_POWER"
    --scipy-lmax-grid "$SCIPY_LMAX_GRID"
    --scipy-ridge-grid "$SCIPY_RIDGE_GRID"
    --shtools-lmax-grid "$SHTOOLS_LMAX_GRID"
  )
  if ((${#COMPARE_CASES[@]} > 0)); then
    cmd+=(--lunar_examples "${COMPARE_CASES[@]}")
  fi
  log "Running dedicated GravMagSphere vs SHTOOLS residual suite"
  "${cmd[@]}"
  [[ "$CLEAN_INTERMEDIATE" == "1" ]] && rm -f "$tmp_csv"
}

run_compare_harmonica() {
  local tmp_csv="/tmp/gravmag_vs_harmonica_xyz.csv"
  local cmd=(
    "$PYTHON_BIN" diagnostics/external_solver_compare.py compare-harmonica-xyz
    --formatted-root "$FORMATTED_ROOT"
    --report-md diagnostics/gravmag_vs_harmonica_xyz.md
    --report-csv "$tmp_csv"
    --residual-dir diagnostics/gravmag_vs_harmonica_xyz_residuals
    --side-by-side-dir diagnostics/gravmag_vs_harmonica_xyz_side_by_side
    --fortran-auto-mode "$FORTRAN_AUTO_MODE"
    --rsphere-km "$RSPHERE_KM"
    --refine-factor "$REFINE_FACTOR"
    --lmax "$LMAX"
    --ntheta-fit "$NTHETA_FIT"
    --nphi-fit "$NPHI_FIT"
    --reg-lambda "$REG_LAMBDA"
    --reg-power "$REG_POWER"
    --harmonica-damping-grid "$HARMONICA_DAMPING_GRID"
    --harmonica-depth-factor-grid "$HARMONICA_DEPTH_FACTOR_GRID"
    --harmonica-block-factor "$HARMONICA_BLOCK_FACTOR"
  )
  if ((${#COMPARE_CASES[@]} > 0)); then
    cmd+=(--lunar_examples "${COMPARE_CASES[@]}")
  fi
  log "Running GravMagSphere vs Harmonica XYZ comparison"
  "${cmd[@]}"
  [[ "$CLEAN_INTERMEDIATE" == "1" ]] && rm -f "$tmp_csv"
}

run_complexlarge_analysis() {
  local run_dir="$repo_root/$COMPLEXLARGE_DIR"
  mkdir -p "$run_dir"

  log "Generating complexlarge direct_sum baseline"
  SKIP_BUILD=1 "$example_root/run_input_to_xyz.sh" direct_sum "$RSPHERE_KM" \
    "$example_root/examples/gravmag_sphere_1body_mag_polygon_inc30_dec210_complexlarge.in" \
    "$run_dir/complexlarge_direct_xyz.txt" \
    --refine-factor "$REFINE_FACTOR"

  log "Generating complexlarge spectral (local correction on, hybrid off)"
  SKIP_BUILD=1 "$example_root/run_input_to_xyz.sh" spectral "$RSPHERE_KM" \
    "$example_root/examples/gravmag_sphere_1body_mag_polygon_inc30_dec210_complexlarge.in" \
    "$run_dir/complexlarge_spectral_xyz.txt" \
    --lmax "$LMAX" \
    --refine-factor "$REFINE_FACTOR" \
    --ntheta-fit "$NTHETA_FIT" \
    --nphi-fit "$NPHI_FIT" \
    --reg-lambda "$REG_LAMBDA" \
    --reg-power "$REG_POWER" \
    --auto-mode "$FORTRAN_AUTO_MODE" \
    --edge-correction 1 \
    --hybrid-mode 0

  log "Generating complexlarge spectral (local correction off, hybrid off)"
  SKIP_BUILD=1 "$example_root/run_input_to_xyz.sh" spectral "$RSPHERE_KM" \
    "$example_root/examples/gravmag_sphere_1body_mag_polygon_inc30_dec210_complexlarge.in" \
    "$run_dir/complexlarge_spectral_nolocal_xyz.txt" \
    --lmax "$LMAX" \
    --refine-factor "$REFINE_FACTOR" \
    --ntheta-fit "$NTHETA_FIT" \
    --nphi-fit "$NPHI_FIT" \
    --reg-lambda "$REG_LAMBDA" \
    --reg-power "$REG_POWER" \
    --auto-mode "$FORTRAN_AUTO_MODE" \
    --no-local-correction \
    --hybrid-mode 0

  log "Running complexlarge boundary-visibility analysis"
  "$PYTHON_BIN" diagnostics/complexlarge_direct_vs_spectral_analysis.py --run-dir "$COMPLEXLARGE_DIR"

  if [[ "$CLEAN_INTERMEDIATE" == "1" ]]; then
    rm -f \
      "$run_dir/complexlarge_direct_xyz.txt" \
      "$run_dir/complexlarge_spectral_xyz.txt" \
      "$run_dir/complexlarge_spectral_nolocal_xyz.txt"
  fi
}

run_hybrid_sweep() {
  log "Running hybrid vertex sweep"
  "$PYTHON_BIN" diagnostics/hybrid_vertex_sweep.py \
    --output-dir diagnostics/hybrid_vertex_sweep \
    --vertices "$HYBRID_VERTICES" \
    --shapes "$HYBRID_SHAPES" \
    --skip-build

  if [[ "$CLEAN_INTERMEDIATE" == "1" ]]; then
    rm -f diagnostics/hybrid_vertex_sweep/*.csv
    rm -rf diagnostics/_hybrid_vertex_sweep_examples diagnostics/_hybrid_vertex_sweep_workspace
  fi
}

need_formatted=0
if have_test external || have_test shtools || have_test harmonica; then
  need_formatted=1
fi

if [[ "$need_formatted" == "1" ]]; then
  if require_modules numpy scipy; then
    run_external_format
  else
    if [[ "$STRICT_DEPS" == "1" ]]; then
      exit 2
    fi
  fi
fi

if have_test external; then
  if require_modules numpy scipy matplotlib pyshtools; then
    run_compare_external
  fi
fi

if have_test shtools; then
  if require_modules numpy scipy matplotlib pyshtools; then
    run_compare_shtools
  fi
fi

if have_test harmonica; then
  if require_modules numpy scipy matplotlib harmonica; then
    run_compare_harmonica
  fi
fi

if have_test complexlarge; then
  if require_modules numpy matplotlib; then
    run_complexlarge_analysis
  fi
fi

if have_test hybrid; then
  if require_modules numpy matplotlib; then
    run_hybrid_sweep
  fi
fi

if have_test quadrature; then
  if require_modules numpy scipy matplotlib; then
    log "Running gauss-legendre, direct_sum and pure spectral comparison"
    "$PYTHON_BIN" diagnostics/gauss_legendre_compare.py
  fi
fi

log "Completed selected comparison tests: ${TESTS_RAW}"
