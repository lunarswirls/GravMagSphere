#!/usr/bin/env bash
set -euo pipefail

# use zero orders to inherit nr, ntheta and nphi from each body's card 3
if [[ $# -lt 2 || $# -gt 7 ]]; then
  echo "usage: $0 radius_km input.in [output.txt] [radial_order] [latitude_order] [longitude_order] [subdivisions]"
  exit 2
fi
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ "${SKIP_BUILD:-0}" != "1" ]]; then
  bash "${script_dir}/build_gravmag_tools.sh"
fi
"${script_dir}/gravmag_sphere_gauss_legendre" "$@"
