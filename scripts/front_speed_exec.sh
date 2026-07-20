#!/usr/bin/env bash
# GitHub reproducibility helper for the redox-fungicide manuscript.
# Place this file in the repository `scripts/` directory and run it from anywhere.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

RESULTS_ROOT="${RESULTS_ROOT:-${REPO_ROOT}/results}"
mkdir -p "${RESULTS_ROOT}"
case "${RESULTS_ROOT}" in
  */) OUTPUT_ROOT="${RESULTS_ROOT}" ;;
  *)  OUTPUT_ROOT="${RESULTS_ROOT}/" ;;
esac

MACHINE_LABEL="${MACHINE_LABEL:-local}"
OS="$(uname -s)"

CXX_BIN="${CXX:-$(command -v g++ || true)}"
if [[ -z "${CXX_BIN}" ]]; then
  echo "ERROR: g++ compiler not found."
  echo "macOS: install Xcode command line tools and Homebrew gcc if needed."
  echo "Ubuntu/Debian: sudo apt install g++ libtbb-dev"
  echo "CentOS/RHEL: sudo yum install gcc-c++ tbb-devel"
  exit 1
fi

INCLUDES="-Isrc"
TBB_LIBS="-ltbb"
if [[ "${OS}" == "Darwin" ]]; then
  if [[ -d "/opt/homebrew/include/tbb" ]]; then
    INCLUDES="-Isrc -I/opt/homebrew/include"
    TBB_LIBS="-L/opt/homebrew/lib -ltbb"
  elif [[ -d "/usr/local/include/tbb" ]]; then
    INCLUDES="-Isrc -I/usr/local/include"
    TBB_LIBS="-L/usr/local/lib -ltbb"
  else
    echo "ERROR: TBB headers not found. Run: brew install tbb"
    exit 1
  fi
else
  if [[ -d "/usr/include/tbb" ]]; then
    INCLUDES="-Isrc -I/usr/include"
    TBB_LIBS="-L/usr/lib/x86_64-linux-gnu -ltbb"
  elif [[ -d "/usr/local/include/tbb" ]]; then
    INCLUDES="-Isrc -I/usr/local/include"
    TBB_LIBS="-L/usr/local/lib -ltbb"
  else
    echo "ERROR: TBB headers not found. Try: sudo apt install libtbb-dev or sudo yum install tbb-devel"
    exit 1
  fi
fi

build_target() {
  local target="$1"
  echo "Repository root : ${REPO_ROOT}"
  echo "Results root    : ${OUTPUT_ROOT}"
  echo "Machine label   : ${MACHINE_LABEL}"
  echo "Compiler        : ${CXX_BIN}"
  echo "Building target : ${target}"
  make "${target}" \
    CXX="${CXX_BIN}" \
    INCLUDES="${INCLUDES}" \
    TBB_LIBS="${TBB_LIBS}" \
    OUTPUT_ROOT="${OUTPUT_ROOT}" \
    MACHINE_LABEL="${MACHINE_LABEL}"
}

# Purpose:
#   Run the recovered-state front-speed diagnostic and verify the spreading
#   speed  c_inf = 2 sqrt( dF (1 - B) / (1 + theta^m) ),  B = beta theta^m.
# Paper item:
#   Main Figure 6 and Table 3.
# Driver:
#   drivers/main_fig6_front_speed.cpp
# Usage:
#   ./scripts/front_speed_exec.sh                 # beta=10, A0 = 0 2.25 5
#   ./scripts/front_speed_exec.sh 8               # secondary check at beta=8
#   ./scripts/front_speed_exec.sh 10 0 2.25 5     # explicit dose list
# Notes:
#   Panels (a) and (b) of Figure 6 use the A0=0 run with R_in = 1, so that
#   A = 0 and R = 1 hold exactly and the F equation is the recovered
#   Fisher-KPP equation from t = 0.  The positive doses feed only panel (c).
#
#   The run uses t_final = 160 with dt_out = 1.0, both longer and denser than
#   the Figure 5 run, because a front speed cannot be fitted from four
#   snapshots ending at t = 100.  Only a scalar level-set count is written,
#   so the total output is a few kilobytes.
#
#   The driver is invoked ONCE with the whole dose list rather than once per
#   dose: it opens its summary file in mode "w", so a per-dose loop would
#   leave only the last dose in summary_fig6.txt.

TARGET="main_fig6_front_speed"
BETA="${1:-10}"
shift || true
if [[ $# -gt 0 ]]; then
  A0_VALUES=("$@")
else
  A0_VALUES=("0" "2.25" "5")
fi

build_target "${TARGET}"
mkdir -p "${OUTPUT_ROOT}/logs"
LOG_FILE="${OUTPUT_ROOT}/logs/front_speed_beta_${BETA}_$(date +%Y%m%d_%H%M%S).log"

echo "Running front-speed driver with beta=${BETA}, A0 = ${A0_VALUES[*]}"
"${REPO_ROOT}/bin/${TARGET}" "${BETA}" "${A0_VALUES[@]}" 2>&1 | tee "${LOG_FILE}"

echo "Done. Outputs are under: ${OUTPUT_ROOT}front_speed/"
echo "Log file: ${LOG_FILE}"
echo "Next: python scripts/plot_fig6.py --beta ${BETA}"
