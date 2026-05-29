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
#   Run the convergence / invariant-region diagnostic driver.
# Paper item:
#   Supplementary Tables S1--S2 and numerical invariant-region diagnostics.
# Driver:
#   drivers/main_table_S1_S2_convergence.cpp
# Usage:
#   ./scripts/convergence_exec.sh          # default A0=0.5
#   ./scripts/convergence_exec.sh 0.5
# Notes:
#   Spatial/temporal refinement still requires editing src/parameters.h and rebuilding,
#   as described in the supplement and in src/parameters.h.

TARGET="main_table_S1_S2_convergence"
A0="${1:-0.5}"

build_target "${TARGET}"
mkdir -p "${OUTPUT_ROOT}/logs"
LOG_FILE="${OUTPUT_ROOT}/logs/convergence_A0_${A0}_$(date +%Y%m%d_%H%M%S).log"

echo "Running convergence driver with A0=${A0}"
"${REPO_ROOT}/bin/${TARGET}" "${A0}" 2>&1 | tee "${LOG_FILE}"

echo "Done. Outputs are under: ${OUTPUT_ROOT}"
echo "Log file: ${LOG_FILE}"
