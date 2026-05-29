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
#   Run beta-sweep time-series simulations for phase-transition verification.
# Paper item:
#   Main Figures 3--4.
# Driver:
#   drivers/main_fig3_fig4_verification.cpp
# Usage:
#   ./scripts/time_exec.sh
#   A0=0.5 BETA_LIST="5 8 10 11 12 15 20" ./scripts/time_exec.sh

TARGET="main_fig3_fig4_verification"
A0="${A0:-0.5}"
BETA_LIST="${BETA_LIST:-5 8 10 11 12 15 20}"

build_target "${TARGET}"
LOG_DIR="${OUTPUT_ROOT}/sweep_logs/sweep_$(date +%Y%m%d_%H%M%S)"
mkdir -p "${LOG_DIR}"

GLOBAL_SUMMARY="${OUTPUT_ROOT}/time_series/theorem_1_5_II/global_summary.txt"
if [[ -f "${GLOBAL_SUMMARY}" ]]; then
  echo "Removing previous global summary: ${GLOBAL_SUMMARY}"
  rm -f "${GLOBAL_SUMMARY}"
fi

TOTAL=$(echo ${BETA_LIST} | wc -w | tr -d ' ')
echo "Running beta sweep for Figures 3--4"
echo "A0=${A0}"
echo "beta values=${BETA_LIST}"
echo "logs=${LOG_DIR}"

COUNT=0
for BETA in ${BETA_LIST}; do
  COUNT=$((COUNT + 1))
  LOG_FILE="${LOG_DIR}/run_beta_${BETA}.log"
  echo "[$COUNT/$TOTAL] A0=${A0}, beta=${BETA}"
  "${REPO_ROOT}/bin/${TARGET}" "${A0}" "${BETA}" > "${LOG_FILE}" 2>&1
  tail -n 3 "${LOG_FILE}" | sed 's/^/  /'
done

echo "Done. Time-series outputs are under: ${OUTPUT_ROOT}/time_series"
echo "Logs are under: ${LOG_DIR}"
if [[ -f "${GLOBAL_SUMMARY}" ]]; then
  echo "Global summary: ${GLOBAL_SUMMARY}"
fi
