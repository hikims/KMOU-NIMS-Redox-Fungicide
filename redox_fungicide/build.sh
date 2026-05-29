#!/usr/bin/env bash
set -euo pipefail

: "${OUTPUT_ROOT:=./results/}"
: "${MACHINE_LABEL:=local}"

echo "Building all drivers"
echo "OUTPUT_ROOT=${OUTPUT_ROOT}"
echo "MACHINE_LABEL=${MACHINE_LABEL}"
make all OUTPUT_ROOT="${OUTPUT_ROOT}" MACHINE_LABEL="${MACHINE_LABEL}"
