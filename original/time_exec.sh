#!/bin/bash
# ===========================================
#  All-in-one Build & Sweep Script
#  Compiles main4.cpp once, then runs over multiple beta values.
#
#  Edit BETA_LIST and A0 below to change the sweep.
#
#  Usage:
#    ./run.sh                 # uses defaults below
#    nohup caffeinate -i ./run.sh > sweep.log 2>&1 &   # background, prevent sleep
# ===========================================

set -e

# =============================================================
# >>>>> EDIT THESE TO CHANGE THE SWEEP <<<<<
# =============================================================
A0="0.5"
BETA_LIST="5 8 10 11 12 15 20"
# =============================================================

# Output paths -- must match PATH macro in main4.cpp
ROOT_DIR="/Users/mkim/Desktop/Research/fungicide_jh"
TIME_SERIES_DIR="$ROOT_DIR/time_series"
THEOREM_DIR="$TIME_SERIES_DIR/theorem_1_5_II"
GLOBAL_SUMMARY="$THEOREM_DIR/global_summary.txt"
LOG_DIR="$ROOT_DIR/sweep_logs/sweep_$(date +%Y%m%d_%H%M%S)"

# -------------------------------------------------------------
# Compiler detection
# -------------------------------------------------------------
OS=$(uname -s)
echo "Detected OS: $OS"

CXX=$(which g++)
if [[ -z "$CXX" ]]; then
    echo "g++ compiler not found!"
    echo "  Ubuntu/Debian: sudo apt install g++"
    echo "  CentOS/RHEL  : sudo yum install gcc-c++"
    exit 1
fi

STD="-std=c++17"
OPT="-O3"
INCLUDE_PATH=""
LIB_PATH=""

if [[ "$OS" == "Darwin" ]]; then
    if [[ -d "/opt/homebrew/include" ]]; then
        INCLUDE_PATH="-I/opt/homebrew/include"
        LIB_PATH="-L/opt/homebrew/lib"
        echo "Using Homebrew TBB path"
    else
        echo "Homebrew TBB not found. Run: brew install tbb"
        exit 1
    fi
else
    if [[ -d "/usr/include/tbb" ]]; then
        INCLUDE_PATH="-I/usr/include"
        LIB_PATH="-L/usr/lib/x86_64-linux-gnu"
        echo "Using Ubuntu TBB path"
    elif [[ -d "/usr/local/include/tbb" ]]; then
        INCLUDE_PATH="-I/usr/local/include"
        LIB_PATH="-L/usr/local/lib"
        echo "Using /usr/local TBB path"
    else
        echo "TBB not found. Try: sudo apt install libtbb-dev"
        exit 1
    fi
fi

# -------------------------------------------------------------
# Build (only if needed)
# -------------------------------------------------------------
SRC="main4.cpp functions.cpp initial.cpp"
OUT="main4"

NEEDS_BUILD=1
if [[ -x "./$OUT" ]]; then
    NEEDS_BUILD=0
    for f in $SRC parameters.h functions.hpp initial.hpp; do
        if [[ -f "$f" && "$f" -nt "./$OUT" ]]; then
            NEEDS_BUILD=1
            break
        fi
    done
fi

if [[ $NEEDS_BUILD -eq 1 ]]; then
    echo "Building project..."
    $CXX $STD $OPT $INCLUDE_PATH $LIB_PATH $SRC -ltbb -o $OUT
    echo "Build completed."
else
    echo "Build is up-to-date, skipping recompile."
fi

# -------------------------------------------------------------
# Prepare output directories and clean previous summary
# -------------------------------------------------------------
mkdir -p "$THEOREM_DIR"
mkdir -p "$LOG_DIR"

if [[ -f "$GLOBAL_SUMMARY" ]]; then
    echo "Removing previous global_summary.txt"
    rm -f "$GLOBAL_SUMMARY"
fi

# -------------------------------------------------------------
# Banner
# -------------------------------------------------------------
TOTAL=$(echo $BETA_LIST | wc -w | tr -d ' ')

echo ""
echo "=================================================="
echo "  Theorem 1.5 II - beta sweep"
echo "=================================================="
echo "  A0          : $A0"
echo "  beta values : $BETA_LIST  ($TOTAL runs)"
echo "  Output root : $ROOT_DIR"
echo "  Logs        : $LOG_DIR"
echo "  Started     : $(date)"
echo "=================================================="
echo ""

START_TIME=$(date +%s)

# -------------------------------------------------------------
# Sweep loop
# -------------------------------------------------------------
COUNT=0
for BETA in $BETA_LIST; do
    COUNT=$((COUNT + 1))
    LOG_FILE="$LOG_DIR/run_beta_${BETA}.log"

    echo "--------------------------------------------------"
    echo "[$COUNT/$TOTAL]  A0=$A0  beta=$BETA"
    echo "  log     : $LOG_FILE"
    echo "  started : $(date '+%H:%M:%S')"
    echo "--------------------------------------------------"

    RUN_START=$(date +%s)

    # Run main4 directly (no exec.sh needed!)
    ./$OUT "$A0" "$BETA" > "$LOG_FILE" 2>&1

    RUN_END=$(date +%s)
    RUN_ELAPSED=$((RUN_END - RUN_START))

    echo "  finished: $(date '+%H:%M:%S')  (took ${RUN_ELAPSED}s)"
    echo "  tail of log:"
    tail -n 3 "$LOG_FILE" | sed 's/^/    /'
    echo ""
done

# -------------------------------------------------------------
# Summary
# -------------------------------------------------------------
END_TIME=$(date +%s)
TOTAL_ELAPSED=$((END_TIME - START_TIME))
HOURS=$((TOTAL_ELAPSED / 3600))
MINUTES=$(( (TOTAL_ELAPSED % 3600) / 60 ))
SECS=$((TOTAL_ELAPSED % 60))

echo "=================================================="
echo "  Sweep complete"
echo "=================================================="
echo "  Total runs    : $TOTAL"
echo "  Total elapsed : ${HOURS}h ${MINUTES}m ${SECS}s"
echo "  Finished      : $(date)"
echo "  Summary       : $GLOBAL_SUMMARY"
echo "=================================================="

if [[ -f "$GLOBAL_SUMMARY" ]]; then
    echo ""
    echo "Contents of global_summary.txt:"
    echo "--------------------------------------------------"
    cat "$GLOBAL_SUMMARY"
    echo "--------------------------------------------------"
fi