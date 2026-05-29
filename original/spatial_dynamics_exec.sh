#!/bin/bash
# ===========================================
#  Auto Build & Run Script for Figure 5 rerun
#  A0 = 0, 2.25, 5 only
#  Supports macOS, CentOS7, Ubuntu
# ===========================================

set -e

OS=$(uname -s)
echo "Detected OS: $OS"

CXX=$(which g++)
if [[ -z "$CXX" ]]; then
    echo "g++ compiler not found!"
    echo "Ubuntu/Debian: sudo apt install g++"
    echo "CentOS/RHEL:   sudo yum install gcc-c++"
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
        echo "Using Ubuntu TBB path (/usr/include, /usr/lib/x86_64-linux-gnu)"
    elif [[ -d "/usr/local/include/tbb" ]]; then
        INCLUDE_PATH="-I/usr/local/include"
        LIB_PATH="-L/usr/local/lib"
        echo "Using /usr/local TBB path"
    else
        echo "TBB not found. Try: sudo apt install libtbb-dev"
        exit 1
    fi
fi

# ===========================================
# Source and executable
# ===========================================
SRC="main_fig5.cpp functions.cpp initial.cpp"
OUT="main_fig5"

echo "Building project..."
$CXX $STD $OPT $INCLUDE_PATH $LIB_PATH $SRC -ltbb -o $OUT

echo "Build completed successfully."

# ===========================================
# Figure 5 rerun settings
# ===========================================
BASE="/Users/mkim/Desktop/Research/fungicide_jh"
RUN_NAME="fig5_rerun_A0_0_2p25_5"

RUN_DIR="$BASE/$RUN_NAME"
SUMMARY_DIR="$RUN_DIR/summary"

mkdir -p "$RUN_DIR"
mkdir -p "$SUMMARY_DIR"

SUMMARY_FILE="$SUMMARY_DIR/summary_fig5_rerun.txt"

echo "# A0 total_F_t20 area_t20 total_F_t50 area_t50 total_F_t100 area_t100" > "$SUMMARY_FILE"

# Only these three doses
A0_VALUES=("0" "2.25" "5")

echo "==========================================="
echo "Running Figure 5 simulations"
echo "A0 values: ${A0_VALUES[*]}"
echo "Output root: $RUN_DIR"
echo "==========================================="

for A0 in "${A0_VALUES[@]}"; do
    echo ""
    echo "-------------------------------------------"
    echo "Running: ./$OUT $A0"
    echo "-------------------------------------------"

    ./$OUT "$A0"
done

echo ""
echo "All Figure 5 reruns finished."
echo "Check output folders under:"
echo "$RUN_DIR"