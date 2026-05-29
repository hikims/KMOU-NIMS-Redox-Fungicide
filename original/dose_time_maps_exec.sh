#!/bin/bash
# ===========================================
# Build & run MATCOM Tier-B B2-only driver
# Machine: macstudio
# Output: /Users/mkim/Desktop/Research/fungicide_jh/TierB/
# ===========================================

set -euo pipefail

BASE_DIR="/Users/mkim/Desktop/Research/fungicide_jh/"
LABEL="macstudio"
OS=$(uname -s)

echo "Detected OS: $OS"
echo "Base directory: $BASE_DIR"
echo "Machine label: $LABEL"

CXX=$(which g++ || true)
if [[ -z "$CXX" ]]; then
    echo "g++ compiler not found!"
    echo "macOS: install Xcode command line tools and Homebrew gcc if needed"
    echo "Ubuntu/Debian: sudo apt install g++ libtbb-dev"
    echo "CentOS/RHEL: sudo yum install gcc-c++ tbb-devel"
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
        echo "Using Homebrew TBB path (/opt/homebrew)"
    elif [[ -d "/usr/local/include" ]]; then
        INCLUDE_PATH="-I/usr/local/include"
        LIB_PATH="-L/usr/local/lib"
        echo "Using Homebrew TBB path (/usr/local)"
    else
        echo "Homebrew TBB not found. Run: brew install tbb"
        exit 1
    fi
else
    if [[ -d "/usr/include/tbb" ]]; then
        INCLUDE_PATH="-I/usr/include"
        LIB_PATH="-L/usr/lib/x86_64-linux-gnu"
        echo "Using Ubuntu/Debian TBB path"
    elif [[ -d "/usr/local/include/tbb" ]]; then
        INCLUDE_PATH="-I/usr/local/include"
        LIB_PATH="-L/usr/local/lib"
        echo "Using /usr/local TBB path"
    else
        echo "TBB not found. Try: sudo apt install libtbb-dev or sudo yum install tbb-devel"
        exit 1
    fi
fi

SRC="main.cpp functions.cpp initial.cpp"
OUT="b2_main_${LABEL}"

mkdir -p "${BASE_DIR}/TierB/summary" "${BASE_DIR}/TierB/logs"
cp -f b2_cases.txt "${BASE_DIR}/TierB/b2_cases.txt"

echo "Building B2-only driver..."
$CXX $STD $OPT $INCLUDE_PATH $LIB_PATH $SRC -ltbb -o "$OUT"

echo "Build completed successfully."
echo "Running B2-only computations..."

LOG_FILE="${BASE_DIR}/TierB/logs/B2_${LABEL}_$(date +%Y%m%d_%H%M%S).log"
if [[ "$OS" == "Darwin" ]]; then
    /usr/bin/time -p ./"$OUT" 2>&1 | tee "$LOG_FILE"
else
    /usr/bin/time -v ./"$OUT" 2>&1 | tee "$LOG_FILE"
fi

echo "B2 computations finished."
echo "Outputs are in: ${BASE_DIR}/TierB/summary"
echo "Log file: $LOG_FILE"
