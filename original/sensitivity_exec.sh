#!/bin/bash
# ===========================================
#  Auto Build & Run Script for TBB + C++17
#  Sensitivity analysis version
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

SRC="main3.cpp functions.cpp initial.cpp"
OUT="main3"

echo "Building project..."
$CXX $STD $OPT $INCLUDE_PATH $LIB_PATH $SRC -ltbb -o $OUT

echo "Build completed successfully."

BASE="/Users/mkim/Desktop/Research/fungicide_jh"
LIST_FILE="$BASE/code/list.txt"

if [[ ! -f "$LIST_FILE" ]]; then
    echo "list.txt not found at: $LIST_FILE"
    exit 1
fi

mkdir -p "$BASE/sensitivity/data"
mkdir -p "$BASE/sensitivity/region"
mkdir -p "$BASE/sensitivity/summary"

echo "Running sensitivity analysis..."
./$OUT

echo "All sensitivity runs finished."
echo "Results saved in:"
echo "$BASE/sensitivity/"