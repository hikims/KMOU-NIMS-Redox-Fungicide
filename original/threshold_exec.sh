#!/bin/bash
# ===========================================
#  Auto Build & Run Script for TBB + C++17
#  Threshold / Sensitivity analysis version
# ===========================================

set -e

OS=$(uname -s)
echo "Detected OS: $OS"

CXX=$(which g++)
if [[ -z "$CXX" ]]; then
    echo "g++ compiler not found!"
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
        echo "Run: brew install tbb"
        exit 1
    fi
else
    if [[ -d "/usr/include/tbb" ]]; then
        INCLUDE_PATH="-I/usr/include"
        LIB_PATH="-L/usr/lib/x86_64-linux-gnu"
    else
        echo "Install TBB: sudo apt install libtbb-dev"
        exit 1
    fi
fi

SRC="main3.cpp functions.cpp initial.cpp"
OUT="main3"

echo "Building project..."
$CXX $STD $OPT $INCLUDE_PATH $LIB_PATH $SRC -ltbb -o $OUT

echo "Build completed successfully."

BASE="/Users/mkim/Desktop/Research/fungicide_jh"
LIST_FILE="$BASE/code/threshold_list.txt"

if [[ ! -f "$LIST_FILE" ]]; then
    echo "list.txt not found at: $LIST_FILE"
    exit 1
fi

# threshold 전용 폴더
mkdir -p "$BASE/threshold_parameter/data"
mkdir -p "$BASE/threshold_parameter/region"
mkdir -p "$BASE/threshold_parameter/summary"

echo "Running threshold analysis..."
./$OUT

echo "======================================"
echo " All runs finished"
echo "======================================"
echo "Results saved in:"
echo "$BASE/threshold_parameter/"