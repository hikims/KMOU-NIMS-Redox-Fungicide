#!/bin/bash
# ===========================================================
#  Auto Build & Run Script for TBB + C++17
#  Supports macOS, CentOS7, Ubuntu
#
#  Usage:
#     ./exec.sh           build, then run ./main for every A0 in list.txt
#  Requires:
#     list.txt            one initial dose A0 per line (see list.txt)
# ===========================================================

set -e

# -----------------------------------------------------------
#  Project root.
#  IMPORTANT: this must be the SAME path as the macro PATH in
#  main.cpp. Change both together if you move the project.
# -----------------------------------------------------------
PROJECT_ROOT="/Users/mkim/Desktop/Research/fungicide_jh"

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

SRC="main.cpp functions.cpp initial.cpp"
OUT="main"

echo "Building project..."
$CXX $STD $OPT $INCLUDE_PATH $LIB_PATH $SRC -ltbb -o $OUT
echo "Build completed successfully."

# -----------------------------------------------------------
#  Make sure list.txt exists
# -----------------------------------------------------------
if [[ ! -f list.txt ]]; then
    echo "ERROR: list.txt not found."
    echo "Create a file 'list.txt' with one A0 value per line."
    exit 1
fi

# -----------------------------------------------------------
#  Output directories.
#  main.cpp itself also creates data/, region/, summary/ under
#  PROJECT_ROOT; this just makes sure summary/ exists up front.
# -----------------------------------------------------------
SUMMARY_DIR="${PROJECT_ROOT}/summary"
mkdir -p "$SUMMARY_DIR"

echo "Running program for every A0 in list.txt ..."

while IFS= read -r line || [[ -n "$line" ]]; do
    # skip blank lines and comments
    [[ -z "$line" || "$line" =~ ^[[:space:]]*# ]] && continue

    A0=$(echo "$line" | tr -d '[:space:]')

    echo "Running: ./$OUT $A0"
    ./$OUT "$A0"

done < list.txt

echo "All runs finished."
echo "Summary file: ${SUMMARY_DIR}/summary_A0_vs_area.txt"
