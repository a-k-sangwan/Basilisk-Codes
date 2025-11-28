#!/usr/bin/env bash
# Simple helper script to compile and run the Basilisk simulation.
# Usage:
#   ./run.sh
#
# Requirements:
#   - Basilisk C installed and `qcc` available in PATH
#   - Python 3 with numpy, scipy, sympy, matplotlib
#   - Files: single.c, genData.py, output_vtu.h all in the same directory

EXEC_NAME="run"

echo "Compiling single.c with qcc..."
qcc main.c -lm -O2  -o "$EXEC_NAME"

echo "Running simulation..."
echo " "
./"$EXEC_NAME"
