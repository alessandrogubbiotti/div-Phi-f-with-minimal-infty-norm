#!/bin/bash

# Name of output executable
OUTPUT="poisson_mc"

# Source file
SRC="poisson_mc.c"

# Compilation command
gcc -O2 "$SRC" -o "$OUTPUT" -lfftw3 -lm -lgsl -lgslcblas

# Check compilation result
if [ $? -eq 0 ]; then
    echo "Compilation successful. Executable: $OUTPUT"
else
    echo "Compilation failed."
fi

