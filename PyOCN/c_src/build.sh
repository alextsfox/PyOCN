#!/bin/bash

set -e

# Get the directory of this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Move to the script directory
cd "$SCRIPT_DIR"

# Compile into position-independent object files
# gcc -fPIC -c ocn.c flowgrid.c status.c rng.c
gcc -fPIC -O3 -flto -Wall -pedantic -std=c99 \
    -c ocn.c flowgrid.c status.c rng.c

# Link into shared library
# gcc -shared -o libocn.so ocn.o flowgrid.o status.o rng.o
gcc -shared -O3 -flto -Wall -pedantic \
    -o libocn_dev.so ocn.o flowgrid.o status.o rng.o

# Move the shared library to the Python package root
mv libocn_dev.so ../libocn_dev.so

# Clean up object files
rm -f ocn.o flowgrid.o status.o rng.o

echo "Built libocn_dev.so in $(cd .. && pwd)"