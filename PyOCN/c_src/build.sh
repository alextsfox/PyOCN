#!/bin/bash

set -e

# Get the directory of this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Move to the script directory
cd "$SCRIPT_DIR"

# Compile into position-independent object files
gcc -fPIC -c ocn.c streamgraph.c status.c rng.c
# gcc -fPIC -O3 -flto -c ocn.c streamgraph.c status.c rng.c

# Link into shared library
gcc -shared -o libocn.so ocn.o streamgraph.o status.o rng.o
# gcc -shared -O3 -flto -o libocn.so ocn.o streamgraph.o status.o rng.o

# Move the shared library to the Python package root
mv libocn.so ../libocn.so

# Clean up object files
rm -f ocn.o streamgraph.o status.o rng.o

echo "Built libocn.so in $(cd .. && pwd)"