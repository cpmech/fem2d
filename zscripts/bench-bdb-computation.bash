#!/bin/bash

set -e

# compile optimized code
touch benchmarks/bdb-computation/bdb_classical.cpp
touch benchmarks/bdb-computation/bdb_expanded_full.cpp
touch benchmarks/bdb-computation/bdb_expanded.cpp
bash all.bash ON

# change to build dir
cd /tmp/build-fem2d/benchmarks/bdb-computation

# run benchmarks
echo "### bdb_expanded #######################"
./bdb_expanded

echo "### bdb_expanded_full ##################"
./bdb_expanded_full

echo "### bdb_classical ######################"
./bdb_classical
