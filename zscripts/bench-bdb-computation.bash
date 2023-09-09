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
./bmark_bdb "expanded"
./bmark_bdb "expanded_full"
./bmark_bdb "classical"
