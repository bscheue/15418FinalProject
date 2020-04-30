#!/usr/bin/env bash

mkdir -p test/seq test/omp
# need to clean out the results of the previous tests
# since the file is opened to append
rm test/seq/* test/omp/*

# small test with default seed
./crun-seq -m 15 -o "test/seq/small_def.txt"
./crun-omp -m 15 -o "test/omp/small_def.txt"

# small test with default seed
./crun-seq -m 15 -o "test/seq/small_418.txt" -s 418
./crun-omp -m 15 -o "test/omp/small_418.txt" -s 418

# larger test with default seed
./crun-seq -m 70 -o "test/seq/larger_def.txt"
./crun-omp -m 70 -o "test/omp/larger_def.txt"

diff -q test/seq test/omp
