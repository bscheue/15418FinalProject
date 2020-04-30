#!/usr/bin/env bash

echo "small seq"
./crun-seq -m 15 -t
echo "small omp"
./crun-omp -m 15 -t
echo "larger seq"
./crun-seq -m 70 -t
echo "larger omp"
./crun-omp -m 70 -t

