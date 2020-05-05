#!/usr/bin/env bash

echo "larger omp ref" >&2
./crun-omp-ref -m 125 -b -t 8
echo "larger omp" >&2
./crun-omp -m 125 -b -t 8

