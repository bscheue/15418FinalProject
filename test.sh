#!/usr/bin/env bash

rm -rf code
# move asst3 code from afs to head node
cp -r ~/AFS/private/15418/finalproject/code .
cd code

# as per @406
module load gcc-4.9.2
module load binutils-2.26

make clean
make all

./submitjob.py

while [ ! -f ~/code/*.sh.o* ]; do sleep 1; done
cat ~/code/*.sh.o* | tail -n8
