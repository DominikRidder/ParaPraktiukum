#!/bin/bash

rm *out*
rm *err*
./compile.sh
sbatch run.ll

watch ls -l
