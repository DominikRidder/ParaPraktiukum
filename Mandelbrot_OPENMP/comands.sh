#!/bin/bash

g++ -fopenmp -o mandel mandelseq2.c ppmwrite.c  

export OMP_SCHEDULE="static,2"
export OMP_NUM_THREADS=2
./mandel -v -w 4096 -h 4096 -x -.59 -.54 -.58 -.53 -i 1024 # zum testen von mandelcol.ppm
convert mandelcol.ppm mandelcol.png

# Timing
timingOn=1
if [[ $timingOn == 1 ]]; then
	cd ./time

	./createPlotDat1.sh
	./createPlotDat2.sh
	./createPlotDat3.sh

	gnuplot plottest1.gpl
	gnuplot plottest2.gpl
	gnuplot plottest3.gpl

	cd ..
fi

