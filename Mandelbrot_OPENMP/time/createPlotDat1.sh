#/bin/bash -X

RUNTIME=(method1 method2 method3 method4)

rm speedtest1.dat
touch speedtest1.dat
echo "#processes static,1 dynamic guided auto" >> speedtest1.dat

for procs in `seq 1 24`
do
    COUNTER=0
    for t in "static,1" "dynamic" "guided" "auto"
    do
	export OMP_SCHEDULE=$t
	export OMP_NUM_THREADS=$procs	
	RET=$(../mandel -v -w 4096 -h 4096 -x -.59 -.54 -.58 -.53 -i 1024)
	ALLTIMES=$(echo $RET | grep -o "run=[ ]*[0-9.]*" | grep -o "[0-9.]*")
	if [[ "X$ALLTIMES" == "X" ]];
	then
	    RUNTIME[$COUNTER]="-       "
	else
	    RUNTIME[$COUNTER]=$(./mean.sh $ALLTIMES)
	fi
	COUNTER=$((COUNTER+1))
    done
    echo "$procs ${RUNTIME[*]}"
    echo "$procs ${RUNTIME[*]}" >> speedtest1.dat
done
