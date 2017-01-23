#/bin/bash -X

RUNTIME=(method1 method2 method3)

rm speedtest1.dat
touch speedtest1.dat
echo "#processes ROW BLOCK MASTER" >> speedtest1.dat

for procs in `seq 1 48`
do
    for t in `seq 0 2`
    do
	RET=$(srun --ntasks=$procs ../mandelseq -v -w 4096 -h 4096 -x -.59 -.54 -.58 -.53 -i 1024 -t $t)
	ALLTIMES=$(echo $RET | grep -o "run=[ ]*[0-9.]*" | grep -o "[0-9.]*")
	if [[ "X$ALLTIMES" == "X" ]];
	then
	    RUNTIME[$t]="-       "
	else
	    RUNTIME[$t]=$(./mean.sh $ALLTIMES)
	fi
    done
    echo "$procs ${RUNTIME[*]}"
    echo "$procs ${RUNTIME[*]}" >> speedtest1.dat
done
