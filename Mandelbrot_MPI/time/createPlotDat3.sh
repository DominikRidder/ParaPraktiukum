#/bin/bash -X

declare -A ALLTIMES ORDER

rm speedtest3.dat
touch speedtest3.dat
echo "process-ID ROW BLOCK MASTER" >> speedtest3.dat

INDEX=0

INDEXROW=0
INDEXBLOCK=0
INDEXMASTER=0
PROCS=48

for t in `seq 0 2`
do
    RET=$(srun --ntasks=$PROCS ../mandelseq -v -w 4096 -h 4096 -x -.59 -.54 -.58 -.53 -i 1024 -t $t)

    INDEX=0
    for i in $(echo $RET | grep -o "PE[0-9]*" | grep -o "[0-9]*") # ORDERING
    do
	ORDER[$INDEX]=$i
	INDEX=$((INDEX+1))
    done

    INDEX=0
    for i in $(echo $RET | grep -o "mpi=[ ]*[0-9.]*" | grep -o "[0-9.]*") # MPI
    do
	NEXT=$(echo "${ORDER[$INDEX]} + $t * 48" | bc -l)
	ALLTIMES[$NEXT]=$i
	INDEX=$((INDEX+1))
    done

    INDEX=0
    for i in $(echo $RET | grep -o "calc=[ ]*[0-9.]*" | grep -o "[0-9.]*") # CALC
    do
	NEXT=$(echo "${ORDER[$INDEX]} + $t * 48" | bc -l)
	ALLTIMES[$NEXT]=$(echo "${ALLTIMES[$NEXT]} + $i" | bc -l)
	INDEX=$((INDEX+1))
    done
done

for procID in `seq 0 $((PROCS-1))`
do
    INDEXROW=$procID
    INDEXBLOCK=$((procID+$PROCS))
    INDEXMASTER=$((procID+$PROCS+$PROCS))
    echo "$((procID)) ${ALLTIMES[$INDEXROW]} ${ALLTIMES[$INDEXBLOCK]} ${ALLTIMES[$INDEXMASTER]}" >> speedtest3.dat
done
cat speedtest3.dat
