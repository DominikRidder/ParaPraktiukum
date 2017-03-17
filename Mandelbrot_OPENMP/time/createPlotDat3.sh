#/bin/bash -X

OPENMP=(procs method1 method2 method3 method4)
MASTER=(procs master)
NEXTLINE=(procs method1 method2 method3 method4 master)

rm speedtest3.dat
touch speedtest3.dat

# Collect OpenMP

COUNTER=0
for i in `cat speedtest1.dat`
do	
	OPENMP[$COUNTER]=$i

	COUNTER=$((COUNTER+1))
done

# Collect Master Slave MPI

COUNTER=0
for i in `cat MasterSlave.dat`
do	
	MASTER[$COUNTER]=$i

	COUNTER=$((COUNTER+1))
done

# Merge

for i in `seq 1 25` # headline + 1-24 threads
do
	for j in `seq 0 4` # procs + 4 openmp methods
	do
		NEXTLINE[$j]=${OPENMP[$(((i-1)*5+j))]}
	done
	NEXTLINE[5]=${MASTER[$(((i-1)*2+1))]}
	echo ${NEXTLINE[@]} >> speedtest3.dat
done
