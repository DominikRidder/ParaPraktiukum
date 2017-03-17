#!/bin/bash

FIELDS=(Prozesses MasterSlave)
COUNTER=0
rm MasterSlave.dat
touch MasterSlave.dat
for i in `cat MasterSlaveRaw.dat`; do
	if [[ $((COUNTER%4)) == 0 ]]; then
		FIELDS[0]=$i
	fi
	if [[ $((COUNTER%4)) == 3 ]]; then
		FIELDS[1]=$i
		echo ${FIELDS[@]} >> MasterSlave.dat
	fi
	COUNTER=$((COUNTER+1))
done
