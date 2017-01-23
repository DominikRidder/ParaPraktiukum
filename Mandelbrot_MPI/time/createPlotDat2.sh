#/bin/bash -X

cp speedtest1.dat speedtest2.dat

RET=$(srun --ntasks=1 ../mandelseq -t 3 -v)
SERIALTIME=$(echo $RET | grep -o "run=[ ]*[0-9.]*" | grep -o "[0-9.]*")

# Update Serial time in plot
sed s/"seq = [0-9.]*"/"seq = $SERIALTIME"/ plottest2.gpl > tmp.txt
mv tmp.txt plottest2.gpl
