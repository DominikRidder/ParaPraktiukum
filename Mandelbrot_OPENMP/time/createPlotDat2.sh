#/bin/bash -X

cp speedtest1.dat speedtest2.dat

export OMP_NUM_THREADS=1
RET=$(../mandel -v -w 4096 -h 4096 -x -.59 -.54 -.58 -.53 -i 1024)
SERIALTIME=$(echo $RET | grep -o "run=[ ]*[0-9.]*" | grep -o "[0-9.]*")

# Update Serial time in plot
sed s/"seq = [0-9.]*"/"seq = $SERIALTIME"/ plottest2.gpl > tmp.txt
mv tmp.txt plottest2.gpl
