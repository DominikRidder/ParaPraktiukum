#/bin/bash -X

COUNT=$#
if [[ "$COUNT" == "0" ]];
then
    COUNT=1
fi
SUM=0

while (( "$#" )); do
    ADD=$1
    SUM=$(echo "$SUM + $ADD" | bc -l)
    shift
done

SUM=$(echo "$SUM / $COUNT" | bc -l | grep -o "........" | head -1)

TEST=$(echo $SUM | grep -o "." | head -1)
if [[ "$TEST" == "." ]];
then
    SUM=$(echo "0$SUM" | grep -o "........" | head -1)
fi

echo $SUM
