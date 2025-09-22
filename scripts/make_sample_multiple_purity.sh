PID=$$
PURITY=$4
for i in $(seq 1 $1)
do	
	echo "Iteration = $i"
	PID=$((PID+1))
	bash make_sample.sh $2 $3 $PURITY pid=$PID
done
