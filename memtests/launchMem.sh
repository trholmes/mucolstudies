python test.py &
PID=$!
source logMemory.sh $PID
python plotMemory.py $!
echo Finished monitoring process $PID

