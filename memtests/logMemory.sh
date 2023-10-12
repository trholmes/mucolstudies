PID=$1
while ps -p $PID &>/dev/null; do
    top -bn1 -p $PID | grep $PID | awk '{print $6}' >> memlog_$PID.txt
    sleep 5
done


