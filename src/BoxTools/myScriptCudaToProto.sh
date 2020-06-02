#!/bin/bash

echo " Script $0 is launched "

FILE=$1

echo " Apply to $FILE"

LIST="PointerAttributes PointerGetAttributes( MemcpyFromSymbol Stream_t Error Success DeviceProp Malloc( Free( HostAlloc( FreeHost( Memcpy( MemcpyAsync( MemcpyDeviceToDevice MemcpyDeviceToHost MemcpyHostToDevice GetLastError GetErrorString StreamCreate StreamDestroy DeviceSynchronize StreamSynchronize SetDevice GetDeviceProperties"

for VAR in $LIST
do
	FROM=cuda$VAR
	TO=proto$VAR 
	count=`grep -c "$FROM" $FILE`
	if [ "$count" -gt 0 ]
	then
		echo " Convert $FROM to $TO -> Ocurrence(s) $count"
		sed -i "s/$FROM/$TO/g" $FILE
		count2=`grep -c "$TO" $FILE`
		if [ "$count" -eq "$count2" ]
		then
			echo " Check Success"
		else
			echo " Check Failled"
		fi
	fi

done
