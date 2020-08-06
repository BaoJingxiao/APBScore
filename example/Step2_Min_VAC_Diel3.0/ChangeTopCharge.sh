#!/bin/bash
InTop=$1
OutTop=$2
Diel=$3
if [ `echo $Diel" > 0" | bc` -eq 1 ];then
    Temp=`echo "scale=10; sqrt( 1 / "$Diel" )" | bc`
else
    Temp=0
fi
LineID1=`cat $InTop | grep -n "%FLAG CHARGE" | awk -F ":" '{print($1+2)}'`
LineID2=`cat $InTop | grep -n "%FLAG ATOMIC_NUMBER" | awk -F ":" '{print($1-1)}'`
cat $InTop | head -n $[$LineID1-1] > $OutTop
cat $InTop | head -n $LineID2 | tail -n $[$LineID2-$LineID1+1] | awk -v C=$Temp -F " " '{for(n=1;n<=NF;n++){if(n!=NF){printf("%16.8E",$n*C)}else{printf("%16.8E\n",$n*C)}}}' >> $OutTop
cat $InTop | tail -n +$[$LineID2+1] >> $OutTop
