#!/bin/bash
#Modify the charge of top file
bash ChangeTopCharge.sh COM.top COM_mod.top 3.0

#Run Minimization
pmemd.cuda -O -i MD_min1.in -p COM_mod.top -c COM.crd      -o Com_min1.out -r Com_min1.rst -ref COM.crd
pmemd.cuda -O -i MD_min2.in -p COM_mod.top -c Com_min1.rst -o Com_min2.out -r Com_min2.rst
rm -f mdinfo COM_mod.top
