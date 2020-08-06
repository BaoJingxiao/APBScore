#!/bin/bash
#Settings
AmberHome="XXX"
Fortran="gfortran"
Flags="-march=native -O3 -s"

#Make
if [ ! -d ../bin ];then
    mkdir ../bin
fi
$Fortran APBScore.f95 $Flags -I $AmberHome/include $AmberHome/lib/libnetcdff.a $AmberHome/lib/libnetcdf.a -o ../bin/APBScore
rm -f *.mod
