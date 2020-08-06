#!/bin/bash
#Setting
PDBName="1eby"
NetCharge=0

#Calculate AM1-BCC Charge Without Opt
antechamber -i $PDBName"_ligand.mol2" -fi mol2 -o $PDBName"_ligand_noOpt_BCC_gaff2.mol2" -fo mol2 -c bcc -nc $NetCharge -rn "LIG" -ek "qm_theory='AM1', maxcyc=0, scfconv=1.d-8, ndiis_attempts=700," -at gaff2 -pl 25 -pf y -dr no 1>/dev/null
rm -f sqm.in sqm.out sqm.pdb

#Get FrcMod File
parmchk2 -i $PDBName"_ligand_noOpt_BCC_gaff2.mol2" -f mol2 -o $PDBName"_ligand_gaff2.frcmod" -s gaff2 -a N 1>/dev/null
