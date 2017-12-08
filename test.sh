#!/bin/bash

source /usr/local/gromacs/bin/GMXRC
#time g_nitrile_hbond -f examples/prot.gro -h

n=`grep CNF examples/prot.gro | grep NH | awk '{print $3}'`
c=`grep CNF examples/prot.gro | grep CT | awk '{print $3}'`

time g_nitrile_hbond -f examples/prot.gro \
    -s examples/v4.tpr \
    -a1 $c \
    -a2 $n \
    -select "not resname CNF and (same residue as within 0.5 of resname CNF and name NH)" \
    -o -oa -or -onwr -op 

time ./g_nitrile_hbond -f examples/prot.xtc \
    -s examples/v4.tpr \
    -a1 $c \
    -a2 $n \
    -select "not resname CNF and (same residue as within 0.5 of resname CNF and name NH)" \
    -o -oa -or -onwr -op #-e 1000


