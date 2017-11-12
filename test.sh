#!/bin/bash

n=`grep CNF examples/prot.gro | grep NH | awk '{print $3}'`
c=`grep CNF examples/prot.gro | grep CT | awk '{print $3}'`
time ./new_g_nitrile_hbond -f examples/prot.gro \
    -s examples/v4.tpr \
    -a1 $c \
    -a2 $n \
    -select "not resname CNF and (same residue as within 0.5 of resname CNF and name NH)" \
    -o -oa -or -onwr 

#time ./new_g_nitrile_hbond -f examples/prot.xtc \
#    -s examples/v4.tpr \
#    -a1 $c \
#    -a2 $n \
#    -select "not resname CNF and (same residue as within 0.5 of resname CNF and name NH)" \
#    -o -oa -or -onwr ##-e 1000 


