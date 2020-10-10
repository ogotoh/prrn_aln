#!/bin/sh

prrn5 -pi pas/ce13a17.fa
aln -yl2 -L -pi nas/CET10B9 pas/ce13a.msa
aln -s pas Multi_A Multi_B
prrn5 -s pas Multi_A Multi_B
prrn5 -s pas -U Multi_A Multi_B
