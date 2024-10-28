#!/bin/bash
mkdir -p $4
bowtie2 -x $3 -1 $1 -2 $2 --un-conc-gz $4 -p 32 -S $4/aligned_rRNA.sam
rm $4/aligned_rRNA.sam