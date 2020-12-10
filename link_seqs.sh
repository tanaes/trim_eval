#!/bin/bash

indir=$1
shift # short for shift 1
for i
do
    f=${indir}/${i}*_R1_*.fastq.gz
    echo $f
    ln -s $f ${i}.R1.fastq.gz

    r=${indir}/${i}*_R2_*.fastq.gz
    echo $r
    ln -s $f ${i}.R2.fastq.gz
done