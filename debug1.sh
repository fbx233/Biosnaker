#!/bin/bash

# 定义变量
rRNA_index="/datapool/home/fengbx/data1/ref_genome/index/rRNA/Homo_sapiens.rRNA"
clean_R1="/datapool/home/fengbx/data3/test/results/test_result/quality_control/P25-D16-1/P25-D16-1.clean.1.fastq.gz"
clean_R2="/datapool/home/fengbx/data3/test/results/test_result/quality_control/P25-D16-1/P25-D16-1.clean.2.fastq.gz"
rm_rRNA_R="/datapool/home/fengbx/data3/test/results/test_result/remove_rRNA/fastq/P25-D16-1/P25-D16-1.rm_rRNA"
log="/datapool/home/fengbx/data3/test/scripts/logs/debug.log"


cd /datapool/home/fengbx/data3/test/results/test_result/remove_rRNA/fastq/
bowtie2 -x ${rRNA_index} -1 ${clean_R1} -2 ${clean_R2} --un-conc-gz ${rm_rRNA_R}.fq.gz  -p 64 -S ${rm_rRNA_R}.sam &>${log}
find /datapool/home/fengbx/data3/test/results/test_result/remove_rRNA/fastq/P25-D16-1/ -name "*.rm_rRNA.fq.1.gz" -exec mv {} .rm_rRNA.1.fq.gz \;