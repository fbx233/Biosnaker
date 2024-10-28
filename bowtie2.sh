#!/bin/bash

# 定义参数
path="/datapool/home/fengbx/data3/test/data/"
output_path="/datapool/home/fengbx/data3/test/results/test_result"
genome_index_path="/datapool/home/fengbx/data1/ref_genome/index/STAR/genome"
rRNA_index="/datapool/home/fengbx/data1/ref_genome/index/rRNA/Homo_sapiens.rRNA"
GTF="/datapool/home/fengbx/data1/ref_genome/index/GTF/Homo_sapiens.GRCh38.110.gtf"

# 确保输出目录存在
mkdir -p $output_path/remove_rRNA/fastq/
mkdir -p $output_path/logs/remove_rRNA/

# 执行bowtie2命令
bowtie2 -x $rRNA_index -1$pathclean_R1 -2 $pathclean_R2 --un-conc-gz$output_path/remove_rRNA/fastq/{sample}/{sample}.rm_rRNA.fq -p 16 -S $output_path/remove_rRNA/fastq/{sample}/{sample}.aligned_rRNA.sam &>$output_path/logs/remove_rRNA/{sample}.log

