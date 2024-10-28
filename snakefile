import yaml
import os
import glob

# Load configuration file
configfile: "envs/config.yaml"

# Get the list of samples
samples = config["samples"]

# Rule to generate the final output file
rule all:
    input:
        config["output_path"] + "read_counts/result/sample_info.csv",
        expand(config["output_path"] + "cellranger/{sample}/outs/filtered_feature_bc_matrix", sample=sc_samples)

rule fastqc:
    input:
        R1=config["path"] + "{sample}/{sample}_R1.fq.gz",
        R2=config["path"] + "{sample}/{sample}_R2.fq.gz"
    output:
        config["output_path"] + "fastqc/{sample}/{sample}_R1_fastqc.html",
        config["output_path"] + "fastqc/{sample}/{sample}_R2_fastqc.html"
    log:
        config["output_path"] + "logs/fastqc/{sample}.log"
    conda:
        "envs/bioinfo.yaml"
    shell:
        """
        mkdir -p {config[output_path]}fastqc/{wildcards.sample}
        fastqc {input.R1} {input.R2} --outdir {config[output_path]}fastqc/{wildcards.sample} --noextract &> {log}
        """

rule fastp:
    input:
        R1=config["path"] + "{sample}/{sample}_R1.fq.gz",
        R2=config["path"] + "{sample}/{sample}_R2.fq.gz"
    output:
        clean_R1=config["output_path"] + "quality_control/{sample}/{sample}.clean.1.fastq.gz",
        clean_R2=config["output_path"] + "quality_control/{sample}/{sample}.clean.2.fastq.gz",
        json=config["output_path"] + "quality_control/{sample}/{sample}.json",
        html=config["output_path"] + "quality_control/{sample}/{sample}.html"
    log:
        config["output_path"] + "logs/fastp/{sample}.log"
    conda:
        "envs/bioinfo.yaml"
    shell:
        """
        mkdir -p {config[output_path]}quality_control/{wildcards.sample}
        fastp -i {input.R1} -I {input.R2} -o {output.clean_R1} -O {output.clean_R2} --thread=16 -l 15 -j {output.json} -h {output.html} &> {log}
        """

rule remove_rRNA:
    input:
        clean_R1=config["output_path"] + "quality_control/{sample}/{sample}.clean.1.fastq.gz",
        clean_R2=config["output_path"] + "quality_control/{sample}/{sample}.clean.2.fastq.gz"
    output:
        rm_rRNA_R1=config["output_path"] + "remove_rRNA/fastq/{sample}/{sample}.rm_rRNA.fq.1.gz",
        rm_rRNA_R2=config["output_path"] + "remove_rRNA/fastq/{sample}/{sample}.rm_rRNA.fq.2.gz"
    params:
        rRNA_index=config["rRNA_index"]
    log:
        config["output_path"] + "logs/remove_rRNA/{sample}.log"
    conda:
        "envs/bioinfo.yaml"
    shell:
        """
        mkdir -p {config[output_path]}remove_rRNA/fastq/{wildcards.sample}
        bowtie2 -x {params.rRNA_index} -1 {input.clean_R1} -2 {input.clean_R2} --un-conc-gz {config[output_path]}remove_rRNA/fastq/{wildcards.sample}/{wildcards.sample}.rm_rRNA.fq.gz -p 32 -S {config[output_path]}remove_rRNA/fastq/{wildcards.sample}/{wildcards.sample}.rm_rRNA.sam &> {log}
        """

rule mapping:
    input:
        rm_rRNA_R1=config["output_path"] + "remove_rRNA/fastq/{sample}/{sample}.rm_rRNA.fq.1.gz",
        rm_rRNA_R2=config["output_path"] + "remove_rRNA/fastq/{sample}/{sample}.rm_rRNA.fq.2.gz"
    output:
        sorted_bam=config["output_path"] + "mapping_expression/{sample}/{sample}Aligned.sortedByCoord.out.bam"
    params:
        genome_index_path=config["genome_index_path"]
    log:
        config["output_path"] + "logs/mapping/{sample}.log"
    conda:
        "envs/bioinfo.yaml"
    shell:
        """
        mkdir -p {config[output_path]}mapping_expression/{wildcards.sample}
        STAR --runThreadN 32 --limitBAMsortRAM 20000000000 --outFilterType BySJout --outFilterMismatchNmax 10 --genomeDir {params.genome_index_path} --readFilesIn {input.rm_rRNA_R1} {input.rm_rRNA_R2} --readFilesCommand 'zcat' --outFileNamePrefix {config[output_path]}mapping_expression/{wildcards.sample}/{wildcards.sample} --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM GeneCounts --outSAMattributes All --outSAMstrandField intronMotif --outBAMcompression 6 --outReadsUnmapped Fastx &> {log}
        samtools sort -T {config[output_path]}mapping_expression/{wildcards.sample}/{wildcards.sample}Aligned.out.sorted -o {output.sorted_bam} {config[output_path]}mapping_expression/{wildcards.sample}/{wildcards.sample}Aligned.out.bam
        samtools index {output.sorted_bam}
        """

rule featureCounts:
    input:
        sorted_bam=config["output_path"] + "mapping_expression/{sample}/{sample}Aligned.sortedByCoord.out.bam"
    output:
        featurecounts_txt=config["output_path"] + "read_counts/{sample}/{sample}.featurecounts.txt",
        featurecounts_all_txt=config["output_path"] + "read_counts/{sample}/{sample}.featurecounts.all.txt"
    params:
        GTF=config["GTF"]
    log:
        config["output_path"] + "logs/featureCounts/{sample}.log"
    conda:
        "envs/bioinfo.yaml"
    shell:
        """
        mkdir -p {config[output_path]}read_counts/{wildcards.sample}
        featureCounts -T 32 -s 0 -p -t CDS -g gene_id -a {params.GTF} -o {output.featurecounts_txt} {input.sorted_bam} &> {log}
        featureCounts -T 32 -s 0 -p -t exon -g gene_id -a {params.GTF} -o {output.featurecounts_all_txt} {input.sorted_bam} &> {log}
        """

rule merge_counts:
    input:
        expand(config["output_path"] + "read_counts/{sample}/{sample}.featurecounts.all.txt", sample=samples)
    output:
        config["output_path"] + "read_counts/result/count.all.txt"
    log:
        config["output_path"] + "logs/merge_counts.log"
    conda:
        "envs/bioinfo.yaml"
    shell:
        """
        mkdir -p {config[output_path]}read_counts/result
        mkdir -p {config[output_path]}deseq2/data
        python {config[count_matrix]} --input_files {input} --output_file {output}
        cp {output} {config[output_path]}deseq2/data
        """
        

rule generate_sample_info:
    input:
        config["output_path"] + "read_counts/result/count.all.txt"
    output:
        config["output_path"] + "read_counts/result/sample_info.csv"
    log:
        config["output_path"] + "logs/sample_info.log"
    conda:
        "envs/bioinfo.yaml"
    shell:
        """
        python3 {config[sample_info]} --input_files {input} --output_file {output} &> {log}
        cp {output} {config[output_path]}deseq2/data
        """
        
rule RNA_seq_DEA:
    input:
        matrix=config["output_path"] + "read_counts/result/count.all.txt"
        group=config["output_path"] + "read_counts/result/sample_info.csv"
    output:
        config["output_path"] + "deseq2/results/"
    log:
        config["output_path"] + "logs/deseq2.log"
    conda:
        "envs/sc.yaml"
    shell:
        """
        mkdir -p {config[output_path]}deseq2/results
        Rscript ./Rscripts/multi_group_DEA.R -i {input.matrix} -s {input.group} -o {output} -p {config[pathway_info]}
        """
Rscript ./Rscripts/multi_group_DEA.R -i ${matrix} -s ${group} -o ${output} -p ${pathway_info} -m ${hm_gene}
Rscript ./Rscripts/PLSDA.R -i ${matrix} -s ${group} -o ${output} -n PD25




rule cellranger_count:
    input:
        fastq_dir = config["sc_fastq_path"] + "{sample}",
    output:
        matrix_dir = config["output_path"] + "cellranger/{sample}/outs/filtered_feature_bc_matrix"
    params:
        transcriptome = config["cellranger_transcriptome"],
        sample = "{sample}",
        outdir = config["output_path"] + "cellranger/{sample}"
    log:
        config["output_path"] + "logs/cellranger/{sample}.log"
    shell:
        """
        mkdir -p {params.outdir}
        mkdir -p {config[output_path]}logs/cellranger/
        /datapool/home/fengbx/cellranger-8.0.1/cellranger count \
            --id={wildcards.sample} \
            --transcriptome={params.transcriptome} \
            --fastqs={input.fastq_dir} \
            --sample={params.sample} \
            --localcores=16 \
            --localmem=64 \
            --nosecondary \
            --output-dir={params.outdir} \
            &> {log}
        """
