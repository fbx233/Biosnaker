## About Biosnaker
Biosnaker simplifies bioinformatic process by building, versioning, and hosting of your docs, automatically.

## Why use this tool？


```txt
├── data
├── results
└── scripts
    ├── bowtie2.sh
    ├── cluster.json
    ├── cluster.yaml
    ├── config
    │   ├── bioinfo.yaml
    │   ├── cluster.json
    │   ├── cluster.yaml
    │   ├── config.yaml
    │   └── jobscript.sh
    ├── config.yaml
    ├── Count_matrix.py
    ├── Count_matrix.sh
    ├── debug1.sh
    ├── deseq2.R
    ├── envs
    │   ├── bioinfo.yaml
    │   ├── py.yaml
    │   └── sc.yaml
    ├── generate_sample_info.py
    ├── jobscript.sh
    ├── logs
    │   └── debug.log
    ├── mamba.sh
    ├── merge_counts.yaml
    ├── modify_files.sh
    ├── py.yaml
    ├── README.md
    ├── refer
    │   ├── genetic_screen_human_metabolism_genes.txt
    │   ├── human_metabolic_genes.csv
    │   └── pathway_class_information.csv
    ├── remove_rRNA.sh
    ├── Rscripts
    │   ├── heatmap.r
    │   ├── multi_group_DEA.R
    │   ├── PLSDA.R
    │   ├── README.md
    │   └── RNA-seq-DEA-function-repos.R
    ├── run_snakemake.pbs
    └── snakefile
```

## First time here?
**Biosnaker** is a [Snakemake workflow](https://snakemake.readthedocs.io/en/stable/index.html), aimed at performing a typical RNA-seq workflow in a reproducible, automated, and partially contained manner. It is implemented such that alternative or similar analysis can be added or removed.

Biosnaker consists of a `Snakefile`, a [`conda`](https://conda.io/docs/) environment file (`envs/environment.yaml`) a configuration file (`config.yaml`) and a set of `R` scripts, to perform quality control, preprocessing and differential expression analysis of RNA-seq data. The output can be integrated with Python/R scripts for further analysis, allowing for more in-depth exploration and sharing of the results.

The config.yaml file contains several important parameter configurations that facilitate the automation of downstream data analysis workflows. The key parameters are as follows:

- **Sample List**: The samples section includes identifiers for a series of samples that are the focus of the data analysis. These identifiers are categorized into different classes, such as P7-DTP, P7-Par, and P7-Reg, with each category further divided into sub-samples like D12, D4, D9, etc. Each sample sequence has multiple replicates, such as "1", "2", "3", "4", to ensure the reliability and reproducibility of the experiments.

- **Data Path**: The path specifies where the raw omics data is stored, set to /public5/home/t6s001510/data1/omics_data/Rawdata/.

- **Output Path**: The output_path is the destination for the analysis results, configured as /public5/home/t6s001510/data1/omics_data/outputs/.

- **Genome Index Path**: The genome_index_path is used to store the genome indexes, supporting alignment for transcriptome data, specified as /public5/home/t6s001510/data1/omics_data/ref/index/STAR/genome.

- **rRNA Index Path**: The rRNA_index provides the reference index for rRNA to filter out rRNA fragments from the data, located at /public5/home/t6s001510/data1/omics_data/ref/index/rRNA/Homo_sapiens.rRNA.

- **GTF File Path**: The GTF parameter specifies the annotation file for precise gene localization and functional annotation, at /public5/home/t6s001510/data1/omics_data/ref/index/GTF/Homo_sapiens.GRCh38.110.gtf.

- **Count Matrix Script**: The count_matrix indicates the script path for generating the gene expression count matrix, located at /public5/home/t6s001510/data1/omics_data/scripts/Count_matrix.py.

- **Sample Information Script**: The sample_info provides the script path for generating sample information, helping organize metadata about the samples, at /public5/home/t6s001510/data1/omics_data/scripts/generate_sample_info.py.

- **Pathway Information Reference File**: pathway_info contains the file path with information about biological pathway classifications, available for reference in the analysis, located at /public5/home/t6s001510/data1/omics_data/scripts/refer/pathway_class_information.csv.

- **Human Metabolic Genes Reference File**: The hm_gene lists human genes involved in metabolic processes, specified at /public5/home/t6s001510/data1/omics_data/scripts/refer/human_metabolic_genes.csv.


```sh
samples:
- Sample1
- Sample2
  
path: "Your_input_data"
sc_fastq_path: "Your_input_scRNA_data"
output_path: "output_path"
genome_index_path: "Your_STAR_genome_path"
rRNA_index: "rRNA_path"
GTF: "gtf_doc_path"
count_matrix: "scripts/Count_matrix.py"
sample_info: "scripts/generate_sample_info.py"
pathway_info: "refer/pathway_class_information.csv"
hm_gene: "refer/human_metabolic_genes.csv"
cellranger_transcriptome: "path_to_cellranger_transcriptome_GRCh38-2020-A"
```

## How-to use it？
By default, the pipeline performs all the steps shown in the below. However, you can turn off any combination of the light-colored steps (e.g `STAR` alignment or `DRIMSeq` analysis) in the `config.yaml` file.

![alt text](SOP.png)

Advanced use: If you prefer other software to run one of the outlined steps (e.g. `limma-voom` over `DESeq2`, or `kallisto` over `Salmon`), you can use the software of your preference provided you have your own script(s), and change some lines within the `Snakefile`. If you think your "custom rule" might be of use to a broader audience, let us know by opening an issue.
