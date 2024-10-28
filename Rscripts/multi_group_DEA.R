.libPaths("/datapool/home/fengbx/mambaforge/envs/fbx1/lib/R/library")
library(optparse)
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(openxlsx)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggrepel)
library(tools)

# 函数定义 (假设这些函数在一个外部文件中)
source("Rscripts/RNA-seq-DEA-function-repos.R")

# 设置命令行参数
option_list <- list(
  make_option(c("-i", "--input_files"), type="character", default=NULL, help="Path to the input gene expression matrix file", metavar="character"),
  make_option(c("-s", "--sample_info_path"), type="character", default=NULL, help="Path to the sample information file", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL, help="Directory to save the output files", metavar="character"),
  make_option(c("-p", "--pathway_list_path"), type="character", default="Data/pathway_class_information.csv", help="Path to the pathway class information file", metavar="character"),
  make_option(c("-m", "--human_metabolic_gene"), type="character", default="Data/human_metabolic_genes.csv", help="Path to the human metabolic gene file", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# 检查必需的输入参数是否存在
if (is.null(opt$input_files) | is.null(opt$sample_info_path) | is.null(opt$output_dir)) {
  print_help(opt_parser)
  stop("Input files, sample information path, and output directory must be provided", call.=FALSE)
}

# 参数设置
input_files <- opt$input_files
sample_info_path <- opt$sample_info_path
pathway_list_path <- opt$pathway_list_path
human_metabolic_gene <- opt$human_metabolic_gene
plot_name <- "P25 RNA-seq"
plot_title <- ""
xlab <- "log2(FoldChange)"
ylab <- "-log10(FDR)"
padjThread <- 0.05
readssumThread <- 10
L2FCThread <- 1.5

# 检查文件扩展名并读取样本信息文件
file_extension <- tools::file_ext(sample_info_path)
if (file_extension == "csv") {
  sample_info <- read.csv(sample_info_path)
} else if (file_extension == "txt") {
  sample_info <- read.delim(sample_info_path)
} else {
  stop("Unsupported file format. Please provide a csv or txt file for sample information.")
}

# 提取组别
groups <- unique(sample_info$group)

# 生成两两比较的组合
comparisons <- combn(groups, 2, simplify = FALSE)

# 创建输出目录
if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = TRUE)
}

# 设置工作目录
setwd(opt$output_dir)

# 确保结果目录存在
result_dir <- file.path(opt$output_dir, "Result")
if (!dir.exists(result_dir)) {
  dir.create(result_dir)
}

# 确保所有必要的文件存在
if (!file.exists(input_files)) {
  stop(paste("Input file does not exist:", input_files))
}
if (!file.exists(sample_info_path)) {
  stop(paste("Sample information file does not exist:", sample_info_path))
}
if (!file.exists(pathway_list_path)) {
  stop(paste("Pathway list file does not exist:", pathway_list_path))
}
if (!file.exists(human_metabolic_gene)) {
  stop(paste("Human metabolic gene file does not exist:", human_metabolic_gene))
}

# 运行多组差异表达分析和KEGG通路富集分析
for (i in seq_along(comparisons)) {
  comparison <- comparisons[[i]]
  output_name <- paste0("T", sprintf("%02d", i))
  
  cat("Running comparison:", comparison, "\n")
  
  # 差异表达分析
  allNmr <- Differential_expression_analysis(
    input_files,
    sample_info_path,
    comparison,
    comparison[1], # 选择第一个组作为对照组
    output_name,
    human_metabolic_gene,
    plot_title,
    xlab,
    ylab,
    padjThread,
    readssumThread,
    L2FCThread
  )
  
  all_gene <- allNmr$all
  mr_gene <- allNmr$mr
  all_DE_gene <- all_gene %>% filter(sig != "none")
  mr_DE_gene <- mr_gene %>% filter(sig != "none")
  
  # KEGG富集分析
  all_kegg_enrich <- KEGGpathwayEnirch(all_DE_gene$SYMBOL, "SYMBOL", "ENTREZID", "org.Hs.eg.db", "hsa", 0.05, 0.05)
  mr_kegg_enrich <- KEGGpathwayEnirch(mr_DE_gene$SYMBOL, "SYMBOL", "ENTREZID", "org.Hs.eg.db", "hsa", 0.05, 0.05)
  
  # 确保结果文件路径正确
  output_all_path <- file.path(result_dir, paste0(output_name, "-all-KEGG-Enrichment-plot.pdf"))
  output_mr_path <- file.path(result_dir, paste0(output_name, "-mr-KEGG-Enrichment-plot.pdf"))
  
  # 绘制KEGG富集图
  KEGG_enrich_plot(all_kegg_enrich, output_all_path)
  KEGG_enrich_metabolism_related_plot(all_kegg_enrich, output_all_path, pathway_list_path)
  KEGG_enrich_plot(mr_kegg_enrich, output_mr_path)
  KEGG_enrich_metabolism_related_plot(mr_kegg_enrich, output_mr_path, pathway_list_path)
}


