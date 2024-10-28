#### 转录组降维分析 #### 
.libPaths("/datapool/home/fengbx/mambaforge/envs/fbx1/lib/R/library")
library(optparse)
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(mixOmics)

# 设置命令行参数
option_list <- list(
  make_option(c("-i", "--input_files"), type="character", default="Data/count.all.txt", help="Path to the input gene expression matrix file", metavar="character"),
  make_option(c("-s", "--sample_info_path"), type="character", default="Data/Sample_information_without_Regrowth-1.csv", help="Path to the sample information file", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL, help="Directory to save the output files", metavar="character"),
  make_option(c("-n", "--name"), type="character", default="Sample", help="All samples name", metavar="character")
  make_option(c("-l", "--log2fc"), type="number", default=1.5, help="log2 fold change", metavar="number")
  make_option(c("-p", "--padj"), type="number", default=0.05, help="adjust P value", metavar="number")
)


opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)


# 转录组的降维分析利用DESeq2分析中的Count转换后的数据
# 基因表达水平在不同样本之间存在一定差异，DESeq2提出两种Count值转换算法, rlog和VST转换。
# rlog适用于样本量<30的数据。
#setwd("~/Desktop/P25/P25-PLSDA/")
rm(list=ls())
# 设置FDR的阈值
padjThread=opt$padj
# 设置读段总数阈值
readssumThread=10
# 设置Foldchange变化阈值
L2FCThread=opt$log2fc
# 设置文件读取路径
input_dir = opt$input_files
sample_info_path = opt$sample_info_path

sample_info <- read.csv(sample_info_path) 
raw_count <- read.table(input_dir,sep='\t',header=T,check.names = FALSE)
row.names(raw_count)<-raw_count[,1]
raw_count <-raw_count[,-1]
countdata <-raw_count
countdata <- countdata[rowSums(countdata)>readssumThread,]
countdata<-countdata[,sample_info$sample]

### 2.1) provide sample information
condition_merge<-factor(sample_info$group)
### 2.2) Associate grouping information with sample names.
colData <- data.frame(row.names=sample_info$sample,condition_merge)

### 3) DESeq analysis
### 3.1) Create a DESeq2Dataset object.
dds <- DESeqDataSetFromMatrix(countdata,colData,design = ~condition_merge)
### 3.2) 基因表达水平在不同样本之间存在一定差异，DESeq2提出两种Count值转换算法，rlog和VST转换。
#### rlog适用于样本量<30的数据。
rld <- rlog(dds)
rlog_data <- assay(rld)
rlog_data <- as.data.frame(rlog_data)

pca_res <- prcomp(t(rlog_data))
pca_data <- data.frame(Sample = rownames(pca_res$x),
                       PC1 = pca_res$x[, 1],
                       PC2 = pca_res$x[, 2],
                       Group = colData)  # 根据你的数据调整

library(ggsci)
p1 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition_merge,fill = condition_merge)) +
  geom_point(alpha = 0.5,size=4) +  # 增加点的透明度
  stat_ellipse(aes(group = condition_merge), 
               type = "t", 
               level = 0.8, 
               size = 1, 
               linetype = "dashed",  # 使用虚线
               position = position_nudge(x = 0.01, y = 0.01)) +  # 轻微移动椭圆位置，以防被点覆盖
  ggtitle("PCA Plot") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.line = element_line(color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )+scale_color_npg()+scale_fill_npg()
p1
# 保存PCA plot
ggsave(paste0(opt$output_dir, opt$name, "PCA.pdf"),plot = p1,width = 8, height = 6, dpi = 600)


# 执行PLS-DA分析
plsda_res <- plsda(t(rlog_data), colData$condition)  # 使用你的分组信息
pdf(paste0(opt$output_dir, opt$name, "PLSDA.pdf"),width=11.69,height=8.27)
# 绘制PLS-DA分数图
plotIndiv(plsda_res,
          ind.names = FALSE,
          legend = TRUE,
          title = "PLS-DA Plot",
          ellipse = TRUE,
          legend.title = "Group")
dev.off()