#!/usr/bin/env Rscript

library("DESeq2")

# Read in the count data and col data
countData <- read.table("count_data.txt", header=TRUE, row.names=1)
colData <- read.table("col_data.txt", header=TRUE, row.names=1)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~condition)

# Run the DESeq2 pipeline
dds <- DESeq(dds)

# Get the results
res <- results(dds)

# Save the results to a CSV file
write.csv(as.data.frame(res), file="differential_expression_results.csv")
