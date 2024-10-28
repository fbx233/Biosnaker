
library(ggheatmap)
library(tidyr)
set.seed(123)
df <- matrix(runif(225,0,10),ncol = 15)
colnames(df) <- paste("sample",1:15,sep = "")
rownames(df) <- sapply(1:15, function(x)paste(sample(LETTERS,3,replace = F),collapse = ""))
df[1:4,1:4]


row_metaData <- data.frame(exprtype=sample(c("Up","Down"),15,replace = T),
                           genetype=sample(c("Metabolism","Immune","None"),15,replace = T))
rownames(row_metaData) <- rownames(df)
col_metaData <- data.frame(tissue=sample(c("Normal","Tumor"),15,replace = T),
                           risklevel=sample(c("High","Low"),15,replace = T))
rownames(col_metaData) <- colnames(df)
exprcol <- c("#EE0000FF","#008B45FF" )
names(exprcol) <- c("Up","Down")
genecol <- c("#EE7E30","#5D9AD3","#D0DFE6FF")
names(genecol) <- c("Metabolism","Immune","None")
tissuecol <- c("#98D352","#FF7F0E")
names(tissuecol) <- c("Normal","Tumor")
riskcol <- c("#EEA236FF","#46B8DAFF")
names(riskcol) <- c("High","Low")
col <- list(exprtype=exprcol,genetype=genecol,tissue=tissuecol,risklevel=riskcol)
text_rows <- sample(rownames(df),3)

p<- ggheatmap(df,cluster_rows = T,cluster_cols = T,scale = "row",
              text_show_rows = text_rows,
              cluster_num = c(3,3),
              tree_color_rows = c("#008B45FF","#631879FF","#008280FF"),
              tree_color_cols = c("#1F77B4FF","#FF7F0EFF","#2CA02CFF"),
              annotation_rows = row_metaData,
              annotation_cols = col_metaData,
              annotation_color = col
)
ggsave(p, file ="~/test.png")
