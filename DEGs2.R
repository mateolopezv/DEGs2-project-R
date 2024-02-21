#DEGS2
# Load the count data
count_table <- read.csv("C:/Users/ignac/Desktop/R/Curso/Basics/raw_counts.tsv", sep = "\t", header = TRUE, row.names = 1)
row.names(count_table)
head(count_table)
#Load the sample data
sample_info <- read.csv("C:/Users/ignac/Desktop/R/Curso/Basics/design.tsv", sep = "\t", header = TRUE, row.names = 1)
row.names(sample_info)
head(sample_info)
#Set factors
factors <- factor(sample_info$Group)
groups <- unique(sample_info$Group)
groups
factors
#Create DESeq object
dds <- DESeqDataSetFromMatrix(countData = count_table, colData = sample_info, design = ~ Group)
#Set reference factor
dds$Group <- relevel(dds$Group, ref = "control")
#Filter genes with low counts (n >= 10)
keep <- rowSums(counts(dds) >= 10) >= min(table(sample_info$Group))
keep2 <- as.matrix(keep)
dds <- dds[keep,]
#Perform statistical test
dds <- DESeq(dds, test = "Wald", sfType = "poscount")
#Get result
deseq_result <- results(dds)
deseq_result
deseq_result <- as.data.frame(deseq_result)

deseq_result$GeneName <- row.names(deseq_result)
names(deseq_result)
deseq_result <- subset(deseq_result,
                       select = c("GeneName", "padj", "pvalue", "stat", "lfcSE", "log2FoldChange", "baseMean"))
names(deseq_result)

#Extract genes with padj < 0.05 and logFoldChange >1 and <-1
deg <- subset(deseq_result, padj<0.05 & abs(log2FoldChange>1))
deg <- deg[order(deg$padj),]
View(deg)

#Gene expression data viz
plotDispEsts(dds, main = "Dispersion Estimates")
#Histogram of p-values
hist(deseq_result$padj)
#Volcano Plot
#Set color
old.pal <- palette(c("#00BFFF", "#FF3030"))
#Set margin size
par(mar=c(4,4,2,1), cex.main=1.5)
#Set title
title = paste(groups[1], "vs", groups[2])
#Plot values
plot(deseq_result$log2FoldChange, -log10(deseq_result$padj), main=title, xlab="log2FC", ylab="-log10(padj)", pch=20, cex=0.7)
with(subset(deseq_result, padj < 0.05 & abs(log2FoldChange) >= 1),
     points(log2FoldChange, -log10(padj), pch=20, cex=0.7, col=(sign(log2FoldChange)+3)/2))
legend("bottomleft", title = paste("Padj<", 0.05, sep=""),
       legend = c("Down","Up"), pch = 20, col = 1:2)

#Heatmap
normalized_counts <- counts(dds, normalized = TRUE)
top_hits <- deseq_result[order(deseq_result$padj),][1:10,]
top_hits <- row.names(top_hits)
cal_z_score <- function(x) {(x - mean(x))/ sd(x)}
zscore_all <- t(apply(normalized_counts,1,cal_z_score))
zscore_subset <- zscore_all[top_hits,]
pheatmap(zscore_subset)






