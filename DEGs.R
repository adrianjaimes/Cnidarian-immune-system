################ Differential expression analysis ################


##__________________________________DESeq2________________________________________##


library(DESeq2)
library(ggplot2)
library(reshape2)
library(dplyr)

counttable <- read.csv("count_matrix.csv",header = TRUE, row.names = 1)
counttable
library(forestmangr)
counttable <- round_df(counttable, rf= "trunc")
counttable

sample_information <- read.table('Sample_Information.csv',header=T)
sample_information

all(rownames(sample_information) %in% colnames(counttable))
all(rownames(sample_information) == colnames(counttable))


filtered_samples <- sample_information[rownames(sample_information) %in% c("a5","a6","a7","a8","a9","a10","a27","a28","a29","a30","a31","a32"), ]
filtered_counts <- counttable[, row.names(filtered_samples)]

filtered_samples$Treatment <- as.factor(filtered_samples$Treatment)
filtered_samples$Genotype <- as.factor(filtered_samples$Genotype)

filtered_samples

#comparison
dds <- DESeqDataSetFromMatrix(countData = filtered_counts,        
                              colData = filtered_samples,
                              design = ~ Treatment)


keep <- rowSums(counts(dds)) >= 10 #Keep rows with at least ten reads
dds <- dds[keep,]
dds
nrow(dds)


################ Transformation-exploration and visualization (PCA, Heatmaps)################

#vst Transformation

vsd <- vst(dds, blind = FALSE)
vsd_1 <- vst(dds, blind = TRUE)
head (assay(vsd), 3)
vsd
colData(vsd)

#rlog Transformation

rld <- rlog(dds, blind = FALSE)
rld
rld_1 <- rlog(dds, blind = TRUE)
head(assay(rld_1), 3)
rld_1
colData(rld_1)

##Extracting the matrix of transformed data:

vsd_data_matrix <- assay(vsd)
vsd_data_matrix_1 <- assay(vsd_1)

rld_data_matrix <- assay(rld)
rld_data_matrix_1 <- assay (rld_1)
rld_data_matrix_1


##Computing the pairwise correlation values between the samples:

vsd_data_correlations <- cor(vsd_data_matrix)
vsd_data_correlations_1 <- cor(vsd_data_matrix_1)

rld_data_correlations <- cor(rld_data_matrix)
rld_data_correlations_1 <- cor(rld_data_matrix_1)


## Drawing a heatmap of sample to sample distances:

library(pheatmap)
library(dplyr)


pheatmap(vsd_data_correlations, annotation = select (filtered_samples, Genotype))
pheatmap(vsd_data_correlations_1, annotation = select (filtered_samples, Genotype))

pheatmap(rld_data_correlations, annotation = select (filtered_samples, Genotype))
pheatmap(rld_data_correlations_1, annotation = select (filtered_samples, Genotype))

library(ggrepel)

pcaData <- plotPCA(rld, intgroup=c("Genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, shape=dds$Treatment, fill=Genotype, color=Genotype)) +
  geom_point(size=4, na.rm = TRUE) +
  scale_shape_manual(values=c(21,24)) +
  scale_fill_manual(values = c("pink2", "yellow")) +  # Fills color inside the shape
  scale_color_manual(values = c("pink2", "yellow")) + # Changes the border color
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()


################Differential expression analysis (Standard Analysis)####################

##__________________________________DESeq2________________________________________##

dds
dds_1 <- DESeq(dds)
dds_1
resultsNames(dds_1)

res <- results(dds_1, contrast = c('Treatment', 'PIC','NACL'))
is.data.frame(res)
res_df <- as.data.frame(res)
is.data.frame(res_df)
res_df

summary(res)
filtered <- res_df %>% filter(res_df$padj < 0.05)
filtered_1 <- filtered %>% filter(abs(filtered$log2FoldChange) >= 1)

dim(res_df)
dim(filtered)
dim(filtered_1)

write.csv(res_df, 'InterestGene_Deseq2_all.csv')
write.csv(filtered_1, 'Deseq2_InterestGene_DE.csv')

upregulated_WT <- filtered_1 %>% filter(log2FoldChange >= 1)
dim(upregulated_WT)
downregulated_WT <- filtered_1 %>% filter(log2FoldChange <= -1)
dim(downregulated_WT)

write.csv(upregulated_WT, "upregulated_genes_deseq2.csv", row.names = TRUE)
write.csv(downregulated_WT, "downregulated_genes_deseq2.csv", row.names = TRUE)

upregulated_genes <- data.frame(gene = rownames(upregulated_WT))
write.csv(upregulated_genes, "upregulated_genes_for_clusterProfiler_deseq2.csv", row.names = FALSE, quote = FALSE)

downregulated_genes <- data.frame(gene = rownames(downregulated_WT))
write.csv(downregulated_genes, "downregulated_genes_for_clusterProfiler_deseq2.csv", row.names = FALSE, quote = FALSE)

##__________________________________EdgeR________________________________________##

library('edgeR')
library('statmod')
library('limma')


dgeFull <- DGEList(filtered_counts, group = filtered_samples$Treatment)
dim(dgeFull)

keep_1 <- rowSums(dgeFull$counts) >= 10
dgeFull <- dgeFull[keep_1, ]
head(dgeFull$counts)

dgeFull <- calcNormFactors(dgeFull, method="TMM")
dgeFull$samples
design <- model.matrix(~filtered_samples$Treatment)
dgeFull <- estimateCommonDisp(dgeFull, design)
dgeFull <- estimateTagwiseDisp(dgeFull)
dgeFull

dgeExactTest  <- exactTest(dgeFull, pair=c("NACL","PIC"))
dgeExactTest

resExactTest <- topTags(dgeExactTest, n = nrow(dgeExactTest$table)) # toptags generates FDR
head(resExactTest$table)
write.table(resExactTest, file = "InterestGene_EdgeR_all.csv", sep = ",")

selectedET <- resExactTest$table$FDR < 0.05 & abs(resExactTest$table$logFC) >= 1
selectedET <- resExactTest$table[selectedET, ]
nrow(selectedET)
head(selectedET)

upregulated_genes_EdgeR <- selectedET[selectedET$logFC >= 1, ]
dim(upregulated_genes_EdgeR)
write.table(upregulated_genes_EdgeR, file = "Upgenes_Edger.csv", sep = ",")

downregulated_genes_EdgeR <- selectedET[selectedET$logFC <= -1, ]
dim(downregulated_genes_EdgeR)
write.table(downregulated_genes_EdgeR, file = "Downgenes_Edger.csv", sep = ",")


##################Intersection EdgeR-Deseq2 or any two o more lists of genes###########################

upregulated_list1 <- rownames(upregulated_WT)  # Upregulated genes from dataset 1 (deseq2)
upregulated_list2 <- rownames(upregulated_genes_EdgeR)  # Upregulated genes from dataset 2 (EdgeR)

upregulated_intersection_WT <- intersect(upregulated_list1, upregulated_list2)
length(upregulated_intersection_WT)

write.csv(upregulated_intersection_WT, "upregulated_genes_Intersection.csv", row.names = FALSE, quote = FALSE)

#-------------------------------------------------------------------------------------------------------

downregulated_list1 <- rownames(downregulated_WT)  # Upregulated genes from dataset 1 (Deseq2)
downregulated_list2 <- rownames(downregulated_genes_EdgeR)  # Upregulated genes from dataset 2 (EdgeR)

downregulated_intersection_WT <- intersect(downregulated_list1 , downregulated_list2)
length(downregulated_intersection_WT)

write.csv(downregulated_intersection_WT, "downregulated_genes_Intersection.csv", row.names = FALSE, quote = FALSE)


#####################################################################################################################


##__________________________________Goterms Analysis________________________________________##

#Clusterprofile

###Overrepresentation analysis

library (clusterProfiler)

gene.vector <- read.csv(file = 'genes.csv',header=TRUE)
head(gene.vector)
genes <- gene.vector$Genes #Converter one column to vector
is.vector(genes)

TermGene  <- read.csv(file = 'GoName_PM.csv',header=TRUE, check.names=FALSE)
TermGene
is.data.frame(TermGene)

TermName  <- read.csv(file = 'Goterms_PM.csv',header=TRUE, check.names=FALSE)
TermName
is.data.frame(TermName)

Results <- enricher(genes, pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, TERM2GENE = TermGene, TERM2NAME = TermName)
write.csv(Results, file = "genes_result.csv", quote = FALSE, row.names = FALSE)

barplot(Results)
barplot(Results, showCategory=15)


#####################################################################################################################