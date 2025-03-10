
#Adipose v Myocardium DEG=======================================================
#===============================================================================
library(dplyr)
setwd("C:\\Users\\User\\Desktop")
counts = read.csv(".\\omics\\datasets\\Full_Counts.csv")
rownames(counts) = counts[,1]
colnames(counts) = substring(colnames(counts), 7, 10)
counts = counts[,c(2:81)]

A <- counts[, grepl("A", names(counts))] #subset counts DF to columns matching "A" or "M"
M <- counts[, grepl("M", names(counts))]
A_M_metadata <- data.frame(DONOR = c(colnames(A),colnames(M)), TISSUE = c(rep("A", 40), rep("M", 40)))

counts = select(counts, contains(unlist(A_M_metadata[,1]))) 
all(colnames(counts) %in% A_M_metadata$DONOR) #identical
all(colnames(counts) == A_M_metadata$DONOR)   #same order
dim(A_M_metadata)
dim(counts)

library(DESeq2)
alpha_p = .05
dds = DESeqDataSetFromMatrix(countData = counts, 
                             colData = A_M_metadata, 
                             design = ~TISSUE)
dds = dds[rowSums(counts(dds))>1]
dds$TISSUE = relevel(dds$TISSUE, ref="A") #Set reference level manually
dds = DESeq(dds)

resultsNames(dds)
result1 = results(dds, alpha = alpha_p, name = resultsNames(dds)[2])
result1_padj = as.data.frame(result1[which(result1$padj < alpha_p),])

library(biomaRt)
cols <- c("ensembl_gene_id","external_gene_name",'description','ensembl_transcript_id','chromosome_name','start_position','end_position','strand','gene_biotype')
ensembl = useMart("ensembl", dataset='hsapiens_gene_ensembl')
rownames(result1_padj) = substring(rownames(result1_padj), 0, 15) #15 is a hardcode length of ENSG id
result1_padj['ensembl_gene_id'] = rownames(result1_padj)
biomart_response = distinct(select(ensembl, keys=rownames(result1_padj), columns=cols, keytype="ensembl_gene_id"), ensembl_gene_id, .keep_all = TRUE)
biomart_response = biomart_response %>% mutate('DEG_CATEGORY' = resultsNames(dds)[2]) 
biomart_response = merge(biomart_response, result1_padj[order(rownames(result1_padj)),], by="ensembl_gene_id", all = T)
write.csv(biomart_response, "A-M_DEG.csv")


#Adipose v Myocardium METABOLOMICS==============================================
#===============================================================================
library(caret)
adipose = read.csv(".\\omics\\datasets\\adipose_targeted.csv")
myocyte = read.csv(".\\omics\\datasets\\myocyte_targeted.csv")
myocyte = myocyte[myocyte$Original %in% adipose$Original, ]
adipose = adipose[adipose$Original %in% myocyte$Original, ]
myocyte = myocyte[,order(colnames(myocyte))]
adipose = adipose[,order(colnames(adipose))]

rownames(adipose) = adipose$Original
rownames(myocyte) = myocyte$Original
myocyte = subset(myocyte, select = -c(`Original`))
adipose = subset(adipose, select = -c(`Original`))
all(rownames(myocyte) %in% rownames(adipose)) #same rownames
all(rownames(myocyte) == rownames(adipose)) #same rownames order

metabo = cbind(adipose, myocyte) 
metabo_norm =  as.data.frame(normalizeBetweenArrays(metabo, method = "quantile")) #Makes the distributions highly uniform.

######################ASSESSING NORMALIZATION###################################
# View data ranges:
range(metabo)
range(metabo_norm)
# Boxplots - Visualize Distribution of each variable
boxplot(metabo, main = "Before Normalization")
boxplot(metabo_norm, main = "After Normalization")
# Density Plots - Visualize metabolite abundances (distribution of intensities)
plot(density(metabo[, 1]), main = "Density Plots Before Normalization")
for (i in 2:ncol(metabo)) {lines(density(metabo[, i]))}
plot(density(metabo_norm[, 1]), main = "Density Plots After Normalization")
for (i in 2:ncol(metabo_norm)) {lines(density(metabo_norm[, i]), col = i)}
# PCA - What do custering differences look like among samples?
pca_before <- prcomp(t(metabo))
plot(pca_before$x[, 1], pca_before$x[, 2], main = "PCA Before Normalization")
pca_after <- prcomp(t(metabo_norm))
plot(pca_after$x[, 1], pca_after$x[, 2], main = "PCA After Normalization")
#Depending on PCA, heatmaps, boxplot, anova: may need to use removeBatchEffect()
################################################################################

len = (length(metabo_norm)/2)
metabo_design = data.frame(category = c(replicate(len, 0), replicate(len, 1)))
metabo_design <- model.matrix(~category, data = metabo_design)

library(limma)
fit <- eBayes(lmFit(metabo_norm, design = metabo_design))
stats_df <- topTable(fit, number = nrow(metabo_norm)) %>% tibble::rownames_to_column("Metabolite")
write.csv(stats_df, "A-M_METABO.csv")

