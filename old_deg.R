#BiocManager::install(version = "3.20")
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
BiocManager::install("GRCh38.p14", force = TRUE)
#install.packages("plotly")

#LOAD DATA AND EXTRACTION=======================================================
#===============================================================================
library(dplyr)
library(DESeq2)
library(plotly)
setwd("C:\\Users\\User\\Desktop\\Efimov")
expression <- read.csv("Full_Counts.csv")
metadata <- read.csv("Validated_Metadata-EAT.csv")

#Remove NAN or Empty String in EAT
dim(metadata)
tail(metadata)
colSums(is.na(metadata))  #sum(is.na(metadata$EAT.CATEGORY)) 
apply(metadata, 2, function(x) sum(x == ""))

#Donor Cleanup, create outcome
metadata = metadata[!(metadata$EAT.CATEGORY %in% ""),]
outcome = metadata[, c("DONOR.ID", "EAT.CATEGORY")]
outcome = outcome[order(outcome$DONOR.ID, decreasing = FALSE),]
row.names(outcome) = NULL
colnames(outcome) = c("donor","outcome")
head(outcome)

#Counts Cleanup - Separate Samples, Harmonize Donor Names, Verify Sample Pairs
rownames(expression) <- expression[,1]
counts = expression[,c(2:81)]

colnames(counts) = substring(colnames(counts), 7, 10)
countsA = select(counts, contains('A'))
countsM = select(counts, contains('M'))
colnames(countsA) = substring(colnames(countsA), 0, 3)
colnames(countsM) = substring(colnames(countsM), 0, 3)
setdiff(colnames(countsM), colnames(countsA))

#Donor Counts Compared Cleanup
setdiff(unlist(outcome['donor']), colnames(countsA))
outcome = outcome[(unlist(outcome['donor']) %in% colnames(countsA)),]
countsA = select(countsA, contains(unlist(outcome['donor'])))
countsM = select(countsM, contains(unlist(outcome['donor'])))
dim(countsA)
dim(countsM)
dim(outcome)

#Recombine A+M DataFrame, Preserve EAT Outcome for A and M:
rownames(outcome) = unlist(outcome['donor'])
outcomeA = outcome
outcomeM = outcome
rownames(outcomeA) <- paste(rownames(outcome),"A",sep="") 
rownames(outcomeM) <- paste(rownames(outcome),"M",sep="") 
outcome = rbind(outcomeA, outcomeM)

colnames(countsA) <- paste(colnames(countsA),"A",sep="") 
colnames(countsM) <- paste(colnames(countsM),"M",sep="") 
counts = cbind(countsA, countsM)

dim(counts)
dim(outcome)
all(colnames(counts) %in% rownames(outcome)) #identical
all(colnames(counts) == rownames(outcome))   #same order

DEGcounts = DEGcounts[, !colnames(DEGcounts) %in% c('D56M','D56A'), drop = F] #D56M, D56A are variance outliers
DEGoutcomes = outcome[!rownames(outcome) %in% c('D56M','D56A'),, drop = F] #D56M, D56A are variance outliers

#Differential Gene Expression Analysis==========================================
#===============================================================================
#Vignette: https://bioconductor.org/packages/3.21/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential-expression-analysis
dds = DESeqDataSetFromMatrix(countData = counts, 
                             colData = outcome, 
                             design = ~outcome)
dds$outcome = relevel(dds$outcome, ref="NORMAL") #Set reference level manually

#keep = rowSums(counts(dds)) >= 30
#dds2 = dds[keep,]
al = .05

dds = DESeq(dds)
#collapseReplicates(dds) #To collapse technical replicates, not biological replicates
resultsNames(dds)
res1 = results(dds, alpha = al, name = 'Intercept')
res2 = results(dds, alpha = al, name = 'outcome_OBESE_vs_NORMAL')
res3 = results(dds, alpha = al, name = 'outcome_OVERWEIGHT_vs_NORMAL')
#res = results(dds, alpha = al, contrast = c('outcome', 'NORMAL','OVERWEIGHT'))
#res = results(dds, alpha = al, contrast = c(1,.5,.5))
#res <- lfcShrink(dds, coef=2)  #is this better for visualization?
#res <- lfcShrink(dds=dds, contrast=c("condition","B","A")) 

res2[which(res2$padj < .05 & res2$log2FoldChange < 0),]

summary(res1)
summary(res2)
summary(res3)
#res3[order(res3$pvalue),] #View ordered results
res2p = res2[which(res2$padj < .05),]
res2p
res3p = res3[which(res3$padj < .05),]
res3p$log2FoldChange


dds2 = DESeqDataSetFromMatrix(countData = counts, 
                              colData = outcome, 
                              design = ~outcome)
dds2$outcome = relevel(dds2$outcome, ref="OBESE") #Set reference level manually
dds2 = DESeq(dds2)
resultsNames(dds2)
res4 = results(dds2, alpha = al, name = 'outcome_OVERWEIGHT_vs_OBESE')
summary(res4)
res4p = res4[which(res4$padj < .05),]
res4p


#MA-Plot - #Most expressed genes in OBESE that are not in NORMAL
plotMA(res2, cex = 0.7, ylim = c(-6,6), main='Most DEGs in OBESE v Normal')
abline(h=c(-1,1), col='red', lwd=3)
plot

#MA-Plot - #Most expressed genes in OVERWEIGHT that are not in NORMAL
plotMA(res3, cex = 0.7, ylim = c(-6,6), main='Most DEGs in OVERWEIGHT v Normal')
abline(h=c(-1,1), col='red', lwd=3)
plot

#MA-Plot - #Most expressed genes in OVERWEIGHT that are not in OBESE
plotMA(res4, cex = 0.7, ylim = c(-6,6), main='Most DEGs in OVERWEIGHT v Obese')
abline(h=c(-1,1), col='red', lwd=3)
plot

#Dispersion
plotDispEsts(dds, main = "Dispersion Plot")

#Subset Counts to DEGs - Run PCA:
DEGcounts = counts[DEGs,] 
DEGcounts = DEGcounts[, !colnames(DEGcounts) %in% c('D56M','D56A'), drop = F] #D56M, D56A are variance outliers
DEGoutcomes = outcome[!rownames(outcome) %in% c('D56M','D56A'),, drop = F] #D56M, D56A are variance outliers
DEGdds = DESeqDataSetFromMatrix(countData = DEGcounts, 
                             colData = DEGoutcomes, 
                             design = ~outcome)

rld = rlogTransformation(DEGdds, blind = FALSE)
PCA <- plotPCA(rld, intgroup = c('outcome'))
fig <- plot_ly(data = PCA$data, x = PCA$data$PC1, y = PCA$data$PC2, color = PCA$data$group)
fig

#biomaRt GENE LOOKUP============================================================
#===============================================================================
#obese = data.frame(EnsembleGene = setdiff(rownames(res2p), rownames(res3p)))
#overweight = data.frame(EnsembleGene = setdiff(rownames(res3p), rownames(res2p)))
#DEGs = rbind(obese, overweight)
#normal = data.frame(EnsembleGene = setdiff(rownames(res4p),as.vector(DEGs['EnsembleGene'])))
#setdiff(as.vector(DEGs['EnsembleGene']),rownames(res4p))

obese = substring(rownames(res2p), 0, 15)
oweight = substring(rownames(res3p), 0, 15)
owob = substring(rownames(res4p), 0, 15)
#DEGs = DEGs %>% mutate('EAT' = 'OVERWEIGHT')
library(biomaRt)
cols <- c("ensembl_gene_id","external_gene_name",'description','ensembl_transcript_id','chromosome_name','start_position','end_position','strand','ucsc','gene_biotype')
ensembl = useMart("ensembl", dataset='hsapiens_gene_ensembl')
#listDatasets(ensembl)
#listAttributes(ensembl)
results = select(ensembl, keys=owob, columns=cols, keytype="ensembl_gene_id")
results = distinct(results, ensembl_gene_id, .keep_all = TRUE)
res3 = results %>% mutate('DEG_CATEGORY' = 'Overweight v Obese')

results = rbind(res1, res2, res3)
results 
data.frame(table(results$ensembl_gene_id))

write.csv(results, "C:\\Users\\User\\Desktop\\Efimov\\results.csv")

