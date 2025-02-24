
library(dplyr)
library(DESeq2)
setwd("C:\\Users\\User\\Desktop\\omics-analysis")
counts = read.csv("Full_Counts.csv")
metadata = read.csv("Validated_Metadata.csv")


#Step 1: SUBSETTING SAMPLES - CONSOLIDATING SAMPLES ACROSS METADATA COLUMNS & COUNTS
#Commenting out lines 13-16 will increase sample size
metadata = metadata[!(metadata$EAT.CATEGORY %in% ""),]
metadata = metadata[!(metadata$HEART.DISEASE %in% ""),]
#metadata = metadata[!(metadata$HIGH.BLOOD.PRESSURE %in% ""),]
#metadata = metadata[!(metadata$CAD %in% ""),]
ct = dim(metadata)[1]
ct2 = (ct+1)
ct3 = (ct2*2)

rownames(counts) = counts[,1]
colnames(counts) = substring(colnames(counts), 7, 10)
counts = counts[,c(2:81)]
counts = counts[, !colnames(counts) %in% c('D56M','D56A'), drop = F] #D25A D56M, D56A are variance outliers


#Step 2: BUILDING DESIGN LABEL (OUTCOME) FOR CATEGORIES FROM CONSOLIDATED METADATA SAMPLES
cat_EAT = (metadata[, c("DONOR.ID", "EAT.CATEGORY")])[order(metadata$DONOR.ID, decreasing = FALSE),]
colnames(cat_EAT) = c("sample","EATcat")
cat_EAT = rbind(cat_EAT, cat_EAT)
cat_EAT[0:ct,1] = paste(cat_EAT[0:ct,1],"A",sep="") #these values change depending on which are included
cat_EAT[ct2:ct3,1] = paste(cat_EAT[ct2:ct3,1],"M",sep="")
row.names(cat_EAT) = NULL
cat_EAT = cat_EAT[order(cat_EAT$sample,decreasing = FALSE), ]
cat_EAT = cat_EAT[cat_EAT$sample %in% colnames(counts), ]
row.names(cat_EAT) = NULL

cat_SEX = (metadata[, c("DONOR.ID", "SEX")])[order(metadata$DONOR.ID, decreasing = FALSE),]
colnames(cat_SEX) = c("sample","SEXcat")
cat_SEX = rbind(cat_SEX, cat_SEX)
cat_SEX[0:ct,1] = paste(cat_SEX[0:ct,1],"A",sep="")
cat_SEX[ct2:ct3,1] = paste(cat_SEX[ct2:ct3,1],"M",sep="")
row.names(cat_SEX) = NULL
cat_SEX = cat_SEX[order(cat_SEX$sample,decreasing = FALSE), ]
cat_SEX = cat_SEX[cat_SEX$sample %in% colnames(counts), ]
row.names(cat_SEX) = NULL

cat_HD = (metadata[, c("DONOR.ID", "HEART.DISEASE")])[order(metadata$DONOR.ID, decreasing = FALSE),]
colnames(cat_HD) = c("sample","HDcat")
cat_HD = rbind(cat_HD, cat_HD)
cat_HD[0:ct,1] = paste(cat_HD[0:ct,1],"A",sep="")
cat_HD[ct2:ct3,1] = paste(cat_HD[ct2:ct3,1],"M",sep="")
row.names(cat_HD) = NULL
cat_HD = cat_HD[order(cat_HD$sample,decreasing = FALSE), ]
cat_HD = cat_HD[cat_HD$sample %in% colnames(counts), ]
row.names(cat_HD) = NULL

cat_HBP = (metadata[, c("DONOR.ID", "HIGH.BLOOD.PRESSURE")])[order(metadata$DONOR.ID, decreasing = FALSE),]
colnames(cat_HBP) = c("sample","HBPcat")
cat_HBP = rbind(cat_HBP, cat_HBP)
cat_HBP[0:ct,1] = paste(cat_HBP[0:ct,1],"A",sep="")
cat_HBP[ct2:ct3,1] = paste(cat_HBP[ct2:ct3,1],"M",sep="")
row.names(cat_HBP) = NULL
cat_HBP = cat_HBP[order(cat_HBP$sample,decreasing = FALSE), ]
cat_HBP = cat_HBP[cat_HBP$sample %in% colnames(counts), ]
row.names(cat_HBP) = NULL

cat_CAD = (metadata[, c("DONOR.ID", "CAD")])[order(metadata$DONOR.ID, decreasing = FALSE),]
colnames(cat_CAD) = c("sample","CADcat")
cat_CAD = rbind(cat_CAD, cat_CAD)
cat_CAD[0:ct,1] = paste(cat_CAD[0:ct,1],"A",sep="")
cat_CAD[ct2:ct3,1] = paste(cat_CAD[ct2:ct3,1],"M",sep="")
row.names(cat_CAD) = NULL
cat_CAD = cat_CAD[order(cat_CAD$sample,decreasing = FALSE), ]
cat_CAD = cat_CAD[cat_CAD$sample %in% colnames(counts), ]
row.names(cat_CAD) = NULL

cat_AGE = (metadata[, c("DONOR.ID", "AGE")])[order(metadata$DONOR.ID, decreasing = FALSE),]
colnames(cat_AGE) = c("sample","AGEcat")
cat_AGE = rbind(cat_AGE, cat_AGE)
cat_AGE[0:ct,1] = paste(cat_AGE[0:ct,1],"A",sep="")
cat_AGE[ct2:ct3,1] = paste(cat_AGE[ct2:ct3,1],"M",sep="")
row.names(cat_AGE) = NULL
cat_AGE = cat_AGE[order(cat_AGE$sample,decreasing = FALSE), ]
cat_AGE = cat_AGE[cat_AGE$sample %in% colnames(counts), ]
row.names(cat_AGE) = NULL
cat_AGE$AGEcat[findInterval(cat_AGE$AGEcat, c(18,50)) == 1L] <- '0'
cat_AGE$AGEcat[findInterval(cat_AGE$AGEcat, c(50,60)) == 1L] <- '1'
cat_AGE$AGEcat[findInterval(cat_AGE$AGEcat, c(60,99)) == 1L] <- '2'


#Step 3: Confirm Counts, Metadata match in Samples and Sample Order
counts = select(counts, contains(unlist(cat_EAT[,1]))) #Subset counts samples to match METADATA
all(colnames(counts) %in% cat_EAT$sample) #identical
all(colnames(counts) == cat_EAT$sample)   #same order
dim(cat_EAT)
dim(counts)



#Differential Gene Expression Analysis =========================================
#===============================================================================
alpha_p = .05
dds = DESeqDataSetFromMatrix(countData = counts, 
                             colData = cat_EAT, 
                             design = ~EATcat)
dds = dds[rowSums(counts(dds))>1]
dds$EATcat = relevel(dds$EATcat, ref="NORMAL") #Set reference level manually
dds = DESeq(dds)
resultsNames(dds)

dds2 = DESeqDataSetFromMatrix(countData = counts, 
                              colData = cat_EAT, 
                              design = ~EATcat)
dds2 = dds2[rowSums(counts(dds2))>1]
dds2$EATcat = relevel(dds2$EATcat, ref="OVERWEIGHT")
dds2 = DESeq(dds2)
resultsNames(dds2)

result1 = results(dds, alpha = alpha_p, name = 'EATcat_OBESE_vs_NORMAL')
result1_padj = result1[which(result1$padj < alpha_p),] 
result2 = results(dds, alpha = alpha_p, name = 'EATcat_OVERWEIGHT_vs_NORMAL')
result2_padj = result2[which(result2$padj < alpha_p),] 
result3 = results(dds2, alpha = alpha_p, name = 'EATcat_OBESE_vs_OVERWEIGHT')
result3_padj = result3[which(result3$padj < alpha_p),]



#PCA, Dispersion, Boxplot, and MA-Plot==========================================
#===============================================================================
plotPCA(vst(dds, blind=TRUE), intgroup=c("EATcat", "EATcat"))

plotDispEsts(dds, main = "Dispersion Plot - Shrunk Estimates in Blue") 

boxplot(log10(assays(dds)[['cooks']]), range=0, las=2, main='Boxplot of Cooks Values')

plotMA(result1, cex = 0.7, ylim = c(-6,6), main='Most DEGs in OBESE v NORMAL')
abline(h=c(-1,1), col='red', lwd=3)



#biomaRt GENE LOOKUP ===========================================================
#===============================================================================
library(biomaRt)
cols <- c("ensembl_gene_id","external_gene_name",'description','ensembl_transcript_id','chromosome_name','start_position','end_position','strand','gene_biotype')
ensembl = useMart("ensembl", dataset='hsapiens_gene_ensembl')
#listDatasets(ensembl)
#listAttributes(ensembl)

obese = substring(rownames(result1_padj), 0, 15) #15 is a hardcode length of ENSG id
oweight = substring(rownames(result2_padj), 0, 15)
obow = substring(rownames(result3_padj), 0, 15)

results = distinct(select(ensembl, keys=obese, columns=cols, keytype="ensembl_gene_id"), ensembl_gene_id, .keep_all = TRUE)
results1 = results %>% mutate('DEG_CATEGORY' = 'OBESE v NORMAL') 
results = distinct(select(ensembl, keys=oweight, columns=cols, keytype="ensembl_gene_id"), ensembl_gene_id, .keep_all = TRUE)
results2 = results %>% mutate('DEG_CATEGORY' = 'OVERWEIGHT v NORMAL') 
results = distinct(select(ensembl, keys=obow, columns=cols, keytype="ensembl_gene_id"), ensembl_gene_id, .keep_all = TRUE)
results3 = results %>% mutate('DEG_CATEGORY' = 'OBESE v OVERWEIGHT') 
results = rbind(results1, results2, results3)

write.csv(results, "C:\\Users\\User\\Desktop\\EAT.csv")



#METABOLOMICS DIFFERENTIAL =====================================================
#===============================================================================
library(limma)
library(caret)
adipose = read.csv("adipose_targeted.csv")
myocyte = read.csv("myocyte_targeted.csv")
myocyte = myocyte[myocyte$Original %in% adipose$Original, ]
adipose = adipose[adipose$Original %in% myocyte$Original, ]
myocyte = myocyte[,order(colnames(myocyte))]
adipose = adipose[,order(colnames(adipose))]
rownames(adipose) = adipose$Original
rownames(myocyte) = myocyte$Original
rownames(myocyte) == rownames(adipose)

myocyte = subset(myocyte, select = -c(`Original`))
adipose = subset(adipose, select = -c(`Original`))
metabo = cbind(adipose, myocyte) #d04, D11, D36, D56, D57 
metabo = metabo[,order(colnames(metabo))]

cat_MET_EAT = cat_EAT[cat_EAT$sample %in% c(colnames(metabo)),]
cat_MET_SEX = cat_SEX[cat_SEX$sample %in% c(colnames(metabo)),]
cat_MET_HD = cat_HD[cat_HD$sample %in% c(colnames(metabo)),]
cat_MET_HBP = cat_HBP[cat_HBP$sample %in% c(colnames(metabo)),]
cat_MET_CAD = cat_CAD[cat_CAD$sample %in% c(colnames(metabo)),]
cat_MET_AGE = cat_AGE[cat_AGE$sample %in% c(colnames(metabo)),]

met_EAT_ob = log(subset(metabo, select = names(metabo) %in% cat_MET_EAT[cat_MET_EAT$EATcat == 'OBESE',]$sample))
met_EAT_ow = log(subset(metabo, select = names(metabo) %in% cat_MET_EAT[cat_MET_EAT$EATcat == 'OVERWEIGHT',]$sample))
met_EAT_nl = log(subset(metabo, select = names(metabo) %in% cat_MET_EAT[cat_MET_EAT$EATcat == 'NORMAL',]$sample))
met_SEX_m = log(subset(metabo, select = names(metabo) %in% cat_MET_SEX[cat_MET_SEX$SEXcat == 'Male',]$sample))
met_SEX_f = log(subset(metabo, select = names(metabo) %in% cat_MET_SEX[cat_MET_SEX$SEXcat == 'Female',]$sample))
met_HD_el = log(subset(metabo, select = names(metabo) %in% cat_MET_HD[cat_MET_HD$HDcat == 'Electrical',]$sample))
met_HD_st = log(subset(metabo, select = names(metabo) %in% cat_MET_HD[cat_MET_HD$HDcat == 'Structural',]$sample))
met_HD_nl = log(subset(metabo, select = names(metabo) %in% cat_MET_HD[cat_MET_HD$HDcat == 'Normal',]$sample))
met_HBP_y = log(subset(metabo, select = names(metabo) %in% cat_MET_HBP[cat_MET_HBP$HBPcat == 'Y',]$sample))
met_HBP_n = log(subset(metabo, select = names(metabo) %in% cat_MET_HBP[cat_MET_HBP$HBPcat == 'N',]$sample))
met_CAD_y = log(subset(metabo, select = names(metabo) %in% cat_MET_CAD[cat_MET_CAD$CADcat == 'Y',]$sample))
met_CAD_n = log(subset(metabo, select = names(metabo) %in% cat_MET_CAD[cat_MET_CAD$CADcat == 'N',]$sample))
met_AGE_0 = log(subset(metabo, select = names(metabo) %in% cat_MET_AGE[cat_MET_AGE$AGEcat == '0',]$sample))
met_AGE_1 = log(subset(metabo, select = names(metabo) %in% cat_MET_AGE[cat_MET_AGE$AGEcat == '1',]$sample))
met_AGE_2 = log(subset(metabo, select = names(metabo) %in% cat_MET_AGE[cat_MET_AGE$AGEcat == '2',]$sample))

lb_EAT_ob = replicate(22, 0)
lb_EAT_ow = replicate(16, 1)
lb_EAT_nl = replicate(24, 2)
lb_SEX_m = replicate(26, 0)
lb_SEX_f = replicate(36, 1)
lb_HD_el = replicate(18, 0)
lb_HD_st = replicate(8, 1)
lb_HD_nl = replicate(36, 2)
lb_HBP_y = replicate(30, 0)
lb_HBP_n = replicate(32, 1)
lb_CAD_y = replicate(10, 0)
lb_CAD_n = replicate(52, 1)
lb_AGE_0 = replicate(20, 0)
lb_AGE_1 = replicate(24, 1)
lb_AGE_2 = replicate(18, 2)


#Perform Differential Metabolomics Analysis: Pairwise using above values========
metabo_label = data.frame(category = c(lb_EAT_ow, lb_EAT_nl))
differential = cbind(met_EAT_ow, met_EAT_nl)
differential = as.data.frame(as.matrix(differential))

des_mat <- model.matrix(~category, data = metabo_label)
fit <- eBayes(lmFit(differential, design = des_mat))
stats_df <- topTable(fit, number = nrow(differential)) %>% tibble::rownames_to_column("Metabolite")
#topTable(fit,sort.by="P")
write.csv(stats_df, "C:\\Users\\User\\Desktop\\EAT_ownl.csv")



volcanoplot(fit,coef=1, highlight=50, names=rownames(scaled), main='Volcano', pch=10, cex=0.25)
plotMDS(scaled, col = as.numeric(labels$category))

