#Differential metabolites were identified using non-parametric statistical methods (Wilcox's rank test). 
#Two filtration criteria (FC > 1.2 and p < 0.05) were used to identify differential metabolites. 
#I just read some paper where they removed outliers by PCA "based on a tolerance of 95% beyond the Hotelling's T2 ellipse"
#install.packages("ggrepel")
#BiocManager::install("limma", force = TRUE)
#devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, build_vignettes = FALSE)
options(max.print=1000)
library(dplyr)
library(plotly)
library(tidyverse)
setwd("C:\\Users\\User\\Desktop\\omics-analysis\\")
adipose = read.csv("adipose_targeted.csv")
myocyte = read.csv("myocyte_targeted.csv")
EAT <- read.csv("EAT.csv")
dim(myocyte)
dim(adipose)
rownames(myocyte) == rownames(adipose)
myocyte <- myocyte[-(391:421),]
rownames(myocyte) == rownames(adipose)

myocyte_cols = colnames(myocyte)
myocyte_cols = sort(myocyte_cols,decreasing=TRUE)
adipose_cols = colnames(adipose)
adipose_cols[1] = 'Original'
colnames(adipose)[1] = 'Original'
adipose_cols = sort(adipose_cols,decreasing=TRUE)

adipose = adipose[, adipose_cols]
myocyte = myocyte[, myocyte_cols]
myocyte_rows = rownames(myocyte)
adipose_rows = rownames(adipose)
differences = myocyte_rows != adipose_rows 
differences

colnames(adipose) = gsub('.{1}$','',colnames(adipose))
colnames(myocyte) = gsub('.{1}$','',colnames(myocyte))
colnames(adipose)[1] = 'Original'
colnames(myocyte)[1] = 'Original'

OB = EAT[EAT$EAT.CATEGORY == 'OBESE',]$DONOR.ID
OW = EAT[EAT$EAT.CATEGORY == 'OVERWEIGHT',]$DONOR.ID
NL = EAT[EAT$EAT.CATEGORY == 'NORMAL',]$DONOR.ID

OB_adip = subset(adipose, select = names(adipose) %in% OB)
OW_adip = subset(adipose, select = names(adipose) %in% OW)
NL_adip = subset(adipose, select = names(adipose) %in% NL)
OB_myoc = subset(myocyte, select = names(myocyte) %in% OB)
OW_myoc = subset(myocyte, select = names(myocyte) %in% OW)
NL_myoc = subset(myocyte, select = names(myocyte) %in% NL)
colnames(OB_adip) = lapply(colnames(OB_adip), paste0, "A")
colnames(OW_adip) = lapply(colnames(OW_adip), paste0, "A")
colnames(NL_adip) = lapply(colnames(NL_adip), paste0, "A")
colnames(OB_myoc) = lapply(colnames(OB_myoc), paste0, "M")
colnames(OW_myoc) = lapply(colnames(OW_myoc), paste0, "M")
colnames(NL_myoc) = lapply(colnames(NL_myoc), paste0, "M")
OB = cbind(OB_adip, OB_myoc)
OW = cbind(OW_adip, OW_myoc)
NL = cbind(NL_adip, NL_myoc)
rownames(OB) = adipose[,1]
rownames(OW) = adipose[,1]
rownames(NL) = adipose[,1]

write.csv(OB, "C:\\Users\\User\\Desktop\\Efimov\\metabolomic\\Obese.csv")
write.csv(OW, "C:\\Users\\User\\Desktop\\Efimov\\metabolomic\\Overweight.csv")
write.csv(NL, "C:\\Users\\User\\Desktop\\Efimov\\metabolomic\\Normal.csv")

##############################################################
library(dplyr)
library(plotly)
library(tidyverse)
setwd("C:\\Users\\User\\Desktop\\")
OB <- read.csv("Obese.csv")
OW <- read.csv("Overweight.csv")
NL <- read.csv("Normal.csv")
rownames(OB) = OB[,1]
rownames(OW) = OW[,1]
rownames(NL) = NL[,1]
rownames(OB) = rownames(NL)
rownames(OW) = rownames(NL)
OB = OB[,-1]
OW = OW[,-1]
NL = NL[,-1]
colSums(is.na(NL))
colSums(is.na(OB))
colSums(is.na(OW))

process <- preProcess(as.data.frame(NL), method='range')
NL_p <- predict(process, as.data.frame(NL))
process <- preProcess(as.data.frame(OB), method='range')
OB_p <- preProcess(as.data.frame(OB), method='range')
process <- preProcess(as.data.frame(OW), method='range')
OW_p <- preProcess(as.data.frame(OW), method='range')
OB_L = log(OB)
OW_L = log(OW)
NL_L = log(NL)

##############################################################
#library(devtools)
library(limma)
library(caret)
#dat <- read.csv("metabodat.csv")
#dat2 = as.data.frame(dat) 
#dat2 = sapply(dat2, as.numeric)
lb_ob = replicate(34, 0)
lb_ow = replicate(24, 1)
lb_nm = replicate(30, 2)
label = data.frame(c(lb_ob, lb_ow))
colnames(label) = 'category'
colnames(OB) = lapply(colnames(OB), paste0, "ob")
colnames(OW) = lapply(colnames(OW), paste0, "ow")
OBOW = cbind(OB, OW)
OBOW = as.data.frame(as.matrix(OBOW))


des_mat <- model.matrix(~category, data = label)
fit <- lmFit(OBOW, design = des_mat)
fit <- eBayes(fit)

stats_df <- topTable(fit, number = nrow(OBOW)) %>%
  tibble::rownames_to_column("Metabolite")
stats_df

top = topTable(fit,sort.by="P")
top

volcanoplot(fit,coef=1, highlight=50, names=rownames(scaled), main='Volcano', pch=10, cex=0.25)
plotMDS(scaled, col = as.numeric(labels$category))

write.csv(stats_df, "C:\\Users\\User\\Desktop\\Efimov\\metabo_norm.csv")









