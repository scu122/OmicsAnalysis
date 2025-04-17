#install.packages(c("dplyr", "caret"))
#BiocManager::install(c("DESeq2", "biomaRt", "BSgenome.Hsapiens.UCSC.hg38"))

library(dplyr)
setwd("C:\\Users\\User\\Desktop\\")
adipose = read.csv(".\\omics\\lipid\\adipose669.csv")
myocyte = read.csv(".\\omics\\lipid\\myocyte642.csv")
metadata = read.csv(".\\omics\\datasets\\Validated_Metadata.csv")


#Step 1: Subset Counts samples to Metadata samples, Truncate row names for consistency
metadata = metadata[!(metadata$EAT.CATEGORY %in% ""),]        #Reduce to available EAT
metadata = metadata[!(metadata$HEART.DISEASE %in% ""),]       #Reduce to available Heart Disease


#Step 2: Merge Adipose and Myocyte into 'Counts', Take only common lipidomic rows, 
names(adipose)[names(adipose) == "Original.ID"] = "Original"
counts = merge(adipose, myocyte, by = "Original")
counts = counts[,order(colnames(counts))]
row.names(counts) = counts$Original
row.names(metadata) = NULL


#Step 3: Prepare to create Design Matrix for categorical differences
obs_count = dim(metadata)[1]
labels = c("DONOR.ID", "EAT.CATEGORY") #"HEART.DISEASE" "HIGH.BLOOD.PRESSURE" "CAD" "SEX" "AGE"

#This function outputs a design matrix, automatically identifies: sample + 'category'
build_design_matrix <- function(row_count, vector_labels, metadata_df) {
  cat_ = (metadata_df[, c(vector_labels[1], vector_labels[2])])[order(metadata_df$DONOR.ID, decreasing = FALSE),]
  colnames(cat_) = c("sample",vector_labels[2])
  cat_ = rbind(cat_, cat_)
  cat_[0:row_count,1] = paste(cat_[0:row_count,1],"A",sep="") 
  cat_[(row_count+1):(row_count*2),1] = paste(cat_[(row_count+1):(row_count*2),1],"M",sep="")
  cat_ = cat_[order(cat_$sample,decreasing = FALSE), ]
  cat_ = cat_[cat_$sample %in% colnames(counts), ]
  row.names(cat_) = NULL
  colnames(cat_)[2] = "category"
  return(cat_)}


#Step 4: Create Design Matrix, subset counts samples to design matrix
counts_design = build_design_matrix(obs_count, labels, metadata) #Create the design matrix from the intended metadata category (Label)
counts = select(counts, contains(unlist(counts_design[,1])))     #Subset Counts samples to match design matrix samples
counts_design$tissue = as.factor(substring(colnames(counts), 4, 4))


#Step 5: Verify counts and design matrix contain identical sample names and order of samples
all(colnames(counts) %in% counts_design$sample) #identical sample names
all(colnames(counts) == counts_design$sample)   #same order of sample names
dim(counts_design)
dim(counts)
any(counts < 0)

#Check the datatype of 'counts' and convert to numeric
str(counts)
counts = counts %>% mutate_at(colnames(counts), as.numeric)
rm(adipose, myocyte) #remove these data objects

#This is a good time to review the counts_design object, as it gets transformed into factors later



#Visualize Variance and SD of Normalized Lipidomic Counts=======================
#===============================================================================
# Calculate the average abundance for each lipid (row-wise mean)
average_abundance = rowMeans(counts)
standard_deviation = apply(counts, 1, sd)
variance = apply(counts, 1, var)

plot(average_abundance, standard_deviation,
     xlab = "Average Lipid Abundance",
     ylab = "Standard Deviation",
     main = "Standard Deviation vs. Average Abundance")
abline(lm(standard_deviation ~ average_abundance), col = "red")

plot(average_abundance, variance,
     xlab = "Average Lipid Abundance",
     ylab = "Variance",
     main = "Variance vs. Average Abundance")
abline(lm(variance ~ average_abundance), col = "blue") 
#If you see the pattern where increasing average abundance has increasingly higher variance, 
#then limma-trend will likely be beneficial for this specific contrast. 

which.max(variance)
which.max(average_abundance)
variance[12]



#Differential Lipidomics Analysis - 1 ==========================================
#===============================================================================
#PAIRWISE ANALYSIS of all samples controlled, and contrast for EAT category
library(limma)
library(caret)
#Step 6: Create factor levels for analysis, Use these factor levels to build a 'model.matrix'
reference = factor(counts_design$category) #tissue  (if doing 'A' v 'M')
levels(reference)

counts_design = model.matrix(~0 + reference, data = counts_design)
colnames(counts_design) = levels(reference)


#Step 7: Build Limma Objects
fit = lmFit(counts, counts_design)
fit_ebayes = eBayes(fit, trend = TRUE) #This is Limma-Trend

# Define the contrasts for the three pairwise comparisons
contrast_matrix = makeContrasts(OBESE_vs_NORMAL = OBESE - NORMAL,
                                OVERWEIGHT_vs_NORMAL = OVERWEIGHT - NORMAL,
                                OVERWEIGHT_vs_OBESE = OVERWEIGHT - OBESE,
                                levels = counts_design)

# Apply the contrasts to the fit object
fit_contrasts = contrasts.fit(fit_ebayes, contrast_matrix)
eb_contrasts  = eBayes(fit_contrasts, trend = TRUE)


#Step 7: Produce the Differential Expression P-Values
# 1. Obese vs. Normal
top_obese_vs_normal = topTable(eb_contrasts, coef = "OBESE_vs_NORMAL", number = Inf) %>% tibble::rownames_to_column("Lipid")
top_obese_vs_normal$CATEGORY = paste(levels(reference)[2],"v", levels(reference)[1])
head(top_obese_vs_normal)
# 2. Overweight vs. Normal
top_overweight_vs_normal = topTable(eb_contrasts, coef = "OVERWEIGHT_vs_NORMAL", number = Inf) %>% tibble::rownames_to_column("Lipid")
top_overweight_vs_normal$CATEGORY = paste(levels(reference)[3],"v", levels(reference)[1])
head(top_overweight_vs_normal)
# 3. Overweight vs. Obese
top_overweight_vs_obese = topTable(eb_contrasts, coef = "OVERWEIGHT_vs_OBESE", number = Inf) %>% tibble::rownames_to_column("Lipid")
top_overweight_vs_obese$CATEGORY = paste(levels(reference)[3],"v", levels(reference)[2])
head(top_overweight_vs_obese)


#Step 8: Visualize the distribution of P-Values by category
hist(top_obese_vs_normal$P.Value)
hist(top_overweight_vs_normal$P.Value)
hist(top_overweight_vs_obese$P.Value)

plot(ecdf(top_obese_vs_normal$P.Value))
plot(ecdf(top_overweight_vs_normal$P.Value))
plot(ecdf(top_overweight_vs_obese$P.Value))


#Step 0: Save P-Values
write.csv(top_obese_vs_normal, "lipid_OBESE_vs_NORMAL.csv")
write.csv(top_overweight_vs_normal, "lipid_OVERWEIGHT_vs_NORMAL.csv")
write.csv(top_overweight_vs_obese, "lipid_OVERWEIGHT_vs_OBESE.csv")








