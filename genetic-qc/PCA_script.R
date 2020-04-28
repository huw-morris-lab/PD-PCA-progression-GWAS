########## PCA FROM GENETIC DATA MERGED WITH HAPMAP ########## 
#Using only CEU HapMap samples
library(tidyverse)

#---Load principal components generated from GCTA---####
#This is the PCA from the PROBAND data merged with HapMap
#Using linkage-pruned SNPs, MAF > 5%, excluding palindromic SNPs and flipping missnps

#Read in eigenvectors
PCA <- as_tibble(read.table("PCA.eigenvec", sep = ""))

#Change column names
PCA <- PCA %>% 
  rename(FID = V1,
                IID = V2,
                PC1 = V3,
                PC2 = V4,
                PC3 = V5,
                PC4 = V6,
                PC5 = V7,
                PC6 = V8,
                PC7 = V9,
                PC8 = V10,
                PC9 = V11,
                PC10 = V12)


#---Read in data about HapMap samples (ethnicity)---####

#Read in population data for HapMap samples
HapMap_pops <- as_tibble(read.table("../reference/hapmap/relationships_w_pops_121708.txt",
                                    header = TRUE))

#---Merge PROBAND PCA data with HapMap population data---####

PCA$FID <- as.factor(PCA$FID)

#Join with PCA data table
PCA <- PCA %>% 
  left_join(HapMap_pops)

PCA <- PCA %>% 
 mutate(group = ifelse(!is.na(population), population, "PPMI"))

PCA <- PCA %>% 
  mutate(group = ifelse(group == 2, "CEU",
                        ifelse(group == "PPMI", "PPMI", NA)))

#---Plot first two Principal Components with HapMap samples---####

ggplot(data = PCA, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 0.9, alpha = 0.7) +
  scale_color_discrete(breaks=c("PPMI","CEU")) +
  theme_bw() +
  ggsave("PCA_PPMI_HapMap_CEUsamples.png")

#---Classify outliers who are > 6 SDs away from mean for any of the first 10 PCs---####

#Individuals who were more than 6 standard deviations away from the mean of the any of the first 10 principal components were removed.
#Each PC mean is calculated, and any individual who is away from the mean on any of the PCs is removed

#Look at the mean for each PC
PCA.PCmeans <- PCA %>% 
  dplyr::mutate(mean_PC1 = mean(PC1),
                mean_PC2 = mean(PC2),
                mean_PC3 = mean(PC3),
                mean_PC4 = mean(PC4),
                mean_PC5 = mean(PC5),
                mean_PC6 = mean(PC6),
                mean_PC7 = mean(PC7),
                mean_PC8 = mean(PC8),
                mean_PC9 = mean(PC9),
                mean_PC10 = mean(PC10))

#Remove individuals who are outliers on any PC

PC.outlierResults <- as.data.frame(matrix(ncol = 11, nrow = nrow(PCA)))
PC.outlierResults[, 1] <- PCA$IID


#For loop to calculate the SDs of each Principal Component (first 10 PCs only)
#This outputs into a results table
for (i in 3:12) {
  mean <- mean(PCA.PCmeans[[i]])
  sd <- sd(PCA.PCmeans[[i]])
  
  PC.outlierResults[, i-1] <- PCA.PCmeans %>% 
    mutate(outlier = ifelse(PCA.PCmeans[[i]] > mean + 6*sd, "outlier",
                            ifelse(PCA.PCmeans[[i]] < mean - 6*sd, "outlier", "keep"))) %>% 
    dplyr::select(outlier)
}

#Rename column names
PC.outlierResults <- PC.outlierResults %>% 
  dplyr::rename(ID = V1,
                PC1_result = V2,
                PC2_result = V3,
                PC3_result = V4,
                PC4_result = V5,
                PC5_result = V6,
                PC6_result = V7,
                PC7_result = V8,
                PC8_result = V9,
                PC9_result = V10,
                PC10_result = V11)

#Now merge the outlier results with the main dataset
PCA.PCmeans <- PCA.PCmeans %>% 
  left_join(PC.outlierResults, by = c("IID" = "ID"))

#If any of the PC results are outliers, flag as outlier
PCA.PCmeans <- PCA.PCmeans %>% 
  mutate(PCA_outlier = ifelse(PC1_result == "outlier" |
                                PC2_result == "outlier" |
                                PC3_result == "outlier" |
                                PC4_result == "outlier" |
                                PC5_result == "outlier" |
                                PC6_result == "outlier" |
                                PC7_result == "outlier" |
                                PC8_result == "outlier" |
                                PC9_result == "outlier" |
                                PC10_result == "outlier", "outlier final", "keep final"))

#Plot first 2 PCs by outlier status
ggplot(data = PCA.PCmeans, mapping = aes(x = PC1, y = PC2, color = PCA_outlier)) +
  geom_point(alpha = 0.2) +
  theme_bw()

#Plot third and fourth PC by outlier status
ggplot(data = PCA.PCmeans, mapping = aes(x = PC3, y = PC4, color = PCA_outlier)) +
  geom_point(alpha = 0.2) +
  theme_bw()

#Count how many outliers
PCA.PCmeans %>% 
  group_by(PCA_outlier) %>% 
  filter(group == "PPMI") %>% 
  dplyr::summarise(count = n())

#---Write list of samples that are population outliers---####
#To remove in PLINK

PCA_outliers <- PCA.PCmeans %>% 
  filter(PCA_outlier == "outlier final") %>% 
  filter(group == "PPMI") %>% 
  select(FID, IID)

#Write text file of FID and IID
write.table(PCA_outliers, "PCA_outliers.txt",
            row.names = FALSE, quote = FALSE, col.names = FALSE)
