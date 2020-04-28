### CONVERT BP:POSITIONS TO RSIDS FOR GWAS RESULTS FOR GRCH38###
#From Cornelis script for GRCh37
#Also using this page https://www.biostars.org/p/313408/
#Created 23/08/2019
#Last updated 10/01/2020
#WD: /data/kronos/mtan/
# download site-list for GRCh38

#wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz
#zcat 00-All.vcf.gz | cut -f 1,2,3,4,5 > GRCh38_rsids_new
#grep -v '#' GRCh38_rsids_new > GRCh38_rsids_new_build151
#Headers are #CHROM  POS     ID      REF     ALT


#Create chr:bp as a single column
#cut -f 1,2 GRCh38_rsids_new_build151 | sed 's/\t/:/'g > temp
#cut -f 3,4,5 GRCh38_rsids_new_build151 > temp2
#paste temp temp2 > GRCh38_rsids_new_build151_final.txt

#Separate multiple alleles into separate columns (this is not working in R as the file is too large)
cat GRCh38_rsids_new_build151_final.txt | cut -f 3,4 | sed 's/,/\t/'g > temp
cat GRCh38_rsids_new_build151_final.txt | cut -f 1,2 > temp2
paste temp2 temp > GRCh38_rsids_new_build151_final_sepalleles.txt

#Just select the first 3 alternate alleles (otherwise fread has issues with extra columns)
cat GRCh38_rsids_new_build151_final_sepalleles.txt | cut -f 1,2,3,4,5,6 > GRCh38_rsids_new_build151_final_sepalleles_final.txt

cd /data/kronos/mtan/PROBAND_Oxford_PPMI_final/PCA_GWAS_rvtests_PPMIbaseline_OFF_annualvisits

# Run the following step as a qsub script because it takes ages
#qsub -pe make 2 -cwd rsids_script.R

# merge in R with GWAS results
#R
#!/usr/bin/Rscript
library(dplyr)
library(tidyverse)
library(data.table)

#Read in GRCh38 rsID data
rsids <- fread("/data/kronos/mtan/GRCh38_rsids_new_build151_final_sepalleles_final.txt", header = F, fill = TRUE)

rsids <- rsids %>%
	rename(ALT_1 = "V4",
			ALT_2 = "V5",
			ALT_3 = "V6")
			
#Read in composite progression GWAS results
GWAS <- fread("allChrs_FILE.assoc", header=F)

#The rsIDs data has multiple rsIDs for one position based on the alleles
#Need to match by position as well as alleles
#Otherwise will have duplicates based on position which might not be the correct rsID - this will cause problems in FUMA

#Join GWAS results with rsIDs by position
merged <- GWAS %>%
	inner_join(rsids, by = c("V6" = "V1")) %>%
	rename(chr_bp = "V6",
			rsid = "V2.y",
			REF = "V3.y",
			chr = "V1",
			pos = "V2.x",
			REF_rvtests = "V3.x",
			ALT_rvtests = "V4",
			N = "V5",
			Beta = "V7",
			SE = "V8",
			Pvalue = "V9") %>%
	mutate(allele_match = ifelse(REF_rvtests == REF & ALT_rvtests == ALT_1, "match",
								ifelse(REF_rvtests == ALT_1 & ALT_rvtests == REF, "match",
								ifelse(REF_rvtests == REF & ALT_rvtests == ALT_2, "match",
								ifelse(REF_rvtests == ALT_2 & ALT_rvtests == REF, "match",
								ifelse(REF_rvtests == REF & ALT_rvtests == ALT_3, "match",
								ifelse(REF_rvtests == ALT_3 & ALT_rvtests == REF, "match",
								"remove")))))))	

#Filter out SNPs where the alleles do not match
merged_filtered <- merged %>%
	filter(allele_match == "match")
	
#Filter out indels as a lot of these have multiple rsIDs that are not in use
#Also the FUMA paper filters out indels so we are going to do it too
merged_filtered_noindels <- merged_filtered %>%
	filter(nchar(REF_rvtests) == 1) %>%
	filter(nchar(ALT_rvtests) == 1)


#Check for duplicates based on position
#Hopefully there are none as we have removed SNPs where the alleles do not match rsID alleles
merged_filtered_noindels$chr_bp[duplicated(merged_filtered_noindels$chr_bp)]

#Count number of SNPs in original merged file
merged %>%
	summarise(count = n())

#Count number of SNPs after filtering out SNPs with alleles that do not match rsID
merged_filtered_noindels %>%
	summarise(count = n())
#5434953 SNPs in final file

fwrite(merged_filtered_noindels, file="/data/kronos/mtan/PROBAND_Oxford_PPMI_final/PCA_GWAS_rvtests_PPMIbaseline_OFF_annualvisits/PCA_GWAS_sumstats_with_rsID_new.txt", col.names = TRUE, quote=FALSE,row.names=F,sep="\t")

rm(merged)
rm(merged_filtered_noindels)
rm(merged_filtered)
rm(GWAS)

q()
n

#Select columns for FUMA - need to have first column as a SNP identifier (chr_bp), so must use awk rather than cut
awk '{ print $6, $3, $4, $5, $7, $8, $9, $10 }' PCA_GWAS_sumstats_with_rsID_new.txt > PCA_GWAS_sumstats_with_rsID_new_FUMA.txt