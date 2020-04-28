### RVTEST GWAS RESULTS ###

# Create merged GWAS results file combining all the chromosomes
#cat *.SingleWald.assoc | grep -v 'N_INFORMATIVE' > allChrs_FILE.assoc

#R
library(dplyr)
library(data.table)
library(qqman)


data <- fread("allChrs_FILE.assoc", header = F)

data <- data %>%
	rename(CHROM = "V1",
		POS = "V2",
		REF = "V3",
		ALT = "V4",
		N_INFORMATIVE = "V5",
		Test = "V6",
		Beta = "V7",
		SE = "V8",
		Pvalue = "V9")

##Calculate lambda

#Get p values
p <- data$Pvalue

#Get number of tests
n <- length(data$Pvalue)

#Calculate lambda
x2obs <- qchisq(as.numeric(as.character(p)), 1, lower.tail = FALSE)
x2exp <- qchisq((1:n)/n, 1, lower.tail = FALSE)
lambda <- median(x2obs)/median(x2exp)
lambda


##Annotate with closest genes (just for hits with p < 0.0001)

#Select just SNPs with P < 0.0001
merged_E4 <- data %>% filter(Pvalue<=0.0001)

#Read in gene coordinates table 
gene.coords <- read.table("NCBI38.gene.loc.txt")

### Loop over each line in results file

res.annotated=as.data.frame(do.call(rbind,lapply(1:nrow(merged_E4),function(x){
  
  ### Get SNP chr and bp and corresponding genes
  snp.chr=merged_E4[x,1]
  snp.bp=merged_E4[x,2]
  gene.coords.chr=gene.coords[which(gene.coords$V2==snp.chr),]
  
  ### Calculate distance between snp bp and all genes, select gene with minimum value, combine columns and print distance between SNP and gene
  ### If distance is 0 then SNP is within gene coordinates. Distance is in BP.
  b=cbind(merged_E4[x,],gene.coords.chr[which.min(abs(snp.bp-((gene.coords.chr$V3+gene.coords.chr$V4)/2))),])
  if(b$POS<b$V3){
    d=cbind(b,as.character(b$V3-b$POS))
  } else if(b$V4<b$POS){
    d=cbind(b,as.character(b$POS-b$V4))
  } else if(b$POS>b$V3 & b$POS<b$V4){
    d=cbind(b,as.character("0"))
  }
  names(d)[16]=c("Distance.To.Gene(BP)")
  d
})))

res.annotated=res.annotated[order(res.annotated$Pvalue),]

res.annotated <- res.annotated %>%
	rename(gene.chr = V2,
			gene.start = V3,
			gene.end = V4,
			gene.name = V6) %>%
	select(-V1)
	
write.table(res.annotated, "res.annotated.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)
q()
n