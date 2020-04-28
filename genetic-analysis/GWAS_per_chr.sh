#!/bin/bash
	PHENO=$1
	CHNUM=$2
	KEEPFILE=$3
	
	/data/kronos/mtan/software/rvtest/executable/rvtest \
	--noweb \
	--hide-covar \
	--out $PHENO.$KEEPFILE.chr$CHNUM \
	--single wald \
	--inVcf /data/kronos/mtan/PROBAND_Oxford_PPMI_final/PROBAND_Oxford_PPMI_merged.snpqc.chr$CHNUM.vcf \
	--pheno $PHENO.pheno.txt \
	--pheno-name pheno \
	--covar $PHENO.pheno.txt \
	--covar-name cohort_PROBAND,cohort_Oxford,cohort_PPMI,PC1,PC2,PC3,PC4,PC5 \
	--peopleIncludeFile $KEEPFILE.txt \
	--qtl
