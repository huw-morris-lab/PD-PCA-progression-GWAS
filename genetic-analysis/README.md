# Scripts for running GWAS:


1. create_pheno_file_script.txt

This reads in the progression scores created in R and makes the phenotype and people file inputs for rvtests


2. rvtests_GWAS_script.sh
  
This runs the GWAS in rvtests looping over each chromosome (calls the GWAS_per_chr.sh file).
  

3. plot_GWAS_rvtest_results_script

This calculates the lambda value for the GWAS and annotates the top results with the closest genes. No longer used for plotting as we used FUMA to create Manhattan and QQ plots.


4. add_RSID_toGWAS_new.R

This annotates GWAS results with rsIDs and creates the files to upload to FUMA.


5. cojo_script.txt

This runs COJO to generate the final list of independent SNPs, and merges with the nearest gene data to create final table.


6. PCs_vs_progression_script.txt

Sensibility check to make sure the clinical progression scores are not associated with any of the genetic principal components PC1-5.
