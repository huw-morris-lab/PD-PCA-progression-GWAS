Scripts for running GWAS:


1. create_pheno_file_script.txt

  This reads in the progression scores created in R and makes the phenotype and people file inputs for rvtests

2. rvtests_GWAS_script.sh
  
  This runs the GWAS in rvtests looping over each chromosome (calls the GWAS_per_chr.sh file).
  
3. plot_GWAS_rvtest_results_script

  This calculates the lambda value for the GWAS and annotates the top results with the closest genes. No longer used for plotting as we used FUMA to create Manhattan and QQ plots.

4. add_RSID_toGWAS_new.R

  This annotates GWAS results with rsIDs and creates the files to upload 
