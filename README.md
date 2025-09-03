# COVID-19 scATAC caQTL Manuscript

This repository contains code accompanying Mu et al., Cell Genomics (Accepted). The scripts are examples of the original analysis, with unnecessary duplicates removed.

## Scripts

The scripts can be viewed roughly in the following order:

1. `01_analyze_scATAC_clean.R`: Basic analysis on combined COVID-19, Benaglio et al., and You et al. clean data

2. `02_fit_fastTopics.R`: Applies the fastTopics algorithm to the scATAC-seq data
   - Implements fastTopics by iteratively adding more topics in fitting so results are comparable between fits with different topic numbers

3. `03_analyze_fastTopics.R`: Analysis and plotting of the fastTopics results

4. `04_glimpse_3study_donor.smk`: Snakemake workflow for genotype imputation with GLIMPSE

5. `05_L1_3study_covid19_poi.R`: Example script for sc-PME caQTL mapping

6. `06_coloc_TransAll_L1_metasoft.R`: Example script for colocalization between RA GWAS and caQTL

7. `07_scPME_interaction.R`: Mapping and visualization of dynamic caQTL
