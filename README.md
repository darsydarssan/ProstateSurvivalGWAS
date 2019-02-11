# ProstateSurvivalGWAS
GWAS study of survival in European ancestry patients with prostate cancer

ProgramName	 - Purpose

genedatsplit.r - Read the SNP data and split into 54 data sets.

datamanip1.r	 - Join the prognestic data, PC data to 1 SNP data. Get rid of the SNPs which has more than 90% missing data. Encoding SNPs as per allele frequency. Getrid of SNPs which has same phenotype. Export the analysis ready data.

datamanip2.r	 - Merge a few data sets produceed by 3.1, snps with constant phenotype and more than 90% missing, these are in one file now.

univariatecox1.r - Fitting the unadjusted model and exported the overall p-value for each SNP

univariatecox2.r - Find the singnifican SNPs  unadjusted

adjustecox1.r	- Adjusted model fit, ajusting for age and 7 PCs

adjustecox2.r	- Find the singnifican SNPs of adjusted models

manhatton1.r	- Manhatton plot for unadjusted p-values

manhatton2.r	 - Manhatton plot for adjusted p-values

kmplot.r	- Kaplan-Mier Plot and Univariate analysis
