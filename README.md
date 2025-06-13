# PLM-R
This is a github repository for sample codes used to generate results for the following paper : " ".

## Environment
The codes were tested in R4.4 with following R packages : getopt_1.20.4, data.table_1.14.10

## General description
genotype_subset_wrapper.R : subset plink format genotype file to each gene. 
assign_vep.R : convert genotype matrix to PLM score matrix for each gene
assign_vep_isoform.R : convert genotype matrix to PLM score matrix for each isoform (transcript).
associate.R : perform linear regression/firth logistic regression for continuous/binary traits


## Parameters

