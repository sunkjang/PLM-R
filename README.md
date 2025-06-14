# PLM-R
This is a github repository for sample codes used to generate results for the following paper : " ".

## Environment
The codes were tested in R4.4 with following R packages : getopt_1.20.4, data.table_1.14.10

## General description
1) Subset genotype file to gene-specific variants \
genotype_subset_wrapper.R : subset plink format genotype file to the file containing varinats in a given set of genes.

2) Assign VEP scores for individuals \
assign_vep.R : convert genotype matrix to PLM score matrix for each gene. \
assign_vep_isoform.R : convert genotype matrix to PLM score matrix for each isoform (transcript). 

3) Perform regression \
associate.R : perform linear regression/firth logistic regression for continuous/binary traits. 


## Parameters
```bash
./genotype_subset_wrapper.R
--chr # chromosome number
--gs_path # path for input gene set file (format for the gene set file can be seen in ./example/ENSG00000205560_gs.txt
--geno_path # path for input exome file (plink prefix)
--plink_path # path for plink executable
--out_path # output path 
```

```bash
./assign_vep.R
--gene_list # a text file containing gene names and paths for subsetted genotype files in earlier step
--id # genotype sample ID flie 
--vep_path # a list of variants and corresonding vep scores 
--out_path # output path

```

```bash
./assign_vep_isoform.R
--ts_list # a text file containing gene names, transcript names, and paths for subsetted genotype files
--id # genotype sample ID flie
--vep_path # a list of variants and corresonding vep scores 
--out_path # output path
```

```bash
./associate.R
 --phen_file # a path for phenotype file  
--phen_type # phenotype type (either "qt" for continuous trait or "bt" for binary trait)
--cov_file # a path for covariate file
--vep_file # a list of variants and corresonding vep scores 
--out_file # output path
 
```

Please see the files in the ./example folder for examples of the file formats.

