# PLM-R
This is a github repository for sample codes used to generate results for the following paper : " ".

## Environment
The codes were tested in R4.4 with following R packages : getopt_1.20.4, data.table_1.14.10

## General description
genotype_subset_wrapper.R : subset plink format genotype file to the file containing varinats in a given set of genes. \
assign_vep.R : convert genotype matrix to PLM score matrix for each gene. \
assign_vep_isoform.R : convert genotype matrix to PLM score matrix for each isoform (transcript). \
associate.R : perform linear regression/firth logistic regression for continuous/binary traits. 


## Parameters
```bash
./genotype_subset_wrapper.R
--chr : chromosome number
--gs_path : path for input gene set file (format for the gene set file can be seen in ./example/ENSG00000205560_gs.txt
--geno_path : path for input exome file (plink prefix)
--plink_path : path for plink executable
--out_path : output path 
```

```bash
./assign_vep.R
--gene_list : a text file containing gene names and paths for subsetted genotype files 
--id : genotype sample ID flie 
--vep_path : a list of variants and corresonding vep scores 
--out_path : output path

```

```bash
./assign_vep_isoform.R
--ts_list : a text file containing gene names, transcript names, and paths for subsetted genotype files
--id : genotype sample ID flie 
--vep_path : a list of variants and corresonding vep scores 
--out_path : output path
```

```bash
./associate.R
 --phen_file :  
--phen_type :
--cov_file :
--vep_file :
--out_file :
 
```

