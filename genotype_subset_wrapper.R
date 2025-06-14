library(data.table)
library(getopt)

####### Subset plink geno file by genes: LoF + missense

optspec = matrix(c(
    'chr', 'c', 1, 'integer', # chromosome name
    'gs_path','s', 1, 'character', # path for input gene set file
    'geno_path', 'g', 1, 'character', # path for input exome file (plink prefix)
    'plink_path', 'p', 1, 'character',
    'out_path', 'o', 1, 'character'), # output folder path for variant/geno subset
     byrow=TRUE, ncol=4)

opt = getopt(optspec)
chr = as.integer(opt$chr)
gs_path = as.character(opt$gs_path)
geno_path = as.character(opt$geno_path)
plink_path = as.character(opt$plink_path)
out_path = as.character(opt$out_path)


## Comment this out and provide your own files
## gs_path = "./example/ENSG00000205560_gs.txt"
## plink_path = "/u/project/zaitlenlab/sunkjang/bin/plink2"
## geno_path =  "./example/example"
## out_path = "./example/output/"

gs = fread( gs_path )
gene_list = unique(gs$Gene )

print(chr)


for(i in 1:length(gene_list)){

    print(i)
    gene = gene_list[i]
    vars = data.table( unique(gs[Gene==gene, Var] ) )
    var_output = paste0(out_path, gene,"_vars.txt")
    geno_output = paste0(out_path, gene,"_subset.txt")
 
    fwrite(vars, file=var_output, sep="\t", col.names=FALSE, quote=FALSE)
    
    system(paste0(plink_path, " --bfile ", geno_path, " --extract ", var_output," --make-bfile --out ",geno_output)) 

    file.remove(var_output)
 
}

