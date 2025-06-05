library(getopt)
library(bigsnpr)
library(data.table)
library(dplyr)


optspec = matrix(c(
    'gene_list', 'g', 1, 'character', 
    'id',        'i', 1, 'character', 
    'vep_path',  'v', 1, 'character', 
    'out_path',  'o', 1, 'character'),
    byrow=TRUE, ncol=4)

opt = getopt( optspec)
ts_list = as.character(opt$ts_list) # gene id, genotype path
vep_path = as.character(opt$vep_path) # Var, Gene, Feature, score
out_path = as.character(opt$out_path)
sample.id = as.character(opt$id) # fam file


### Comment this out and provide your own paths
vep.path = paste0("./example/example_score.txt") 
ts_list = paste0("./example/ts_list")
sample.id = paste0("./example/example.fam")
out.path = paste0("./example/example_isoform_ukbb_assigned.txt.gz")

### Functions
assign_vep = function(x){
    x[is.na(x)]=0
    res = max( vep.sub.ord$score[which(x>0)] ) # maximum; higher score = higher pathogenicity
    return(res) # if one carries more than 1 missense variant
    }


### Load VEP scores: make sure that higher scores mean higher pathogenicity
vep = fread( vep.path ) 
col.sel = c("Var","Gene","Feature","score")
vep = vep[, ..col.sel]
vep = vep[!is.na(score)]
vep = vep[score!="."]
vep[,"score":=-as.numeric( score )] # flip sign: high score = higher pathogenicity
    

### Load gene list : genotype file should have the same individuals
ts_list = fread(ts_list,header=FALSE) 

### Load sample ID 
sample.id = fread( sample.id )
sample.id = sample.id[,"V2"]
colnames(sample.id)="IID"

### Loop starts
vep.assigned.all = sample.id

for( i in 1:nrow(ts_list)){

    print(i)

    # Load gene subset
    gene_id = ts_list$V1[i]
    tmpfile = tempfile()
    rds = snp_readBed(ts_list$V3, backingfile = tmpfile)
    geno = snp_attach(rds)
    
    # Subset vep 
    ts_id = ts_list$V2
    vep.sub = vep[Gene==gene_id &Feature== ts_id ]

    # Subset overlapping vars
    rvs = geno$map$marker.ID
    var.ind = which(rvs %in% vep.sub$Var)

    if( length(var.ind)>0){  # run if more than 0 variant is left
        
    rds2 = snp_subset( geno, ind.col = var.ind)
    geno.sub = snp_attach(rds2)
    sample.id = data.table( geno.sub$fam$sample.ID )
    colnames(sample.id)="IID"

    rm(geno)
    rm(rds)
    gc()

    # Order vep variants
    rvs2 = geno.sub$map$marker.ID
    vep.sub.ord = vep.sub[match(rvs2, Var)]
    
    # Identify carriers 
    carr = apply( as.matrix(geno.sub$genotypes[]), 1, function(x) sum(x, na.rm=TRUE)  )
    carr.ind = which(carr>0)
    carr.ids = sample.id$IID[carr.ind]

    # Assign VEP scores for carriers     
    vep.assigned = data.table( apply( as.matrix( geno.sub$genotypes[carr.ind,] ), 1, function(x) assign_vep(x) ) )
    vep.assigned$IID = carr.ids

    # bind
    gene_ts_id = paste0(gene_id,"_",ts_id)        
    vep.assigned2 = merge(sample.id[,"IID"], vep.assigned, by="IID", all.x=TRUE)
    colnames(vep.assigned2)[2]= gene_ts_id
    vep.assigned2[is.na(get(gene_ts_id)), paste(gene_ts_id):=0]    
    vep.assigned.all = merge(vep.assigned.all, vep.assigned2, by="IID", all.x=TRUE)
            
    } else {

        print("no variants left")

    }
}

fwrite( vep.assigned.all, file=out.path, row.names=FALSE,quote=FALSE, sep="\t")

