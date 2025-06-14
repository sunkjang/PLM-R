library(data.table)
library(getopt)
library(bigstatsr)
library(SPAtest)
library(scales)

optspec = matrix(c(
    'phen_file', 'p', 1, 'character',
    'phen_type', 't',1,'character', # either qt or bt
    'cov_file', 'c',1,'character',
    'vep_file','v',1,'character',
    'out_file','o',1,'character'),  
    byrow=TRUE, ncol=4)

opt = getopt(optspec)
phen_file = opt$phen_file
cov_file = opt$cov_file
vep_file = opt$vep_file
out_file = opt$out_file
phen_type = opt$phen_type


### Comment this out and provide your own paths
## phen_file = "./example/phen_qt.txt"
## phen_file = "./example/phen_bt.txt"
## cov_file = "./example/cov.txt"
## vep_file = "./example/example_ukbb_assigned.txt.gz"
## out_file = "./example/example_reg.txt"

### Load files 
phen = fread(phen_file)
cov = fread(cov_file)
vep = fread(vep_file)

### Merge pheno, cov, ems1b scores 
dat = Reduce( function(x,y) merge(x,y,by="IID"), list(phen, cov, vep))

### Get gene list    
pred = colnames(vep)[-1]
cov_list = colnames(cov)[-1]

### Run lm using bigstatsr
if(phen_type=="qt"){

y= as.matrix( dat[,2]) # response
covar = dat[, ..cov_list]
tmp_f = tempfile()

X = FBM(dat[, ..pred], nrow=nrow(dat), ncol=length(pred), type="double", backingfile = tmp_f ) # conver mx to FBM

test = big_univLinReg(X,y, covar.train=as.matrix(covar))
test$p.value = predict(test, log10=FALSE)
colnames(test)=c("beta","SE","t","pval")
test$gene_id = pred    
test$LOG10P= -log(test$pval,10)
test = data.table(test)
    
# save
fwrite(test, file=out_file, row.names=FALSE, sep="\t", append=TRUE)
}    


if(phen_type=="bt"){

y= as.matrix( dat[, 2])
covar = dat[, ..cov_list]
tmp_f = tempfile()

X = FBM(dat[, ..pred], nrow=nrow(dat), ncol=length(pred), type="double", backingfile = tmp_f ) # conver mx to FBM


    
# Run logistic regression using fastR    
test = big_univLogReg(X,y, covar.train=as.matrix(covar))
test$p.value = predict(test, log10=FALSE)
colnames(test)=c("beta","SE","niter","z","pval")
test$gene_id = pred    
test$test_type = "logistic"
    
# Run Firth logistic regression for pvalue < 0.01
test = data.table(test)
test_sub = test[pval<0.01]
test_sub = test_sub[!is.na(niter)] # remove non-converged ones

# Fit null model
null_form = as.formula(paste0("pheno ~ ", paste0(c(cov_list), collapse=" + ")))
res_null = ScoreTest_wSaddleApprox_NULL_Model(null_form, data=dat)

# Run spa
    if( nrow(test_sub)>0){

       res_spa_all = NULL

       for(j in 1:nrow(test_sub)){

       print(paste0("SPA-",j))

       gene_id_f = test_sub$gene_id[j]
       dat[, "gene_rescaled":=rescale( get(gene_id_f), to=c(0,2))]
       res = ScoreTest_SPA( dat$gene_rescaled, dat[, ..phen_name], obj.null=res_null, method="SPA", beta.out=TRUE, beta.Cutoff=1 )

       res_spa = c(res$beta, res$SEbeta, res$Is.converge, NA, res$p.value, gene_id_f,  "SPA" ) # firth beta, firth SE, converge, SPA pval, score Pval, test_type
       res_spa = data.table(t(res_spa))
       res_spa_all= rbind(res_spa_all, res_spa)

   }
       colnames(res_spa_all) = colnames(test_sub)
       res_spa_all = data.table(res_spa_all)
       test_all_spa = data.table( rbind(test, res_spa_all) )
    }   else {
       test_all_spa = test
    }


fwrite(test_all_spa, file=out_file, row.names=FALSE, sep="\t")

}




        
     

























