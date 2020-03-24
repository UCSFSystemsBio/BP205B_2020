source('/wynton/home/students/snanda/rds/bp205/analysis/common_functions.R')
library(VISION)
options(mc.cores=25)

## read in the data
pdx_exprmat <- h5read('/wynton/scratch/bp205/processed//PDX_adata.h5ad')
## Median normalization
norm.factor = median(colSums(pdx_exprmat))
pdx_exprmat <- t( t(pdx_exprmat) / colSums(pdx_exprmat)) * norm.factor

pdx_meta <- as.data.frame(data.table::fread('/wynton/scratch/bp205/processed/PDX_adata_R/obs.csv',sep = ','))
pdx_meta$SampleType <- factor(pdx_meta$SampleType,levels = paste0('Sample',1:26))
pdx_genes <- as.data.frame(data.table::fread('/wynton/scratch/bp205/processed/PDX_adata_R/var.csv',sep = ','))
rownames(pdx_genes) <- pdx_genes$V1
## Rename the exprmat
rownames(pdx_exprmat) <- rownames(pdx_genes)
colnames(pdx_exprmat) <- 1:ncol(pdx_exprmat)
###
pdx_umap <- h5read_umap('/wynton/scratch/bp205/processed/PDX_adata.h5ad')
rownames(pdx_umap) <- 1:ncol(pdx_exprmat)

##  Checkign normalized gene distributions

## Load and generate the signature files
hcc_mk_DE_results <- list.files('/wynton/scratch/bp205/signatures/BP205B_2020/DE/HCC_MK/FINAL/',
                                recursive = T,full.names = T,pattern = '.csv')
mda_DE_results <- list.files('/wynton/scratch/bp205/signatures/BP205B_2020/DE_results/Final/',
                             recursive = T,full.names = T,pattern = '.csv')
mda_DE_signatures <- create_signature_set('MDA',mda_DE_results,pdx_genes$V1,abs = F)
hcc_mk_DE_signatures <- create_signature_set('HCC_MK',hcc_mk_DE_results,pdx_genes$V1,abs = F)

all_DE_signatures <- c(mda_DE_signatures,hcc_mk_DE_signatures)
all_DE_signatures <- all_DE_signatures[sapply(all_DE_signatures,function(x){length(x@sigDict)})!=0]


## Load signature files by adding a sign to them 
pdx_vis <- VISION::Vision(data=pdx_exprmat,
                          signatures=all_DE_signatures,
                          sig_gene_threshold=0.005,
                          meta=pdx_meta,projection_methods=c('UMAP','tSNE30'),
                          pool=FALSE,name='PDX_agg_med_norm_all')

pdx_vis_result <- VISION::analyze(pdx_vis)
pdx_vis_result<- addProjection(pdx_vis_result,"UMAP2",pdx_umap)
save(pdx_vis_result,file = '/wynton/home/students/snanda/rds/bp205/analysis/pdx/PDX_agg_med_norm_all.rda',version = 2)

load('/wynton/home/students/snanda/rds/bp205/analysis/pdx/PDX_agg_med_norm_all.rda')

