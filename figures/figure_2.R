## Figure 2.R
source('/wynton/home/students/snanda/rds/bp205/analysis/common_functions.R')
library(pheatmap)
library(RColorBrewer)

select_and_score_metagene <- function(df,metagene){
  
  if(is.character(metagene)){
    genes <- metagene[metagene %in% colnames(df)]
    sign <- rep(1,length(genes))
    names(sign) <- genes
  }else{
    genes <- names(metagene)[names(metagene) %in% colnames(df)]
    sign <- sign(metagene[names(metagene) %in% colnames(df)])
  }
  scored <- df %>%
    select_at(vars(c('cell',genes))) %>%
    pivot_longer(all_of(genes)) %>%
    group_by_at(vars('cell')) %>%
    summarize(activation=sum(sign[name]*value)) %>% ungroup
  return(scored)
}

score_metagene_list <- function(df,meta,metagene_list){
  df$cell <- 1:nrow(df)
  scored <- map(metagene_list,select_and_score_metagene,df=df)
  
  scored_df <- scored %>% dplyr::bind_cols() %>% select_at(vars(c('cell','activation',paste0('activation',1:(length(metagene_list)-1)))))
  colnames(scored_df)[-1] <- names(metagene_list)
  
  scored_mat <- scored_df %>% dplyr::select(-cell) %>% as.matrix %>% t()
  colnames(scored_mat) <- meta$sampleNames
  
  scored_mat_scaled <-t(scale(t(scored_mat)))
  # scored_mat_scaled[scored_mat_scaled > 5] <- 5
  # scored_mat_scaled[scored_mat_scaled < -5] <- -5
  return(scored_mat_scaled)
}


## Load the data
hcc_mk_cds_raw <- R_data_directory_to_cds('/wynton/scratch/bp205/processed/HCC_MK_normalized_cc_adata_R/')
mk_celltype_classification<- data.table::fread('/wynton/scratch/bp205/signatures/BP205B_2020/HCC_MK_celltype_classification.csv')
colData(hcc_mk_cds_raw)[['Subpopulation']] <- mk_celltype_classification$status
colData(hcc_mk_cds_raw)[['leiden']] <- as.factor(mk_celltype_classification$leiden)

hcc_mk_expr <- t(as.matrix(exprs(hcc_mk_cds_raw)))
hcc_mk_expr_df <- as_tibble(hcc_mk_expr)
hcc_mk_meta <- as.data.frame(colData(hcc_mk_cds_raw))


mda_cds_raw <- R_data_directory_to_cds('/wynton/scratch/bp205/processed/MDA_normalized_adata_R/')
mda_celltype_classification<- data.table::fread('/wynton/scratch/bp205/signatures/BP205B_2020/MDA_celltype_classification.csv')
mda_cds_raw <- mda_cds_raw[,mda_celltype_classification$V1]
colData(mda_cds_raw)[['Subpopulation']] <- mda_celltype_classification$status
colData(mda_cds_raw)[['leiden']] <- as.factor(mda_celltype_classification$sorted_leiden)

mda_expr <- t(as.matrix(exprs(mda_cds_raw)))
mda_expr_df <- as_tibble(mda_expr)
mda_meta <- as.data.frame(colData(mda_cds_raw))

# Also need to load pseudotime CDS, can load via Save rather than reading in but need saves


## Load the signatures
hcc_mk_DE_results <- list.files('/wynton/scratch/bp205/signatures/BP205B_2020/DE/HCC_MK/FINAL/',
                                recursive = T,full.names = T,pattern = '.csv')
mda_DE_results <- list.files('/wynton/scratch/bp205/signatures/BP205B_2020/DE_results/Final/',
                             recursive = T,full.names = T,pattern = '.csv')

hcc_mk_sig_set <- load_signature_set('HCC_MK',hcc_mk_DE_results,colnames(hcc_mk_expr),abs = T,sign = F)
mda_sig_set <- load_signature_set('MDA',mda_DE_results,colnames(mda_cds_expr),abs = T,sign = F)

all_sig_set <- c(mda_sig_set,hcc_mk_sig_set)
all_sig_set <- all_sig_set[sapply(all_sig_set,length) != 0]


rel_sigs_mda_names <- c('MDA_0_vs_proxy_parental_UP_FC_1','MDA_1_vs_proxy_parental_UP_FC_1','MDA_3_vs_proxy_parental_UP_FC_1','MDA_4_vs_proxy_parental_UP_FC_1',
                   'MDA_11_vs_proxy_parental_UP_FC_1','MDA_11_vs_proxy_parental_DOWN_FC_1')
                   
rel_sigs_hcc_mk_names <- c('HCC_MK_0_vs_8_UP_FC_1','HCC_MK_2_vs_8_UP_FC_1','HCC_MK_4_vs_8_UP_FC_1','HCC_MK_7_vs_8_UP_FC_1')

rel_sigs_mda <- all_sig_set[rel_sigs_mda_names]
rel_sigs_hcc_mk <- all_sig_set[rel_sigs_hcc_mk_names]



##### Visualizations ##### ##### ##### ##### ##### ##### 
# Plot the score of each signature across the dataset
##### ##### ##### ##### ##### ##### ##### ##### ##### 
# MDA
mda_sigs_scored <- score_metagene_list(mda_expr_df,mda_meta,rel_sigs_mda)
mda_sigs_scored[mda_sigs_scored > 5] <- 5
mda_sigs_scored[mda_sigs_scored < -5] <- -5

f2_mda_sigs_scored_heatmap <- pheatmap(mda_sigs_scored[,order(mda_meta$leiden)],
         scale = 'none',
         color = rev(colorRampPalette(brewer.pal(8, "RdBu"))(25)),
         annotation_col = mda_meta[order(mda_meta$leiden),'leiden',drop=F],
         cluster_cols = F,
         cluster_rows = F,
         treeheight_row = 0,
         show_colnames = F,border_color = NA)
ggsave(plot=f2_mda_sigs_scored_heatmap,filename = '/wynton/home/students/snanda/rds/bp205/analysis/figures/figure_2/f2_mda_sigs_scored_heatmap.svg',width = 10,height = 3)
ggsave(plot=f2_mda_sigs_scored_heatmap,filename = '/wynton/home/students/snanda/rds/bp205/analysis/figures/figure_2/f2_mda_sigs_scored_heatmap.pdf',width = 10,height = 3)

## HCC
hcc_mk_sigs_scored <- score_metagene_list(hcc_mk_expr_df,hcc_mk_meta,rel_sigs_hcc_mk)

f2_hcc_mk_sigs_scored_heatmap <- pheatmap(hcc_mk_sigs_scored[,order(hcc_mk_meta$leiden)],
                                       scale = 'none',
                                       color = rev(colorRampPalette(brewer.pal(8, "RdBu"))(25)),
                                       annotation_col = hcc_mk_meta[order(hcc_mk_meta$leiden),'leiden',drop=F],
                                       cluster_cols = F,
                                       cluster_rows = F,
                                       treeheight_row = 0,
                                       show_colnames = F,border_color = NA)

ggsave(plot=f2_hcc_mk_sigs_scored_heatmap,filename = '/wynton/home/students/snanda/rds/bp205/analysis/figures/figure_2/f2_hcc_mk_sigs_scored_heatmap.svg',width = 10,height = 2)
ggsave(plot=f2_hcc_mk_sigs_scored_heatmap,filename = '/wynton/home/students/snanda/rds/bp205/analysis/figures/figure_2/f2_hcc_mk_sigs_scored_heatmap.pdf',width = 10,height = 2)


############ ##### ##### ##### 
## Plot the score of the top 10% most DE'd of each signature across the dataset
##### ##### ##### ##### ##### ##### 
# rel_sigs_mda_top20 <- map(rel_sigs_mda,function(x){
#   x[ntile(x,n = 10) %in% c(10,9)]
# })
# 
# 
# mda_sigs_scored_genes_top20 <- t(mda_expr[,unlist(unname(sapply(rel_sigs_mda_top20,names)))])
# mda_sigs_scored_genes_top20 <-t(scale(t(mda_sigs_scored_genes_top20)))
# 
# mda_sigs_scored_genes_top20 <- sign(mda_sigs_scored_genes_top20)*log10(abs(mda_sigs_scored_genes_top20))
# 
# 
# rownames(mda_sigs_scored_genes_top20) <- make.unique(rownames(mda_sigs_scored_genes_top20))
# 
# f2_mda_sigs_scored_genes_top20_heatmap <- pheatmap(mda_sigs_scored_genes_top20[,order(mda_meta$leiden)],
#                                              scale = 'none',
#                                              #color = colorRampPalette(brewer.pal(8, "PiYG"))(25),
#                                              annotation_col = mda_meta[order(mda_meta$leiden),'leiden',drop=F],
#                                              annotation_row =  data.frame(Signature=rep(names(rel_sigs_mda_top20),sapply(rel_sigs_mda_top20,length)),
#                                                                           row.names = rownames(mda_sigs_scored_genes_top20)),
#                                              cluster_cols = F,
#                                              cluster_rows = F,
#                                              treeheight_row = 0,
#                                              show_colnames = F,
#                                              show_rownames = F,
#                                              border_color = NA)

# ggsave(plot=f2_mda_sigs_scored_genes_top20_heatmap,filename = '/wynton/home/students/snanda/rds/bp205/analysis/figures/figure_2/f2_mda_sigs_scored_genes_top20_heatmap.svg',width = 20,height = 10)


###########################################################
## key genes visualization
###########################################################
key_genes <-list(key_genes=c('MT-ND4', 'NDUFB2', 'NDUFV2', 'MT-ND3', 'NDUFV1', 'MT-CO2', 'MT-CYB', 'FAM162A',
                             "HLA-DRA","HLA-DRB5","HLA-DRB1","HLA-DPA1","HLA-DPB1",'PTTG1',
                             'TRAC','CDC20', 'CCNB1', 'AURKA', 'RACGAP1',
                             'IFI27',
                             'DDIT4', 'NDRG1', 'HERPUD1',
                             'PSMD14', 'PSMA5',
                             'LGALS1',
                             'TUBA1B' ,
                             'BIRC5',
                             'H2AZ1',
                             'CKB', 'GAL', 'MALAT1'))

mda_key_genes <- colnames(mda_expr)[colnames(mda_expr) %in% key_genes$key_genes]

mda_key_genes_scored <- t(mda_expr[,mda_key_genes])
mda_key_genes_scored <- t(scale(t(mda_key_genes_scored)))
mda_key_genes_scored[mda_key_genes_scored >5 ] <- 5
mda_key_genes_scored[mda_key_genes_scored < -5] <- 5


f2_mda_key_genes_heatmap <- pheatmap(mda_key_genes_scored[,order(mda_meta$leiden)],
                                             scale = 'none',
                                              color = rev(colorRampPalette(brewer.pal(8, "RdBu"))(25)),
                                              annotation_col = mda_meta[order(mda_meta$leiden),'leiden',drop=F],
                                             cluster_cols = F,
                                             cluster_rows = T,
                                             treeheight_row = 0,
                                             show_colnames = F,
                                             show_rownames = T,
                                             border_color = NA)
ggsave(plot=f2_mda_key_genes_heatmap,filename = '/wynton/home/students/snanda/rds/bp205/analysis/figures/figure_2/f2_mda_key_genes_heatmap.svg',width = 8,height = 5)
ggsave(plot=f2_mda_key_genes_heatmap,filename = '/wynton/home/students/snanda/rds/bp205/analysis/figures/figure_2/f2_mda_key_genes_heatmap.pdf',width = 8,height = 5)



### HCC
hcc_mk_key_genes <- colnames(hcc_mk_expr)[colnames(hcc_mk_expr) %in% key_genes$key_genes]

hcc_mk_key_genes_scored <- t(hcc_mk_expr[,hcc_mk_key_genes])
hcc_mk_key_genes_scored <- t(scale(t(hcc_mk_key_genes_scored)))
hcc_mk_key_genes_scored[hcc_mk_key_genes_scored >5 ] <- 5
hcc_mk_key_genes_scored[hcc_mk_key_genes_scored < -5] <- 5


f2_hcc_mk_key_genes_heatmap <- pheatmap(hcc_mk_key_genes_scored[,order(hcc_mk_meta$leiden)],
                                     scale = 'none',
                                     color = rev(colorRampPalette(brewer.pal(8, "RdBu"))(25)),
                                     annotation_col = hcc_mk_meta[order(hcc_mk_meta$leiden),'leiden',drop=F],
                                     cluster_cols = F,
                                     cluster_rows = T,
                                     treeheight_row = 0,
                                     show_colnames = F,
                                     show_rownames = T,
                                     border_color = NA)
ggsave(plot=f2_hcc_mk_key_genes_heatmap,filename = '/wynton/home/students/snanda/rds/bp205/analysis/figures/figure_2/f2_hcc_mk_key_genes_heatmap.svg',width = 8,height = 5)
ggsave(plot=f2_hcc_mk_key_genes_heatmap,filename = '/wynton/home/students/snanda/rds/bp205/analysis/figures/figure_2/f2_hcc_mk_key_genes_heatmap.pdf',width = 8,height = 5)




###########################################################
## Common genes / Signatures fold change
###########################################################
common_genes_mda <- unique(unlist(unname(sapply(rel_sigs_mda,names))))
rel_sigs_mda_key_genes <-common_genes_mda[common_genes_mda %in% key_genes$key_genes]


common_genes_mat_mda <- matrix(0,nrow=length(common_genes_mda),ncol=length(rel_sigs_mda))
rownames(common_genes_mat_mda) <- common_genes_mda
colnames(common_genes_mat_mda) <- names(rel_sigs_mda)

for(i in 1:length(rel_sigs_mda)){
  common_genes_mat_mda[names(rel_sigs_mda[[i]]),names(rel_sigs_mda)[i]] <- rel_sigs_mda[[i]]
}
common_genes_mat_mda[common_genes_mat_mda > 10] <- 10
common_genes_mat_mda <- t(scale(t(common_genes_mat_mda)))
#rownames(common_genes_mat_mda)[!(common_genes_mda %in% rel_sigs_mda_key_genes)] <- '-'

common_genes_annot_mda <- data.frame(row.names = common_genes_mda)
common_genes_annot_mda$label <- ''
common_genes_annot_mda[rel_sigs_mda_key_genes,'label'] <- rel_sigs_mda_key_genes

f2_common_genes_signatures_mda_heatmap<- pheatmap(common_genes_mat_mda,
                                     scale = 'none',
                                     color = rev(colorRampPalette(brewer.pal(8, "RdBu"))(25)),
                                     cluster_cols = F,
                                     cluster_rows = T,
                                     annotation_row = common_genes_annot_mda,
                                     treeheight_row = 0,
                                     show_colnames = T,
                                     show_rownames = F,
                                     border_color = NA)


ggsave(plot=f2_common_genes_signatures_mda_heatmap,filename = '/wynton/home/students/snanda/rds/bp205/analysis/figures/figure_2/f2_common_genes_signatures_mda_heatmap.svg',width = 4,height = 11)
ggsave(plot=f2_common_genes_signatures_mda_heatmap,filename = '/wynton/home/students/snanda/rds/bp205/analysis/figures/figure_2/f2_common_genes_signatures_mda_heatmap.pdf',width = 4,height = 11)



### HCC
common_genes_hcc_mk <- unique(unlist(unname(sapply(rel_sigs_hcc_mk,names))))
rel_sigs_hcc_mk_key_genes <-common_genes_hcc_mk[common_genes_hcc_mk %in% key_genes$key_genes]


common_genes_mat_hcc_mk <- matrix(0,nrow=length(common_genes_hcc_mk),ncol=length(rel_sigs_hcc_mk))
rownames(common_genes_mat_hcc_mk) <- common_genes_hcc_mk
colnames(common_genes_mat_hcc_mk) <- names(rel_sigs_hcc_mk)

for(i in 1:length(rel_sigs_hcc_mk)){
  common_genes_mat_hcc_mk[names(rel_sigs_hcc_mk[[i]]),names(rel_sigs_hcc_mk)[i]] <- rel_sigs_hcc_mk[[i]]
}
common_genes_mat_hcc_mk[common_genes_mat_hcc_mk > 10] <- 10
common_genes_mat_hcc_mk <- t(scale(t(common_genes_mat_hcc_mk)))
#rownames(commn_genes_mat_hcc_mk)[!(common_genes_hcc_mk %in% rel_sigs_hcc_mk_key_genes)] <- '-'

common_genes_annot_hcc_mk <- data.frame(row.names = common_genes_hcc_mk)
common_genes_annot_hcc_mk$label <- ''
common_genes_annot_hcc_mk[rel_sigs_hcc_mk_key_genes,'label'] <- rel_sigs_hcc_mk_key_genes

f2_common_genes_signatures_hcc_mk_heatmap<- pheatmap(common_genes_mat_hcc_mk,
                                                  scale = 'none',
                                                  color = rev(colorRampPalette(brewer.pal(8, "RdBu"))(25)),
                                                  cluster_cols = F,
                                                  cluster_rows = T,
                                                  annotation_row = common_genes_annot_hcc_mk,
                                                  treeheight_row = 0,
                                                  show_colnames = T,
                                                  show_rownames = F,
                                                  border_color = NA)

ggsave(plot=f2_common_genes_signatures_hcc_mk_heatmap,filename = '/wynton/home/students/snanda/rds/bp205/analysis/figures/figure_2/f2_common_genes_signatures_hcc_mk_heatmap.svg',width = 3,height = 11)
ggsave(plot=f2_common_genes_signatures_hcc_mk_heatmap,filename = '/wynton/home/students/snanda/rds/bp205/analysis/figures/figure_2/f2_common_genes_signatures_hcc_mk_heatmap.pdf',width = 3,height = 11)





