library(monocle3)
library(data.table)

R_data_directory_to_cds <- function(R_data_dir_path){
  sample_sheet <- as.data.frame(fread(paste0(R_data_dir_path,'/obs.csv')))
  colnames(sample_sheet)[1] <- 'sampleNames'
  rownames(sample_sheet) <- sample_sheet$sampleNames
  
  gene_annotation <- as.data.frame(fread(paste0(R_data_dir_path,'/var.csv')))
  colnames(gene_annotation)[1:2] <- c('gene_short_name','featureNames')
  rownames(gene_annotation) <- gene_annotation$gene_short_name
  
  expr <- t(as.matrix(fread(paste0(R_data_dir_path,str_sub(basename(R_data_dir_path),1,-3),'.csv'),header = F)))
  colnames(expr) <- sample_sheet$sampleNames
  rownames(expr) <- gene_annotation$gene_short_name
  
  cds <- new_cell_data_set(expr,sample_sheet,gene_annotation)
  return(cds)
  
}

hcc_mk_cds <- R_data_directory_to_cds('/wynton/scratch/bp205/processed/HCC_MK_adata_R/')
mk_celltype_classification<- data.table::fread('/wynton/scratch/bp205/signatures/BP205B_2020/HCC_MK_celltype_classification.csv')
colData(hcc_mk_cds)[['status']] <- mk_celltype_classification$status
colData(hcc_mk_cds)[['leiden']] <- mk_celltype_classification$leiden

hcc_mk_cds <- preprocess_cds(hcc_mk_cds)

#plot_pc_variance_explained(hcc_mk_cds)

hcc_mk_cds2 <- reduce_dimension(hcc_mk_cds,max_components = 2)
hcc_mk_cds2 <- cluster_cells(hcc_mk_cds2,cluster_method = 'leiden',reduction_method = 'UMAP')
hcc_mk_cds2 <- learn_graph(hcc_mk_cds2,verbose = T)

plot_cells(hcc_mk_cds2,
           color_cells_by = 'status',
           label_cell_groups = F,
           label_leaves = T,
           label_branch_points = T,
           cell_size = 0.5)+ylim(-8,5)
        

## Learn marker genes of Leiden clsuters

hcc_mk_leidin_markers <- top_markers(hcc_mk_cds2,group_cells_by = 'status',verbose = T,cores = 18,reference_cells = 1000)

hcc_mk_leidin_markers_top <- hcc_mk_leidin_markers %>%
  filter(fraction_expressing >=0.10) %>%
  group_by(cell_group) %>% 
  top_n(1,pseudo_R2)


plot_cells(hcc_mk_cds2,genes =hcc_mk_leidin_markers_top$gene_id)


# plot_genes_by_group(hcc_mk_cds2,
#                     unique(hcc_mk_leidin_markers$gene_id),
#                     group_cells_by="leiden",
#                     ordering_type="maximal_on_diag",
#                     max.size=3)


### Order the cells in psuedotime
hcc_mk_cds2 <- order_cells(hcc_mk_cds2)

plot_cells(hcc_mk_cds2,
           color_cells_by = 'pseudotime',
           label_cell_groups = F,
           label_leaves = T,
           label_branch_points = T,
           cell_size = 0.5)+ylim(-8,5)

#### Graph test for gens varying w/ pseudotime

hcc_mk_pt_genes <- graph_test(hcc_mk_cds2, neighbor_graph="principal_graph", cores=18)
hcc_mk_pt_genes_signif <- hcc_mk_pt_genes %>% filter(q_value <0.01) %>% arrange(q_value,desc(morans_I))


hcc_mk_pt_subset <- hcc_mk_cds2[hcc_mk_pt_genes_signif$gene_short_name,]



hcc_mk_pt_modules <- find_gene_modules(hcc_mk_pt_subset, resolution=1e-2)




########
cell_group_df <- tibble::tibble(cell=row.names(colData(hcc_mk_cds2)), 
                                cell_group=partitions(hcc_mk_cds2)[colnames(hcc_mk_cds2)])

agg_mat <- aggregate_gene_expression(neurons_cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("Partition ", colnames(agg_mat))

pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=6)


### visualizing in pseudotime
plot_genes_in_pseudotime(hcc_mk_pt_subset[1:5],color_cells_by = 'status',min_expr = 0.5,)+scale_y_continuous()





