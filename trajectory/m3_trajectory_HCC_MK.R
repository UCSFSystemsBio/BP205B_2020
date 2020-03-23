source('/wynton/home/students/snanda/rds/bp205/analysis/common_functions.R')
library(monocle3)

hcc_mk_cds_raw <- R_data_directory_to_cds('/wynton/scratch/bp205/processed/HCC_MK_normalized_cc_adata_R/')
###############################################################################
mk_celltype_classification<- data.table::fread('/wynton/scratch/bp205/signatures/BP205B_2020/HCC_MK_celltype_classification.csv')
colData(hcc_mk_cds_raw)[['Subpopulation']] <- mk_celltype_classification$status
colData(hcc_mk_cds_raw)[['leiden']] <- as.factor(mk_celltype_classification$leiden)

hcc_mk_umap <- as.matrix(data.table::fread('/wynton/scratch/bp205/signatures/BP205B_2020/HCC_MK_UMAP_coords.txt',header = F))
rownames(hcc_mk_umap) <- colnames(hcc_mk_cds_raw)

hcc_mk_cds <- preprocess_cds(hcc_mk_cds_raw,norm_method = 'none')
#plot_pc_variance_explained(hcc_mk_cds)
# hcc_mk_cds2 <- reduce_dimension(hcc_mk_cds,max_components = 2,verbose = T)
hcc_mk_cds2 <- hcc_mk_cds
reducedDims(hcc_mk_cds2)[['UMAP']]  <- hcc_mk_umap

hcc_mk_cds2 <- cluster_cells(hcc_mk_cds2,cluster_method = 'leiden',reduction_method = 'UMAP')
## insert hcc_mk_umap into the object
hcc_mk_cds2 <- learn_graph(hcc_mk_cds2,close_loop = F,verbose = T,learn_graph_control = list(maxiter=500))

## Set root node as Y_25
hcc_mk_cds2 <- order_cells(hcc_mk_cds2)

hcc_mk_plt <- plot_cells(hcc_mk_cds2,
           color_cells_by = 'leiden',
           label_cell_groups = F,
           label_leaves = T,
           label_branch_points = T,trajectory_graph_color = 'black',
           cell_size = 2) #+ scale_color_manual(values = sample(div_pal))
ggsave(plot=hcc_mk_plt,
       '/wynton/home/students/snanda/rds/bp205/analysis/trajectory/HCC_MK_trajectory.pdf',width=5,height = 5)


hcc_mk_principal_curve <- generate_principal_curve(hcc_mk_cds2)
ggsave(plot=hcc_mk_principal_curve,
       '/wynton/home/students/snanda/rds/bp205/analysis/trajectory/HCC_MK_trajectory.svg',width=5,height = 5)

save(hcc_mk_cds2,file = '/wynton/home/students/snanda/rds/bp205/analysis/trajectory/HCC_MK_normalized_cc_adata_PT_std.rds')


# #############################################################################################################################
# ## Learn marker genes of Leiden clsuters
# hcc_mk_leidin_markers <- top_markers(hcc_mk_cds2,group_cells_by = 'status',verbose = T,cores = 18,reference_cells = 1000)
# 
# hcc_mk_leidin_markers_top <- hcc_mk_leidin_markers %>%
#   filter(fraction_expressing >=0.10) %>%
#   group_by(cell_group) %>% 
#   top_n(1,pseudo_R2)
# 
# 
# plot_cells(hcc_mk_cds2,genes =hcc_mk_leidin_markers_top$gene_id)
# 
# # plot_genes_by_group(hcc_mk_cds2,
# #                     unique(hcc_mk_leidin_markers$gene_id),
# #                     group_cells_by="leiden",
# #                     ordering_type="maximal_on_diag",
# #                     max.size=3)
# 
# 
# ### Order the cells in psuedotime
# hcc_mk_cds2 <- order_cells(hcc_mk_cds2)
# 
# plot_cells(hcc_mk_cds2,
#            color_cells_by = 'pseudotime',
#            label_cell_groups = F,
#            label_leaves = T,
#            label_branch_points = T,
#            cell_size = 0.5)+ylim(-8,5)
# 
# #### Graph test for gens varying w/ pseudotime
# hcc_mk_pt_genes <- graph_test(hcc_mk_cds2, neighbor_graph="principal_graph", cores=20)
# hcc_mk_pt_genes_signif <- hcc_mk_pt_genes %>% filter(q_value <0.01) %>% arrange(q_value,desc(morans_I))
# 
# hcc_mk_pt_subset <- hcc_mk_cds2[hcc_mk_pt_genes_signif$gene_short_name,]
# 
# # 
# # 
# # hcc_mk_pt_modules <- find_gene_modules(hcc_mk_pt_subset, resolution=1e-2)
# # 
# # 
# 
# ########
# cell_group_df <- tibble::tibble(cell=row.names(colData(hcc_mk_cds2)), 
#                                 cell_group=partitions(hcc_mk_cds2)[colnames(hcc_mk_cds2)])
# 
# agg_mat <- aggregate_gene_expression(neurons_cds, gene_module_df, cell_group_df)
# row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
# colnames(agg_mat) <- stringr::str_c("Partition ", colnames(agg_mat))
# 
# pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
#                    scale="column", clustering_method="ward.D2",
#                    fontsize=6)
# 
# 
# ### visualizing in pseudotime
# plot_genes_in_pseudotime(hcc_mk_pt_subset[13],color_cells_by = 'status',min_expr = 0.5,)#+scale_y_continuous()
# 
# 
# 
# 

