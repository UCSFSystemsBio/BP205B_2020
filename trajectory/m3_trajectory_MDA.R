source('/wynton/home/students/snanda/rds/bp205/analysis/common_functions.R')
library(monocle3)

mda_cds_raw <- R_data_directory_to_cds('/wynton/scratch/bp205/processed/MDA_normalized_adata_R/')
###############################################################################
mda_celltype_classification<- data.table::fread('/wynton/scratch/bp205/signatures/BP205B_2020/MDA_celltype_classification.csv')

mda_cds_raw <- mda_cds_raw[,mda_celltype_classification$V1]
colData(mda_cds_raw)[['Subpopulation']] <- mda_celltype_classification$status
colData(mda_cds_raw)[['leiden']] <- as.factor(mda_celltype_classification$sorted_leiden)

mda_umap <- as.matrix(data.table::fread('/wynton/scratch/bp205/signatures/BP205B_2020/MDA_UMAP_coord.csv',header = F))
rownames(mda_umap) <- colnames(mda_cds_raw)

mda_cds <- preprocess_cds(mda_cds_raw,norm_method = 'none')
#plot_pc_variance_explained(mda_cds)
# mda_cds2 <- reduce_dimension(mda_cds,max_components = 2,verbose = T)
# m3_umap <- reducedDims(mda_cds2)[['UMAP']]  

reducedDims(mda_cds2)[['UMAP']]  <- mda_umap


#1059 works, 2519 equivalent
#8542 has good graphs for everyone except for primary parental
#rand_seed <- sample(1:10000,size=1)

mda_cds2 <- cluster_cells(mda_cds2,cluster_method = 'leiden',reduction_method = 'UMAP',random_seed = 1059)
# plot_cells(mda_cds2,
#            color_cells_by = 'cluster',
#            label_cell_groups = F,
#            label_leaves = T,
#            label_branch_points = T,
#            cell_size = 1.5) + scale_color_manual(values = div_pal)

mda_cds2@clusters[['UMAP']]$partitions <- as.character(clusters(mda_cds2))
mda_cds2@clusters[['UMAP']]$partitions <- as.factor(mda_cds2@clusters[['UMAP']]$partitions)

mda_cds3_mda <- learn_graph(mda_cds2,verbose = T,close_loop = F,learn_graph_control = list(maxiter=500))
# plot_cells(mda_cds3_mda,
#            color_cells_by = 'Subpopulation',
#            label_cell_groups = F,
#            label_leaves = T,
#            label_branch_points = T,
#            cell_size = 2) #+ scale_color_manual(values = sample(div_pal))

mda_cds3_mda_o <- order_cells(mda_cds3_mda)

mda_leiden_coloring <- c('#1f77b4','#ff7f0e','#3bd7ff','#d62728','#aa40fc','#8c564b',
                         '#e377c2','#b5bd61','#17becf','#aec7e8','#ffbb78',
                         '#98df8a','#ff9896','#c5b0d5','#c49c94','#f7b6d2','#dbdb8d')

mda_plt <- plot_cells(mda_cds3_mda_o,
                         color_cells_by = 'leiden',
                         label_cell_groups = F,
                         label_leaves = T,
                         label_branch_points = F,trajectory_graph_color = 'black',
                         cell_size = 0.5 ) +
  scale_color_manual(values = mda_leiden_coloring)+
  theme(legend.title = element_text(color='white'))


ggsave(plot=mda_plt,
       '/wynton/home/students/snanda/rds/bp205/analysis/trajectory/MDA_trajectory.pdf',width=6,height = 5)

ggsave(plot=mda_plt,
       '/wynton/home/students/snanda/rds/bp205/analysis/figures/figure_2/f2_MDA_trajectory.pdf',width=6,height = 5)


mda_principal_curve <- generate_principal_curve(mda_cds2)
ggsave(plot=mda_principal_curve,
       '/wynton/home/students/snanda/rds/bp205/analysis/trajectory/MDA_trajectory.svg',width=5,height = 5)

save(mda_cds2,file = '/wynton/home/students/snanda/rds/bp205/analysis/trajectory/MDA_normalized_cc_adata_PT_std.rds')






# 
# 
# ## visualizing pseudotime population change
# mda_pseudotime_mda <-as_tibble(colData(mda_cds3_mda_o)) 
# mda_pseudotime_mda$pt <- pseudotime(mda_cds3_mda_o)
# mda_pseudotime_mda$partition <- mda_cds3_mda_o@clusters$UMAP$partitions
# 
# 
# 
# mda_pseudotime_mda %>% ggplot(aes(x=pt,y=leiden)) + geom_point()
# 
# 
# # 
# # mda_pseudotime_mda %>% ggplot(aes(x=leiden,y=pt)) + geom_boxplot() #+ geom_jitter()
# # mda_pseudotime_m3 %>% ggplot(aes(x=leiden,y=pt)) + geom_boxplot() + geom_jitter()
# # #
# # 
# # dplyr::left_join(mda_pseudotime_mda %>% filter(partition!=4) %>%  group_by(leiden) %>% summarize(m=median(pt)),
# #                   mda_pseudotime_m3 %>% filter(partition!=1) %>% group_by(leiden) %>% summarize(m=median(pt)),'leiden') %>%
# #   ggplot(aes(x=m.x,y=m.y,color=leiden)) + geom_point() +
# #   geom_abline(slope = 1,intercept = 0) + scale_color_manual(values = sample(div_pal))
# #      
# #      
# #      
# 
# #################################################
# 
# 
# mda_pseudotime_mda_2 <- mda_pseudotime_mda %>%  mutate(bin=ntile(pt,50)) 
# mda_pseudotime_mda_2 %>% ggplot(aes(x=bin,fill=leiden)) + geom_bar()+ scale_fill_manual(values = sample(div_pal))
# 
# 
# #####
# 
# 
# mda_pseudotime_mda$node <- mda_cds3_m3_o@principal_graph_aux[['UMAP']]$pr_graph_cell_proj_closest_vertex[,1]
# 
# f <- mda_pseudotime_mda %>% group_by(node,leiden) %>% summarize(n=n()) %>% ungroup %>% group_by(node) %>% mutate(m=sum(n),p=n/m)
# 
# f %>% ggplot(aes(x=node,y=p,fill=leiden)) + geom_bar(stat = 'identity')
# 
# 
# mda_cds3_m3_o@principal_graph_aux[['UMAP']]$pseudotime
# 
#  
# ica_space_df <- t(mda_cds3_mda_o@principal_graph_aux[[reduction_method]]$dp_mst) %>% 
#   as.data.frame() %>% 
#   dplyr::select_(prin_graph_dim_1 = x, 
#                  prin_graph_dim_2 = y) %>% 
#   dplyr::mutate(sample_name = rownames(.), sample_state = rownames(.))
# 
# 
# ica_space_df %>% ggplot(aes(x=prin_graph_dim_1,y=prin_graph_dim_2 )) + geom_point() + geom_text(aes(label=sample_name))
# 
# 
# edge_df <- cds@principal_graph[[reduction_method]] %>% 
#   igraph::as_data_frame() %>% dplyr::select_(source = "from", target = "to") %>%
#   dplyr::left_join(ica_space_df %>% 
#                      dplyr::select_(source = "sample_name", source_prin_graph_dim_1 = "prin_graph_dim_1",
#                                     source_prin_graph_dim_2 = "prin_graph_dim_2"), 
#                    by = "source") %>% dplyr::left_join(ica_space_df %>%  dplyr::select_(target = "sample_name", target_prin_graph_dim_1 = "prin_graph_dim_1", target_prin_graph_dim_2 = "prin_graph_dim_2"),  by = "target")
# 
# 
# 
# 
# 
# 
# #   data.frame(pt=pseudotime(mda_cds2)) %>% rownames_to_column() %>% as_tibble()
# # mda_pseudotime$partition <- mda_cds2@clusters$UMAP$partitions
# # mda_pseudotime$leiden <- colData(mda_cds2)$leiden
# # mda_pseudotime <- mda_pseudotime %>% group_by(partition) %>% mutate(bin = ntile(pt,30)) %>% arrange(desc(leiden))
# # 
# # mda_pseudotime %>% ggplot(aes(pt,fill=leiden)) + geom_histogram(aes(y=stat(count/sum(count)))) + facet_wrap(~partition,scales = 'free',ncol = 1) + scale_fill_manual(values = div_pal)
# # 
# # 
# # # 
# # 
# 
# 
