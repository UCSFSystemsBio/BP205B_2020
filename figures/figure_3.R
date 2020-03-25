source('/wynton/home/students/snanda/rds/bp205/analysis/common_functions.R')
library(RColorBrewer)
outdir <- '/wynton/home/students/snanda/rds/bp205/analysis/figures/figure_3/'


## PDX visualization for various tSNES, as well as leiden clustering + barplot of composition
indiv_genes <- c('PTTG1', 'IFI27', 'DDIT4', 'PLK1', 'MALAT1', 'GAL',key_genes$key_genes)
## rel sig names

## Load the vision object
vision_results <- get(load('/wynton/home/students/snanda/rds/bp205/analysis/pdx/PDX_agg_med_norm_all.rda'))

pdx_signature_scores <- vision_results@SigScores %>%
  as.data.frame %>% rownames_to_column(var = 'cell') %>%
  select_at(vars(c('cell',rel_sig_names)))

pdx_expr_df <-t(vision_results@exprData) %>% as.data.frame %>% rownames_to_column(var = 'cell') 
pdx_expr_df <- pdx_expr_df %>%
                select_at(vars(c('cell',indiv_genes[indiv_genes %in% colnames(pdx_expr_df)])))

pdx_vision_meta <- vision_results@metaData
vision_tSNE <- as_tibble(vision_results@Projections$tSNE30)

vision_sig_df <- dplyr::bind_cols(pdx_signature_scores,pdx_vision_meta)
vision_gene_df <- dplyr::bind_cols(pdx_expr_df,pdx_vision_meta)


###### tSNES
## Rendering for tSNE
vision_tSNE_df <- dplyr::bind_cols(vision_tSNE,vision_sig_df)

f3_pdx_subpopulation_tSNE <- vision_tSNE_df %>%
  ggplot(aes(x=`tSNE30-1`,y=`tSNE30-2`,color=Subtype)) + geom_point(size=0.1) +
  scale_color_brewer(palette = 'Set1')+
  theme_classic()+
  labs(x='tSNE30 1',y='tSNE30 2',color='Subpopulation') +
  guides(colour = guide_legend(override.aes = list(size=3)))

ggsave(plot = f3_pdx_subpopulation_tSNE,filename = paste0(outdir,'/f3_pdx_subpopulation_tSNE.pdf'),width = 8,height = 6)

f3_pdx_cluster_tSNE <- vision_tSNE_df %>%
  ggplot(aes(x=`tSNE30-1`,y=`tSNE30-2`,color=VISION_Clusters)) + geom_point(size=0.1) +
  scale_color_manual(values = div_pal)+
  theme_classic()+
  labs(x='tSNE30 1',y='tSNE30 2',color='Leiden Cluster')+
  guides(colour = guide_legend(override.aes = list(size=3)))

ggsave(plot = f3_pdx_cluster_tSNE,filename = paste0(outdir,'/f3_pdx_cluster_tSNE.pdf'),width = 8,height = 6)

f3_pdx_subpopulation_barplot <- vision_tSNE_df %>% group_by(VISION_Clusters,Subtype) %>% summarize(n=n()) %>% mutate(freq = n/sum(n))%>%
  ggplot(aes(x=VISION_Clusters,y=freq,fill=Subtype)) + geom_bar(stat='identity',width = 0.5)+
  scale_fill_brewer(palette = 'Set1')+
  scale_y_continuous(expand = c(0,0))+
  expand_limits(y=1.03)+
  theme_classic() +
  labs(x='Leiden Cluster',y='Fraction',fill='Subpopulation')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  guides(colour = guide_legend(override.aes = list(size=3)))

ggsave(plot = f3_pdx_subpopulation_barplot,filename = paste0(outdir,'/f3_pdx_subpopulation_barplot.pdf'),width = 8,height = 6)


## tSNES
plot_vision <- function(vis_df,pdx_sig,outdir){
  print(pdx_sig)
  
 vision_plt <- vis_df %>% ggplot(aes(x=`tSNE30-1`,y=`tSNE30-2`,color=eval(parse(text=pdx_sig)))) + geom_point(size=0.1) +
    scale_color_viridis()+
    labs(color='Activity')+
    ggtitle(label = pdx_sig)+
    theme_classic()+
   guides(colour = guide_legend(override.aes = list(size=3)))
 
 ggsave(plot = vision_plt,filename = paste0(outdir,'/f3_pdx_VISION_',pdx_sig,'.pdf'),width=8,height = 6)
 
}

lapply(names(rel_sigs),plot_vision,vis_df = vision_tSNE_df,outdir=outdir)


## Marker plots


### Plotting individual gene expression
plot_gene_expression <- function(gene,melted_df,outdir){
  print(gene)
  melted_df <- melted_df %>% filter(variable==gene)
  
  gene_expression_plt_1 <- melted_df %>% 
    ggplot(aes(x=Subtype,y=value,fill=Subtype)) + 
    geom_violin() + facet_grid(~VISION_Clusters,scales = 'free_y')+
    theme_classic()+
    labs(x='Subpopulation',y='Expression')+
    theme(axis.text.x = element_text(angle = 45,hjust = 1))
  
  ggsave(plot = gene_expression_plt_1 ,
         filename = paste0(outdir,'/f3_gene_expression_',gene,'.pdf'),width = 12,height = 3)
  
  gene_expression_plt_1_logscale <- gene_expression_plt_1 + scale_y_log10()
  
  ggsave(plot = gene_expression_plt_1_logscale ,
         filename = paste0(outdir,'/f3_gene_expression_',gene,'_log10.pdf'),width = 12,height = 3)
  
  scaled_exp <- melted_df %>%
    group_by(VISION_Clusters,Subtype,variable) %>% 
    summarize(`Mean Expression` = mean(value),`Fraction Expressing`=sum(value>0)/n())
  
  marker_plt <- scaled_exp %>% 
    ggplot(aes(x=Subtype,y=VISION_Clusters,color=`Mean Expression`)) + 
    geom_point(aes(size=`Fraction Expressing`))+
    scale_color_gradient(low='white',high='red')+
    theme_classic() +
    labs(x='Subpopulation',y='Leiden Cluster')+
    theme(axis.text=element_text(color='black'))
  
  ggsave(plot = marker_plt ,
         filename = paste0(outdir,'/f3_gene_expression_',gene,'_markers.pdf'),width = 5,height = 6)

}

## Individual gene boxplots stratified by multiple features
vision_gene_df_melted <- vision_gene_df %>% reshape2::melt(id.vars=colnames(vision_gene_df)[!(colnames(vision_gene_df) %in% indiv_genes)])

lapply(colnames(vision_gene_df)[colnames(vision_gene_df) %in% indiv_genes],
       plot_gene_expression,melted_df=vision_gene_df_melted,outdir=outdir)




