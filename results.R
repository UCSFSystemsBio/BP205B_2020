## assemble a big data-set containing the signatures object, collated for each signature, and the survival object, and see what it is 
source('/wynton/home/students/snanda/rds/bp205/analysis/common_functions.R')

km_results <- fread('/wynton/home/students/snanda/rds/bp205/analysis/survival/results/survival_results.csv') %>%
  dplyr::filter(str_detect(signature,'FC_1'))

vision_results <- get(load('/wynton/home/students/snanda/rds/bp205/analysis/pdx/PDX_agg_med_norm_all.rda'))

geary <- vision_results@LocalAutocorrelation$Signatures %>% rownames_to_column()

diffs <- vision_results@ClusterComparisons$Signatures[c('Tissue','Subtype','VISION_Clusters')]
diffs_unwrapped <- do.call(cbind,lapply(diffs,function(d){
  do.call(cbind,d)#mapply(names(d),d,FUN=function(x,y){rownames(y) <- paste0(x,'_',rownames(y));y},SIMPLIFY = F))
})) %>% rownames_to_column()

vision_df <- left_join(geary,diffs_unwrapped,'rowname')

not_matched <- km_results$signature[!(km_results$signature %in% vision_df$rowname)]
all_results <- inner_join(km_results,vision_df,c('signature' = 'rowname'))

data.table::fwrite(all_results,'/wynton/home/students/snanda/rds/bp205/analysis/all_results.csv',sep = ',')


#####




all_results <- all_results %>% mutate_if(str_detect(colnames(all_results),'pval|pValue|FDR'),function(x){
  p <- -1*log10(x)
  p[is.infinite(p)] <- 300
  return(p)
  })


all_results <- all_results %>% filter(KM_pval > 3, KM_HZ!=1,pmin(n_high,n_low) >=150 )


all_results_data <- as.matrix(all_results[,-1])
rownames(all_results_data) <- all_results$signature

library(uwot)

all_results_umap <- uwot::umap(all_results_data,verbose = T,)
all_results_umap_df <- as_tibble(all_results_umap)
all_results_umap_df$signature <- all_results$signature


all_results_umap_df$cluster <- as.factor(kmeans(all_results_umap,3)$cluster)
all_results_umap_df %>% ggplot(aes(x=V1,y=V2,color=cluster)) + geom_point()

# pr <- prcomp(all_results_data)
# 

all_results_plt_df <- all_results_umap_df %>% left_join(all_results,'signature') %>% dplyr::select(-V1,-V2)

all_results_plt_df %>% reshape2::melt() %>% ggplot(aes(x=variable,y=value,fill=cluster))+geom_boxplot()+ scale_y_log10()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

+facet_wrap(~cluster,scales = 'free_y',ncol=1) 

avg_profiles <- as_tibble(do.call(rbind,lapply(1:6,function(c){

  km_c <- all_results_data[all_results_umap_df %>% filter(cluster==c) %$% label,,drop=F]
  avg_profile <- colMeans(km_c)
  return(avg_profile)
})))



km_results %>% filter(KM_HZ < 1)

c5 <- km_results_umap_df %>% filter(cluster==5) %$% label
km_results_c5 <- km_results_data[c5,]
c5_avg_profile <- colMeans(km_results_c5)



km_results_c1 <- km_results_data[km_results_umap_df %>% filter(cluster==1) %$% label,]
c1_avg_profile <- colMeans(km_results_c1)


p <- princomp(km_results_data)

p$loadings[,1:2] %>% as_tibble %>% ggplot(aes(x=Comp.1,y=Comp.2)) + geom_point()
