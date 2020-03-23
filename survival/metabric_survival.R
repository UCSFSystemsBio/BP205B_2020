source('/wynton/home/students/snanda/rds/bp205/analysis/common_functions.R')
library(survival)
library(survminer)

read_metabric_data <- function(fix_names=F){
  data <- fread('/wynton/scratch/bp205/METABRIC/METABRIC_MedCen_Collapsed_GEOannot(n=1992).merge.txt')
  
  if(fix_names==TRUE){
    checked <- suppressWarnings(HGNChelper::checkGeneSymbols(data$V1)) %>% 
      mutate(to_keep = !(str_detect(Suggested.Symbol,'//') | is.na(Suggested.Symbol)))
    
    data2  <- data[checked$to_keep,]
    data2$V1 <- checked$Suggested.Symbol[checked$to_keep]
    
    duplicated_names <- names(which(table(data2$V1) > 1))
    duplicated_to_drop <- unlist(lapply(duplicated_names,function(name){
      which(data2$V1 %in% name)[-1]
    }))
    data2 <- data2[-duplicated_to_drop,,drop=FALSE]
  }else{
    data2 <- data
  }
  
  data3 <- dcast(melt(data2,id.vars = 'V1'),formula = as.character(variable)~V1)
  setDF(data3)
  
  colnames(data3)[1] <- 'patient'
  
  return(data3)
}

read_metabric_metadata <- function(censor = 120){
  metadata <- readxl::read_excel('/wynton/scratch/bp205/METABRIC/brca_metabric_clinical_data.xlsx') %>%
    dplyr::rename(patient=`Patient ID`,sampleID=`Sample ID`,
           vital_status = `Patient's Vital Status`,
           overall_survival_months=`Overall Survival (Months)`) %>%
    mutate(patient= str_replace(patient,'-',''),
           sampleID = str_replace(sampleID,'-',''),
           vital_status = ifelse(vital_status=='Died of Disease',1,0),
           vital_status = ifelse(overall_survival_months>censor,0,vital_status),
           overall_survival_months = ifelse(overall_survival_months>censor,censor,overall_survival_months),
           overall_survival_years = overall_survival_months/12,
           Subtype = factor(`Pam50 + Claudin-low subtype`,levels =  c('Normal','Basal','Her2','LumA','LumB','claudin-low','NC')),
           Cellularity = factor(Cellularity,levels=c('low','moderate','high'))
           )
  return(metadata)
}

select_and_score_metagene <- function(df,metagene,facets=c()){
  default_facets<- c('patient','overall_survival_months','overall_survival_years','vital_status')
  added_facets <- facets
  facets <- c(default_facets,added_facets)

  if(is.character(metagene)){
    genes <- metagene[metagene %in% colnames(df)]
    sign <- rep(1,length(genes))
    names(sign) <- genes
  }else{
    genes <- names(metagene)[names(metagene) %in% colnames(df)]
    sign <- metagene[names(metagene) %in% colnames(df)]
  }
  scored <- df %>%
    select_at(vars(c(genes,facets))) %>%
    pivot_longer(all_of(genes)) %>%
    group_by_at(vars(facets)) %>%
    summarize(activation=sum(sign[name]*value)) %>%
    #summarize(activation=sum(value*sign[name])) %>%
    ungroup %>% mutate_at(vars(added_facets),as.factor)
  return(scored)
}

stratify_activation <- function(df,ntiles=2,threshold=T){
  
  if(ntiles == 2 && threshold==TRUE){
    ## Do the iteration over the tresholds to optimize Up/Down partiioning
    thresholds <- sort(df$activation)
    thresholds <- thresholds[10:(length(thresholds)-10)]
    
    thresh_surv <- Surv(df$overall_survival_years,df$vital_status)
    message('Running threshold optimization')
    pvals <- unlist(lapply(thresholds,function(t){
      updn <- as.numeric(df$activation>t)
      diff <- survdiff(thresh_surv~updn)
      pv <- (pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE))
      return(pv)
    }))
    final_threshold <-thresholds[which.max(-log10(pvals))]
    strata <- factor(ifelse(df$activation>final_threshold,'High','Low'),levels=c('Low','High'))
    
  }else if(ntiles==3){
    strata <- factor(c('Hi','Med','Lo')[dplyr::ntile(df$activation,n_tiles)],levels=c('Lo','Med','Hi'))
  }else{
    ## For ntiles==2, this is median partitioning
    strata <- as.factor(dplyr::ntile(df$activation,ntiles))
  }
  df$strata <- strata
  return(df)
}

#########################################################################
### Implement wrapper function for running all signatures and scoring them 
score_stratify_fit_metagene <- function(df,metagene_name,metagene,ntiles=2,threshold=T,other_facets,outdir){
  print(metagene_name)
  # Create the directories if not present
  lapply(c('km_plots','result_out','hazard_plots'),function(path){
    dir.create(paste0(outdir,'/',path,'/'),showWarnings = F)
    return(NULL)
  })

  km_plot_path <- paste0(outdir,'/km_plots/',metagene_name,'_km_plot.pdf')
  result_out_path <- paste0(outdir,'/result_out/',metagene_name,'_result_out.csv')
  hazard_plot_path <- paste0(outdir,'/hazard_plots/',metagene_name,'_hazard_plot.pdf')
  
  ## Generate the dataset by scoring the metagene
  test <- df %>%
    select_and_score_metagene(metagene,facets=other_facets) %>%
    stratify_activation(ntiles = 2,threshold = T)
  
  ## Run Kaplan-Meir analysis
  message('Running survival analysis')
  km_fit <- survminer::surv_fit(survival::Surv(overall_survival_years, vital_status) ~ strata, data=test)
  km_summary <- summary(coxph(Surv(overall_survival_years, vital_status) ~ strata,data=test))
  #km_summary$coefficients[2]
  km_pval <- c(survminer::surv_pvalue(km_fit)$pval,
               km_summary$conf.int[c(3,1,4)],
               km_fit$n)
  names(km_pval) <- c('KM_pval','KM_HZ_0.95_lower','KM_HZ','KM_HZ_0.95_upper','n_low','n_high')
  
  km_plt <- arrange_ggsurvplots(list(ggsurvplot(km_fit,
                                                pval=T,conf.int = T,risk.table = T,color = 'strata',palette = c('#3BA805','#E8634A'),data = test)),
                                print=FALSE,ncol=1,nrow=1
  )
  ggsave(plot = km_plt,filename = km_plot_path,width = 8,height = 6)
  
  ## Run COX anlaysis using continous activation as covariate
  message('Running cox fit')
  test$activation_Z <- c(scale(test$activation))
  cox_fit <- coxph(Surv(overall_survival_years, vital_status) ~ activation_Z +`Chemotherapy` + Subtype + Cellularity, data=test,model = T)
  cox_pval <- c(summary(cox_fit)$coef[,c(2,5)],recursive=T)
  names(cox_pval) <- c(paste0('Cox_HZ_',rownames(summary(cox_fit)$coef)),paste0('Cox_pval_',rownames(summary(cox_fit)$coef)))
  cox_plt <- suppressWarnings(ggforest(cox_fit,data=test))
  
  ggsave(plot = cox_plt,filename = hazard_plot_path, width=10,height=7)
  
  ## aggegate result
  result <- c(km_pval,cox_pval)
  write.table(result,file=result_out_path,sep = ',')
  
  rm(test)
  return(result)
}

##################
## Run the analysis
##################
data <- read_metabric_data(fix_names = T)
metadata <- read_metabric_metadata()
df <- inner_join(data,metadata,"patient")


hcc_mk_DE_results <- list.files('/wynton/scratch/bp205/signatures/BP205B_2020/DE/HCC_MK/FINAL/',
                                recursive = T,full.names = T,pattern = '.csv')
mda_DE_results <- list.files('/wynton/scratch/bp205/signatures/BP205B_2020/DE_results/Final/',
                             recursive = T,full.names = T,pattern = '.csv')

hcc_sig_set <- load_signature_set('HCC',hcc_mk_DE_results,colnames(df),abs = F)
mda_sig_set <- load_signature_set('MDA',mda_DE_results,colnames(df),abs = T)
all_sig_set <- c(mda_sig_set,hcc_sig_set)
all_sig_set <- all_sig_set[sapply(all_sig_set,length) != 0]

f <- list.files('/wynton/home/students/snanda/rds/bp205/analysis/survival/results/result_out//') %>% str_replace('_result_out\\.csv','')
# 
all_sig_set <- all_sig_set[!(names(all_sig_set) %in% f)]


other_facets <- c('Cancer Type' , 'Cellularity' , 'Chemotherapy' ,'ER Status' , 'HER2 Status' , 'PR Status','Tumor Stage' , 'Age at Diagnosis' , 'Subtype')

outdir <- '/wynton/home/students/snanda/rds/bp205/analysis/survival/results/'
rm(data,hcc_sig_set,mda_sig_set,metadata)

# ## testing:###############################################################
# metagene <- all_sig_set$`HCC_6_vs_p-H_UP_FC_1`
# ## Generate the dataset by scoring the metagene
# test <- df %>%
#   select_and_score_metagene(metagene,facets=other_facets) %>%
#   stratify_activation(ntiles = 2,threshold = T)
# 
# ## Run Kaplan-Meir analysis
# message('Running survival analysis')
# km_fit <- survminer::surv_fit(survival::Surv(overall_survival_years, vital_status) ~ strata, data=test)
# km_summary <- summary(coxph(Surv(overall_survival_years, vital_status) ~ strata,data=test))
# km_pval <- c(survminer::surv_pvalue(km_fit)$pval,
#              km_summary$conf.int[c(3,1,4)],
#              km_fit$n)
# names(km_pval) <- c('KM_pval','KM_HZ_0.95_lower','KM_HZ','KM_HZ_0.95_upper','n_low','n_high')
# 
# km_plt <- arrange_ggsurvplots(list(ggsurvplot(km_fit,                                                           ## High , low
#                                               pval=T,conf.int = T,risk.table = T,surv.co,palette = c('#3BA805','#E8634A'),data = test)),
#                               print=FALSE,ncol=1,nrow=1
# )
# 
# 
# ## Run COX anlaysis using continous activation as covariate
# message('Running cox fit')
# test$activation_Z <- c(scale(test$activation))
# cox_fit <- coxph(Surv(overall_survival_years, vital_status) ~ activation_Z +`Chemotherapy` + Subtype + Cellularity, data=test,model = T)
# cox_pval <- c(summary(cox_fit)$coef[,c(2,5)],recursive=T)
# names(cox_pval) <- c(paste0('Cox_HZ_',rownames(summary(cox_fit)$coef)),paste0('Cox_pval_',rownames(summary(cox_fit)$coef)))
# cox_plt <- suppressWarnings(ggforest(cox_fit,data=test))

# 
# 
# 
###############################

models <-mapply(FUN=score_stratify_fit_metagene,
                            metagene_name=names(all_sig_set),
                            metagene = all_sig_set,
                            ntiles=2,threshold=T,outdir = outdir,SIMPLIFY = F,
                            MoreArgs = list(df=df,other_facets = other_facets))

models_files <- list.files('/wynton/home/students/snanda/rds/bp205/analysis/survival/results/result_out//',full.names = T)

models_read <- do.call(cbind,lapply(models_files,function(d){
  read.table(d,header = T,sep = ',')
}))

colnames(models_read) <- str_remove(basename(models_files),'_result_out\\.csv')
models_df <- t(models_read) %>% as.data.frame %>% rownames_to_column() %>% dplyr::rename(signature=rowname)

data.table::fwrite(models_df,paste0(outdir,'survival_results.csv'),sep = ',')

