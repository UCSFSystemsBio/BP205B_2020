## Load shared libraries
library(data.table)
library(tidyverse)
library(ggplot2)
library(MASS)
library(ggplot2)
library(gridExtra)
library(readxl)
library(HGNChelper)

create_signature_set <- function(dataset,DE_results,ref_gene_set,abs=T){
  filenames <- stringr::str_sub(basename(DE_results),1,-5)
  
  allSigSets <- unlist(lapply(1:length(DE_results),function(f){
    gmt <- data.table::fread(DE_results[f],sep=',')
    sigData <- sign(gmt$logfoldchange)
    names(sigData) <- gmt$genes
    matched <- names(sigData) %in% ref_gene_set
    
    print(paste0('Unmatched ',f,' ',filenames[f],' ',sum(!matched)/nrow(gmt)))
    
    sigSets <- unlist(lapply(c(0,1,2),function(thresh){
      fold_data <- sigData[abs(gmt$logfoldchange)>=thresh & matched]
      
      ## If the metagene is all -1s, then flip to positive
      if(abs && all(fold_data==-1)){
        fold_data <- abs(fold_data)
      }
      sig_set_name <- paste0(dataset,'_',filenames[f],'_FC_',thresh)
      
      a <- VISION::createGeneSignature(name = sig_set_name,sigData = fold_data)
      
      retlist <- list(a)
      names(retlist) <- sig_set_name
      
      return(retlist)
    }))
    return(sigSets)
  }))
  return(allSigSets)
}

load_signature_set <- function(dataset,DE_results,ref_gene_set,abs=T){
  filenames <- stringr::str_sub(basename(DE_results),1,-5)
  
  allSigSets <- unlist(lapply(1:length(DE_results),function(f){
    gmt <- data.table::fread(DE_results[f],sep=',')
    sigData <- sign(gmt$logfoldchange)
    names(sigData) <- gmt$genes
    matched <- names(sigData) %in% ref_gene_set
    
    print(paste0('Unmatched ',f,' ',filenames[f],' ',sum(!matched)/nrow(gmt)))
    
    sigSets <- unlist(lapply(c(0,1,2),function(thresh){
      fold_data <- sigData[abs(gmt$logfoldchange)>=thresh & matched]
      
      ## If the metagene is all -1s, then flip to positive
      if(abs && all(fold_data==-1)){
        fold_data <- abs(fold_data)
      }
      
      sig_set_name <- paste0(dataset,'_',filenames[f],'_FC_',thresh)
      
      a <- fold_data
      retlist <- list(a)
      names(retlist) <- sig_set_name
      return(retlist)
      
    }),recursive = F)
    return(sigSets)
  }),recursive = F)
  return(allSigSets)
}


h5read <- function(h5path){
  hfile <- hdf5r::h5file(filename = h5path, mode = 'r')
  
  m <- as.matrix(Matrix::sparseMatrix(i=hfile[['X']][['indices']][]+1,
                                      p=hfile[['X']][['indptr']][],
                                      x=hfile[['X']][['data']][]
  ))
  hdf5r::h5close(hfile)
  return(m)
}

h5read_umap <- function(h5path){
  hfile <- hdf5r::h5file(h5path,mode = 'r')
  
  umap <- as.data.frame(t(hfile[['obsm']][['X_umap']][,]))
  
  hdf5r::h5close(hfile)
  return(umap)
  
}

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

generate_principal_curve <- function(cds,reduction_method='UMAP'){
  reduction_method <- 'UMAP'
  x <- 1
  y <- 2
  ica_space_df <- t(cds@principal_graph_aux[[reduction_method]]$dp_mst) %>% 
    as.data.frame() %>% 
    dplyr::select_(prin_graph_dim_1 = x, 
                   prin_graph_dim_2 = y) %>% 
    dplyr::mutate(sample_name = rownames(.), sample_state = rownames(.))
  
  
  
  dp_mst <- cds@principal_graph[[reduction_method]]
  edge_df <- dp_mst %>% 
    igraph::as_data_frame() %>% dplyr::select_(source = "from", target = "to") %>%
    dplyr::left_join(ica_space_df %>% 
                       dplyr::select_(source = "sample_name", source_prin_graph_dim_1 = "prin_graph_dim_1",
                                      source_prin_graph_dim_2 = "prin_graph_dim_2"), 
                     by = "source") %>% dplyr::left_join(ica_space_df %>%  dplyr::select_(target = "sample_name", target_prin_graph_dim_1 = "prin_graph_dim_1", target_prin_graph_dim_2 = "prin_graph_dim_2"),  by = "target")
  
  
  g <- ica_space_df %>% ggplot(aes(x=prin_graph_dim_1,y=prin_graph_dim_2)) + geom_point()
  g <- g + geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                                   y = "source_prin_graph_dim_2",
                                   xend = "target_prin_graph_dim_1", 
                                   yend = "target_prin_graph_dim_2"), 
                        size = 0.1, 
                        color = 'blue',
                        linetype = "solid", 
                        na.rm = TRUE, 
                        data = edge_df)
  
  return(g)  
  
}


get_node_df <- function(cds,reduction_method='UMAP'){
  
    ica_space_df <- t(cds@principal_graph_aux[[reduction_method]]$dp_mst) %>% 
      as.data.frame() %>% 
      dplyr::select_(prin_graph_dim_1 = x, 
                     prin_graph_dim_2 = y) %>% 
      dplyr::mutate(sample_name = rownames(.), sample_state = rownames(.))
    
    
    return(ica_space_df)
}

