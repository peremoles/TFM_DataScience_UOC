library(abind)
library(parallel)
library(reshape2)
library(dplyr)
library(tidyr)

source("~/Documents/spatial-rna-AJ/spGSVA_funcs.R")


compute_gradients <- function(spe,res){
  df<-as.data.frame(cbind(t(assay(res)),int_colData(spe)$spatialData[,2:3]))
  df <- df %>% arrange(array_row, array_col)
  genesets<-colnames(df)[-c(length(colnames(df))-1,length(colnames(df)))]
  
  gsva_list <-lapply(genesets, function(x){
    prov<-df[,c(which(genesets==x),145,146)]
    matrix_df <- prov %>% pivot_wider(names_from = array_col, values_from = x)
    matrix_df[is.na(matrix_df)] <- 0
    matrix_df<-as.matrix(matrix_df)
    rownames(matrix_df) <- matrix_df[, 1]
    matrix_df <- matrix_df[, -1]
    sorted_row_names <- sort(as.numeric(rownames(matrix_df)))
    sorted_col_names <- sort(as.numeric(colnames(matrix_df)))
    matrix_df <- matrix_df[as.character(sorted_row_names),as.character(sorted_col_names)]
    matrix_df
  })
  
  gkh <- gaussianKernelHex(5, sigma.x=.5, sigma.y=.5)
  
  # for each gene perform convolution, and calculate gradient using Sobel kernels
  
  system.time(graddata_gsva <- mclapply(gsva_list, function(M){ computeGradientHexGen(convolution2D(M, gkh), radius=2) }, mc.cores = 6))
  
  
  # separate angle and magnidute of gradients
  angledata_gsva <- lapply(graddata_gsva, function(x) { 
    x <- x$angle
    x[is.na(x)] <- 0
    x
  })
  
  magdata_gsva <- lapply(graddata_gsva, function(x) { 
    x <- x$magnitude
    x[is.na(x)] <- 0
    colnames(x) <- colnames(gsva_list[[1]])
    rownames(x) <- rownames(gsva_list[[1]])
    x
  })
  
  coords <- int_colData(spe)$spatialData[,2:3]
  coords$spot<-rownames(coords)
  
  long_magdata<-lapply(magdata_gsva, function(x){
    m<-melt(x)
    f<-merge(m,coords, by.x = c("Var1", "Var2"), by.y = c("array_row", "array_col"))
    rownames(f) <- f$spot
    f$spot <- NULL
    f$Var1 <- NULL
    f$Var2 <- NULL
    colnames(f)<-"GSVA_gradient_magnitude"
    f
  })
  
  spots <- rownames(coords)
  magdata_matrix <- do.call(cbind, lapply(long_magdata, function(mat) mat[,1]))
  rownames(magdata_matrix) <- spots
  colnames(magdata_matrix) <-genesets
  return(magdata_matrix)
}