library(RColorBrewer)
library(ggplot2)

# This function gets a SpatialExperiment object with the clustering labels and a SpatialExperiment object with the gsva results 
# and returns a dataframe where the rows are the spots and the columns are the gsva scores for each geneset and its cluster label

get_cluster_gsva_df <- function(spe, res) {
  clusters_df <- colData(spe)[,10,drop=FALSE]
  gsva_df<-t(as.data.frame(assays(res)$es))
  df <- cbind(gsva_df,clusters_df)
  return(df)
}
# This function gets a SpatialExperiment object with the clustering labels and a SpatialExperiment object with the gsva results 
# and returns the clusters more spatially correlated to gsva scores comparing the inside and outside of the cluster using t.test, and their t-statistics.
gsva_clusters_ttests <- function(spe,res,geneset){
  df <- get_cluster_gsva_df(spe,res)
  tests <- sapply(1:length(unique(df$label)), function(x){
    cluster_to_test <- x
    df$label_binary <- ifelse(df$label==cluster_to_test,1,0)
    df_0 <- subset(df, label_binary == 0)
    df_1 <- subset(df, label_binary == 1)
    
    test <-t.test(df_0[[geneset]], df_1[[geneset]], alternative = "less", var.equal = FALSE)
    test
  })
  return(list(clusters = which(p.adjust(unlist(tests[3,]),method = "bonferroni")<0.05), t_statistics = unlist(tests[1,which(tests[3,]<0.05)])))
}

# This function gets a SpatialExperiment object with the clustering labels and a SpatialExperiment object with the gsva results 
# and returns the geneset more correlated to a each cluster.

most_cor_geneset <- function(spe,res){
  df <- get_cluster_gsva_df(spe,res)
  most_corr <- sapply(1:length(unique(df$label)), function(x){
    cluster_to_test <- x
    df$label_binary <- ifelse(df$label==cluster_to_test,1,0)
    correlation_results <- sapply(df[1:(length(df)-2)], function(y) cor.test(y, df$label_binary,"greater"))
    most_corr <- names(which.max(correlation_results[1,]))
    most_corr
    })
  }



gsva_clusters_ttests <- function(spe,res,geneset){
  df <- get_cluster_gsva_df(spe,res)
  tests <- sapply(1:length(unique(df$label)), function(x){
    cluster_to_test <- x
    df$label_binary <- ifelse(df$label==cluster_to_test,1,0)
    
    # Plot comparison
    ggplot(data = df, aes(x=as.factor(label_binary), y= geneset)) + 
      geom_boxplot()
    df_0 <- subset(df, label_binary == 0)
    df_1 <- subset(df, label_binary == 1)
    
    test <-t.test(df_0[[geneset]], df_1[[geneset]], alternative = "less", var.equal = FALSE)
    test
  })
  return(list(clusters = which(p.adjust(unlist(tests[3,]),method = "bonferroni")<0.05), t_statistics = unlist(tests[1,which(tests[3,]<0.05)])))
}

boxplot_cluster_gsva <-function(spe,res, geneset, cluster){
  df <- get_cluster_gsva_df(spe,res)
  df$label_binary <- ifelse(df$label==cluster,1,0)
  ggplot(data = df, aes(x=as.factor(label_binary), y= .data[[geneset]])) + geom_boxplot() + scale_x_discrete(labels = c("Outside cluster", "Inside cluster")) + 
    labs(title=paste0("Boxplot comparison for geneset ",geneset," inside and outside cluster ",cluster),x = "",y = "GSVA score")
}

