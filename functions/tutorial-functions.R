# A) Stouffer integration
doStouffer=function( mat ) {
    # mat should have features on rows and cells on columns
    res <- rowSums(mat) / sqrt( ncol(mat) )
    return(res)
}


# B) Stouffer integration on a cluster basis
StoufferClusters=function( mat, clusters ) 
{
 # mat is a matrix with features on rows and cells on columns
 # clusters is a named vector of clusters, where each element is a specific cell
  
  my_list <- list()
  for( a_clust_id in seq_len(nlevels(clusters)) )
  {
    index <- colnames(mat) %in% names(clusters[clusters == levels(clusters)[a_clust_id]])
    print(sum(index))
    my_list[[a_clust_id]] <- doStouffer( mat = mat[ , index ] )
  }
  names(my_list) <- seq_len(nlevels(clusters))
  mat_stouffered <- do.call(cbind,my_list)	
  colnames(mat_stouffered) <- levels(clusters)
  return( mat_stouffered )
}



# C) Plotting heatmaps
#Function to draw heatmap of genes in seurat object dat, with cluster labels clust
#seurat object must have labels "tissue" and "patient" 
PAHeatmap=function(genes,expression_to_plot,Clusters,colors){
  require(pheatmap)
  annotation_columns<-data.frame(Clusters)
  rownames(annotation_columns)<-colnames(expression_to_plot)
  colnames(annotation_columns)<-c("clusters")
  
  anno_colors <- list(clusters = colors)
  names(anno_colors[[1]]) <- levels(annotation_columns$cluster)
  
  cellType_Order<-unique(Clusters)
  idx<-c()
  for (i in 1:length(cellType_Order)){
    idx<-c(idx, which(annotation_columns[,1]==cellType_Order[i]))
  }
  
  expression_to_plot=expression_to_plot[,idx]
  annotation_columns=data.frame(annotation_columns[idx,])
  rownames(annotation_columns)<-colnames(expression_to_plot)
  colnames(annotation_columns)=c("clusters")
  
  quantile_breaks <- function(xs, n = 10) {
    breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
    breaks[!duplicated(breaks)]
  }
  mat_breaks <- quantile_breaks(as.matrix(apply(expression_to_plot,1,function(x){(x-mean(x))/sd(x)})), n = 30)
  
  p<-pheatmap(expression_to_plot, cluster_rows=FALSE, 
              show_rownames=TRUE, cluster_cols=FALSE, 
              annotation_col=annotation_columns, 
              color = colorRampPalette(colors = c('blue', 'white', 'red'))(length(mat_breaks)),
              fontsize_row = 8, show_colnames = FALSE, annotation_colors = anno_colors,
              scale="row")
  print(p)
  
}




