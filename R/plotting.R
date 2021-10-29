

plot.WGCNA.dendro = function(scWGCNA.data, history=F){
  
  my.trees = scWGCNA.data$moduletrees.history
  
  if (history = F) {
    my.trees = my.trees[[length(my.trees)]]
    
    WGCNA::plotDendroAndColors(my.trees[[1]], my.trees[[2]], "Modules",
                               dendroLabels = NULL,
                               cex.dendroLabels = 0.6,
                               addGuide = TRUE,
                               main = "Gene dendrogram and module colors",
                               guideAll = F)
    
    par(mfrow=c(1,1))
    
  }
  
  if (history = T) {
    
    my.plots = lapply(my.trees, function(t){
      
      # Plot the dendrogram and colors underneath
      WGCNA::plotDendroAndColors(t[[1]], t[[2]], "Modules",
                                 dendroLabels = NULL,
                                 cex.dendroLabels = 0.6,
                                 addGuide = TRUE,
                                 main = "Gene dendrogram and module colors",
                                 guideAll = F)

    })
    
    
    
  }
  
  
  
}

##################

plot.scWGCNA.multiexpression = function(s.cells, scWGCNA.data, modules = "all", reduction="tsne", ncol){
  
  my.cols = levels(as.factor(scWGCNA.data[["dynamicCols"]]))
  toplot = data.frame(Seurat::Embeddings(s.cells[[reduction]]))
  avgexp = scWGCNA.data[["sc.MEList"]]$averageExpr
  
if (modules == "all" | length(modules) > 1){
  
  if (modules == "all") {
    modules = 1:length(my.cols)
  } else{
    
    if (is.character(modules)) {
      modules = which(my.cols == modules)
    }
    
  }
  
  my.plots = list()
  
  for (i in modules) {
    
    my.plots[[my.cols[i]]] = ggplot2::ggplot(toplot[order(avgexp[,i]),],
                                             ggplot2::aes_string(x=colnames(toplot)[1], 
                                                                 y=colnames(toplot)[2])) +
      ggplot2::geom_point(ggplot2::aes_string(color=avgexp[order(avgexp[,i]),i]), size=1) +
      ggplot2::scale_size(range = c(1, 1)) +
      ggplot2::theme_void() +
      ggplot2::theme(legend.position="none") +
      ggplot2::scale_colour_gradientn(colours = c("gray90", "gray90", my.cols[i], my.cols[i])) +
      ggplot2::labs(colour=my.cols[i])
    
  }
  
  gridExtra::grid.arrange(grobs=xx, ncol=4)
  
} else{
    
  if (is.character(modules)) {
    p.module = which(my.cols == modules)
  }
  
  ggplot2::ggplot(toplot[order(avgexp[,p.module]),],
                  ggplot2::aes_string(x=colnames(toplot)[1], 
                                      y=colnames(toplot)[2])) +
    ggplot2::geom_point(ggplot2::aes_string(color=avgexp[order(avgexp[,p.module]),p.module]), size=1) +
    ggplot2::scale_size(range = c(1, 1)) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position="none") +
    ggplot2::scale_colour_gradientn(colours = c("gray90", "gray90", my.cols[p.module], my.cols[p.module])) +
    ggplot2::labs(colour=my.cols[p.module])
  
  }
  
}

############################

plot.scWGCNA.network = function(scWGCNA.data, gnames){
  
  if (is.null(scWGCNA.data[["networks"]])) {
    stop("Networks have not been generated, please run scWGCNA.networks() on your object first")
  }
  
  
  
}





  
