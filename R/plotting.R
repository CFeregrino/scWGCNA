
#' Plots the gene / module WGCNA tree, or trees
#' 
#' This function will help us plot the gene / module WGCNA tree. 
#' @param scWGCNA.data scWGCNA.data. An scWGCNA.data object, as calculated by run.scWGCNA().
#' @param tree numeric. Which of all the trees in the iterations should be plotted?
#' @return Plots the last / current gene / module WGCNA dendrogram. If a tree is specified, it returns that tree.
#' @importFrom graphics par
#' @export
#' @examples
#' # Plot the WGCNA tree
#' scWGCNA.plotdendro(scWGCNA.data = MmLimbE155.scWGCNA)
#' 
#' # Plot the first WGCNA tree of the iteration
#' scWGCNA.plotdendro(MmLimbE155.scWGCNA, 1)
#' 

scWGCNA.plotdendro = function(scWGCNA.data, tree=length(scWGCNA.data$moduletrees.history)){
  
  my.trees = scWGCNA.data$moduletrees.history
  
  if (is.numeric(tree)) {
    
      # Plot the dendrogram and colors underneath
      WGCNA::plotDendroAndColors(my.trees[[tree]][[1]], my.trees[[tree]][[2]], "Modules",
                                 dendroLabels = NULL,
                                 cex.dendroLabels = 0.6,
                                 addGuide = TRUE,
                                 main = "Gene dendrogram and module colors",
                                 guideAll = F)
      


    }
    
    
    
  }
  
#' Plots mean expression of the modules on the single cell data.
#' 
#' This function will help us plot the mean expression of the different modules calculated by scWGCNA
#' @param s.cells Seurat object. The single cell data used to run scWGCNA.
#' @param scWGCNA.data scWGCNA.data. An scWGCNA.data object, as calculated by run.scWGCNA().
#' @param modules numeric or "all". Which modules should be plotted? Could be a single module, by number or color, a vector of several modules, or "all" to plot them all in a grid
#' @param reduction string. The reduction to use for plots. It must be present in the @@reductions slot of sc.data. Default is "tsne"
#' @param ncol numeric. If plotting more than 1 module, how many columns should the grid present?
#' @return A ggplot, in case we're plotting only one module, or a grid of plots if plotting several
#' @export
#' @examples
#' # Plot the expression of all modules
#' scWGCNA.plotexpression(s.cells = my.small_MmLimbE155, scWGCNA.data = MmLimbE155.scWGCNA, modules = "all", ncol=3, reduction="tsne")
#' 

scWGCNA.plotexpression = function(s.cells, scWGCNA.data, modules = "all", reduction="tsne", ncol=2){
  
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
      ggplot2::theme_classic() +
      ggplot2::theme(legend.position="none",
                     axis.title = ggplot2::element_blank(),
                     axis.text = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank()) +
      ggplot2::scale_colour_gradientn(colours = c("gray90", "gray90", my.cols[i], my.cols[i])) +
      ggplot2::labs(colour=my.cols[i]) +
      ggplot2::ggtitle(paste0(i," ",my.cols[i]))
    
  }
  
  gridExtra::grid.arrange(grobs=my.plots, ncol=ncol)
  
} else{
    
  if (is.character(modules)) {
    p.module = which(my.cols == modules)
  } else{p.module = as.numeric(modules)}
  
  ggplot2::ggplot(toplot[order(avgexp[,p.module]),],
                  ggplot2::aes_string(x=colnames(toplot)[1], 
                                      y=colnames(toplot)[2])) +
    ggplot2::geom_point(ggplot2::aes_string(color=avgexp[order(avgexp[,p.module]),p.module]), size=1) +
    ggplot2::scale_size(range = c(1, 1)) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_gradientn(colours = c("gray90", "gray90", my.cols[p.module], my.cols[p.module])) +
    ggplot2::labs(colour=my.cols[p.module]) +
    ggplot2::ggtitle(paste0("Mean expression of module ",p.module," ",my.cols[p.module]))
  
  }
  
}

#' Plots a fruchtermanreingold network of the co-expression modules.
#' 
#' This function will help us plot a fruchtermanreingold network of each of the co-expression modules calculated by scWGCNA.
#' @param scWGCNA.data scWGCNA.data. An scWGCNA.data object, as calculated by run.scWGCNA().
#' @param module numeric . Which module should be plotted?
#' @param gnames Data frame. If you're using gene IDs and no symbols, you might wanna provide a list of gene names for plotting. Two columns: 1= ids present in expression matrix, 2= names to appear in plots. Rownames= same as 1st row
#' @param ... Parameters passed to GGally::ggnet2()
#' @return A plot showing the relationships between the genes in a module. The layot represents a fruchtermanreingold network. The size of the nodes represent the relative membership of the gene, to the module. The size of the edges represent the relative topology overlap between two genes.
#' @importFrom grDevices col2rgb
#' @export
#' @examples
#' 
#' # Calculate the networks in the scWCGNA object
#' MmLimbE155.scWGCNA = scWGCNA.networks(MmLimbE155.scWGCNA)
#' 
#' # Plot one of the modules as network
#' scWGCNA.plotnetwork(MmLimbE155.scWGCNA, module=1)
#' 
#' # If you're using IDs instead of names, and have a separate gene names data.frame
#' scWGCNA.plotnetwork(MmLimbE155.scWGCNA, module=1, gnames = my.small_MmLimbE155@@misc$gnames)

scWGCNA.plotnetwork = function(scWGCNA.data, module = 1, gnames=NULL, ...){
  
  if (is.null(scWGCNA.data[["networks"]])) {
    stop("Networks have not been generated, please run scWGCNA.networks() on your object first")
  }
  
  my.cols = levels(as.factor(scWGCNA.data[["dynamicCols"]]))
  
  if (is.character(module)) {
    module = which(my.cols == module)
    }
  
  lcols = my.cols[module]
  if (apply(col2rgb(my.cols[module]), 2,
            function(x) (x[1]*0.299 + x[2]*0.587 + x[3]*0.114)) < 75) {
    lcols = "white"
  } else {lcols = "black"}
  
  mynet = scWGCNA.data[["networks"]][[module]]
  
  if (is.null(gnames)) {
    gnames = network::network.vertex.names(mynet)
  } else (gnames = gnames[network::network.vertex.names(mynet),2])
  
  set.seed(42)
  
  GGally::ggnet2(mynet,
                 mode = "fruchtermanreingold",
                 layout.par = list(repulse.rad=network::network.size(mynet)^1.1,
                                         area=network::network.size(mynet)^2.3), # Give space in the middle
                 node.size = network::get.vertex.attribute(mynet,'membership01'), max_size = 20, #The size of the nodes
                 node.color = my.cols[module], # The color of the module
                 edge.size = "weight02",
                 edge.color = "black",
                 edge.alpha = network::get.edge.attribute(mynet,'weight01'),
                 ...
  ) +
    ggplot2::theme(legend.position="none") +
    ggplot2::geom_label(ggplot2::aes(label=gnames),
                        fill = my.cols[module],
                        alpha = 0.3,
                        color=lcols,
                        fontface = "bold")
  
}

