
#' Plots the gene / module WGCNA tree, or trees
#' 
#' This function will help us plot the gene / module WGCNA tree. 
#' @param scWGCNA.data scWGCNA.data. An scWGCNA.data object, as calculated by run.scWGCNA().
#' @param tree numeric. Which of all the trees in the iterations should be plotted?
#' @return Plots the last / current gene / module WGCNA dendrogram. If a tree is specified, it returns that tree.
#' @importFrom graphics par
#' @export
#' @examples
#' 
#' # A pre-calculated WGCNA data list
#' class(MmLimbE155.scWGCNA)
#' 
#' # Plot the last, and current WGCNA tree
#' scW.p.dendro(scWGCNA.data = MmLimbE155.scWGCNA)
#' 
#' # Plot the first WGCNA tree of the iteration
#' scW.p.dendro(MmLimbE155.scWGCNA, 1)
#' 

# The function, taking the history within the scW data
scW.p.dendro = function(scWGCNA.data, tree=length(scWGCNA.data$moduletrees.history)){
  
  # Take all the trees that were recorded
  my.trees = scWGCNA.data$moduletrees.history
  
  
  # Plot the dendrogram and colors underneath
      WGCNA::plotDendroAndColors(my.trees[[tree]][[1]], my.trees[[tree]][[2]], "Modules",
                                 dendroLabels = NULL,
                                 cex.dendroLabels = 0.6,
                                 addGuide = TRUE,
                                 main = "Gene dendrogram and module colors",
                                 guideAll = F)
      



    
    
    
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
#' # A pre-analyzed Seurat object, subsampled
#' my.small_MmLimbE155
#' MmLimb.sc = my.small_MmLimbE155
#' 
#' # A pre-calculated WGCNA data list
#' class(MmLimbE155.scWGCNA)
#' MmLimb.scW = MmLimbE155.scWGCNA
#' 
#' # Plot the expression of all modules
#' scW.p.expression(s.cells=MmLimb.sc, scWGCNA.data=MmLimb.scW, modules="all", reduction="tsne")
#' 

# The funciton, taking single cell object, a scW data, which modules to plot, and which reduction to use.
scW.p.expression = function(s.cells, scWGCNA.data, modules = "all", reduction="tsne", ncol=2){
  
  # Take the colors from the scW data
  my.cols = levels(as.factor(scWGCNA.data[["dynamicCols"]]))
  # Take the desired reduction
  toplot = data.frame(Seurat::Embeddings(s.cells[[reduction]]))
  # The sc expression from the scW data
  avgexp = scWGCNA.data[["sc.MEList"]]$averageExpr

  # If we want to print all the modules, OR more than 1  
if (modules == "all" | length(modules) > 1){
  
  # If it says "all", let's change that for actual numbers
  if (modules == "all") {
    modules = 1:length(my.cols)
  } else{
    # If names of modules were provided, we also change them to numbers
    if (is.character(modules)) {
      modules = which(my.cols == modules)
    }
    
  }
  
  # A list to keep all the plots
  my.plots = list()
  
  # We go trough the modules we have chosen
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
  
  # Plot them as a grid
  gridExtra::grid.arrange(grobs=my.plots, ncol=ncol)
  
} else{
  
  # Here, if we only want one module plotted
  # Again, change whatever might be, to numbers
  if (is.character(modules)) {
    p.module = which(my.cols == modules)
  } else{p.module = as.numeric(modules)}
  
  # Plot it!
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
#' scW.p.network(MmLimbE155.scWGCNA, module=1)
#' 
#' # If you're using IDs instead of names, and have a separate gene names data.frame
#' scW.p.network(MmLimbE155.scWGCNA, module=1, gnames = my.small_MmLimbE155@@misc$gnames)

# The function, taking which module, and probable gnames
scW.p.network = function(scWGCNA.data, module = 1, gnames=NULL, ...){
  
  # Check if networks have been generated
  if (is.null(scWGCNA.data[["networks"]])) {
    stop("Networks have not been generated, please run scWGCNA.networks() on your object first")
  }
  
  # Take the colors from the modules
  my.cols = levels(as.factor(scWGCNA.data[["dynamicCols"]]))
  
  # If they are named by name, make it into numbers
  if (is.character(module)) {
    module = which(my.cols == module)
    }
  
  # Check if the labels of plots should be black or white, depending on how dark module colors are
  lcols = my.cols[module]
  # Is it too light?
  if (apply(col2rgb(my.cols[module]), 2,
            function(x) (x[1]*0.299 + x[2]*0.587 + x[3]*0.114)) < 75) {
    lcols = "white"
  } else {lcols = "black"}
  
  # From the pre-calculated networks, take the one you're gonna plot
  mynet = scWGCNA.data[["networks"]][[module]]
  
  # In case we need to change the names to be plotted
  if (is.null(gnames)) {
    gnames = network::network.vertex.names(mynet)
  } else (gnames = gnames[network::network.vertex.names(mynet),2])
  
  # We want to always have the same plot
  set.seed(42)
  
  # Do the plotting
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
    # The labels
    ggplot2::geom_label(ggplot2::aes(label=gnames),
                        fill = my.cols[module],
                        alpha = 0.3,
                        color=lcols,
                        fontface = "bold")
  
}

#' Plots a barplot with the fraction of the expressed and orthologous genes for each module of a scWGCNA comparison.
#' 
#' This function will plot a barplot showing what's the fraction of genes in each module that is (if its the case) present in the provided 1-2-1 orthologues table, and the genes that are expressed in all test samples, and therefore used for the comparison.
#' @param scWGCNA.comp.data scWGCNA comparative data as calculated by scWGNA.compare().
#' @return A ggplot barplot
#' @export
#' @examples
#' 
#' # A pre-calculated list scWGCNA comparative data, calculated with scWGCNA.compare
#' class(MmvGg.comparative)
#' 
#' # Plot the fraction of genes used for the comparison.
#' scW.p.modulefrac(MmvGg.comparative)

# The function, taking only the pre-calculated data
scW.p.modulefrac = function(scWGCNA.comp.data){
  
  # Take the precalculated number of genes lost
  to.plot = scWGCNA.comp.data$misc$modulefrac
  
  # Make them into decimal fractions
  to.plot[,-1] = to.plot[,-1] / to.plot[,2]
  
  # If we had genes lost to orthology
  if (ncol(to.plot) == 4) {
    # Make the difference, to plot them stacked
    to.plot$diff = to.plot$ortho - to.plot$expr
    # Get rid of that column, no needed
    to.plot = to.plot[,-3]
  }
  
  # melt it, to have it ggplot friendly
  to.plot = reshape::melt(to.plot[,c(1,3:ncol(to.plot))])
  
  # Plot it!
  ggplot2::ggplot(to.plot, 
                  ggplot2::aes(x=module, 
                               y=value, 
                               fill = variable,
                               color=module)) +
    ggplot2::geom_bar(position=ggplot2::position_stack(reverse = TRUE), 
                      stat = "identity",
                      size=2) + 
    ggplot2::scale_color_manual(values = as.character(to.plot$module)) + 
    ggplot2::scale_fill_manual(values=c("gray","gray93"),
                               labels = c("expressed in\nall samples","orthologous\nnon-expressed")) + 
    ggplot2::theme_light() + 
    ggplot2::guides(color="none") +
    ggplot2::labs(y="fraction", 
                  fill="Genes in module") + 
    ggplot2::theme(legend.text = ggplot2::element_text(margin = ggplot2::margin(t = 10, b=10)))
  
}

#' Plots dotplot showing the z-score values for preservation and other aspects of it
#' 
#' This function will plot a dotplot, with the zscore values for global preservation, as well as density and connectivity for each module and each test sample that was compared. Can also plot the median rank.
#' @param scWGCNA.comp.data scWGCNA comparative data as calculated by scWGNA.compare().
#' @param to.plot character or character vector. Which aspects of the preservation should be plotted? Options: "preservation", "median.rank", "density", "connectivity".
#' @param test.samples character or character vector. Which test samples to plot. Default is all
#' @return Either a single ggplot dotoplot showing the desired aspect of preservation. If several were requested, a gridExtra of the different ggplot dotplots.
#' @export
#' @examples
#' 
#' # S pre-calculated list scWGCNA comparative data, calculated with scWGCNA.compare
#' class(MmvGg.comparative)
#' 
#' # Plot the overall preservation and median rank.
#' scW.p.preservation(scWGCNA.comp.data = MmvGg.comparative, to.plot=c("preservation", "median.rank"))

# The function takinf the comparative data and what will be plotted.
scW.p.preservation = function(scWGCNA.comp.data, to.plot=c("preservation", "median.rank"), test.samples = NULL){
  
  # leave grey and gold modules out
  modColors = rownames(scWGCNA.comp.data$preservation$observed[[1]][[2]])
  plotMods = !(modColors %in% c("grey", "gold"))
  
  # The variable where we keep the data to plot
  plotData = data.frame(Rank = numeric(), Zsum= numeric(), Density = numeric(), Connectivity = numeric(),
                        Size = numeric(), Cols = character(), Sample = character())
  
  # Fill in with data from the samples
  for (test in 1:(length(scWGCNA.comp.data$quality$observed[[1]])-1)) {
    plotData = rbind(plotData,data.frame(Rank = scWGCNA.comp.data$preservation$observed[[1]][[test+1]][, 2],
                                         Zsum= scWGCNA.comp.data$preservation$Z[[1]][[test+1]][, 2],
                                         Density = scWGCNA.comp.data$preservation$Z[[1]][[test+1]][, 3],
                                         Connectivity = scWGCNA.comp.data$preservation$Z[[1]][[test+1]][, 4],
                                         Size = scWGCNA.comp.data$preservation$observed[[1]][[test+1]][, 1],
                                         Cols = rownames(scWGCNA.comp.data$preservation$observed[[1]][[test+1]]),
                                         Sample = rep(scWGCNA.comp.data$misc$testnames[test],nrow(scWGCNA.comp.data$preservation$observed[[1]][[test+1]]))))
  }
  
  # We kick out the golden and grey modules
  plotData = plotData[rep(plotMods, length(scWGCNA.comp.data$quality$observed[[1]])-1),]
  
  # For the density and connectivity, we will set the same limits
  my.lim = (max(plotData[,3:4])-min(plotData[,3:4]))*0.025
  my.lim = c(min(plotData[,3:4])-my.lim, max(plotData[,3:4])+my.lim)
  
  #For the median rank, we will show transparency for modules without evidence of preservation
  my.alpha = rep(1,dim(plotData)[1])
  my.alpha[which(plotData$Zsum < 2)] = 0.75
  plotData$alpha = my.alpha
  
  # The plotting of the Zsummary and Median Rank, against the module size
  my.p=list()
  
  # The labels for the plots
  my.plotData = plotData
  my.plotData$Sample = factor(my.plotData$Sample, levels = scWGCNA.comp.data$misc$testnames)
  my.plotData$mylabels = my.plotData$Cols
  my.plotData$mylabels[which(plotData$Zsum < 2)] = NA

  plotData = my.plotData #?????????
  
  # Do we plot only certain test samples?
  if (!is.null(test.samples)) {
    plotData = plotData[which(plotData$Sample %in% test.samples),]
  }
  
  # A global theme, for all plots
  my.theme = list(ggplot2::geom_text(nudge_x = 0.08, hjust=0, size=2.5, check_overlap = T),
    ggplot2::scale_fill_manual(values= as.character(unique(plotData$Cols))),
    ggplot2::scale_x_continuous(trans='log2', expand = c(0.1, 0)),
    ggplot2::theme_classic())
  
  # If we have too many samples
  if (length(unique(plotData$Sample)) > 5) {
    my.theme = c(my.theme, ggplot2::scale_shape_manual(name= "Test\nsample:",values= c(0:25),
                                                       labels=scWGCNA.comp.data$misc$testnames))
  } else{
    my.theme = c(my.theme, ggplot2::scale_shape_manual(name= "Test\nsample:",values= c(21:25),
                                                       labels=scWGCNA.comp.data$misc$testnames))
  }
  
  # Some other shared theme items
  my.gp = list(ggplot2::geom_point(ggplot2::aes(fill=Cols, shape=Sample), size=3.5),
               ggplot2::geom_hline(yintercept = c(2,10), linetype="dashed", color=c("red", "limegreen")))
  my.gpa = ggplot2::geom_point(ggplot2::aes(fill=Cols, shape=Sample, alpha=alpha), size=3.5)
  
  # Check which plots to make
  if (any(to.plot %in% "preservation")) {
    
    my.p[["preservation"]] = ggplot2::ggplot(plotData, ggplot2::aes(x=Size, y=Zsum, label=Cols)) +
      my.gp +
      my.theme +
      ggplot2::labs(x ="Module size", y = "Zsummary")
    
  }
  
  if (any(to.plot %in% "median.rank")) {
    
    my.p[["median.rank"]] = ggplot2::ggplot(plotData, ggplot2::aes(x=Size, y=Rank, label=mylabels)) + 
      my.gpa +
      my.theme +
      ggplot2::scale_y_continuous(trans = "reverse", expand = c(0.1,0)) +
      ggplot2::labs(x ="Module size", y = "Median Rank")
    
  }
  
  if (any(to.plot %in% "density")) {
    
    my.p[["density"]] = ggplot2::ggplot(plotData, ggplot2::aes(x=Size, y=Density, label=Cols)) +
      my.gp +
      my.theme +
      ggplot2::labs(x ="Module size", y = "Density") +
      ggplot2::ylim(my.lim)
    
  }
  
  if (any(to.plot %in% "connectivity")) {
    
    my.p[["connectivity"]] = ggplot2::ggplot(plotData, ggplot2::aes(x=Size, y=Connectivity, label=Cols)) +
      my.gp +
      my.theme +
      ggplot2::labs(x ="Module size", y = "Connectivity") +
      ggplot2::ylim(my.lim)
    
  }
  
  # Arrange the plots as needed
  my.p = my.p[to.plot]
  
  # If put the legend in the last plot
  my.p[[length(my.p)]] = my.p[[length(my.p)]] + 
    ggplot2::theme(legend.text = ggplot2::element_text(size = 8), legend.title = ggplot2::element_text(size=10)) +
    ggplot2::guides(fill = "none", alpha="none", color="none")
  
  # If we have more than one, remove the legend from the non-last plots
  if (length(my.p) > 1) {
    my.p[-length(my.p)] = lapply(my.p[-length(my.p)], function(gp){
      gp + ggplot2::theme(legend.position="none")
    })
    gridExtra::grid.arrange(grobs=my.p, ncol=2)
  } else{my.p}
  
}