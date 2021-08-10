#' Foo
#'
#' @param x foo
#'
#' @return fooo
#' @noRd
#'
my.foo=function(x){
  cat("\n\nIMPORTANT NOTE!!!\nYou have run this analysis witht the option less=TRUE. This means that the analysis will try to reduce the number of modules detected, based on their expression pattern. If modules have very similar expression profile (distance < 0.25), they will be merged. Moreover, if a module seems to be highly expressed in only 1-3 cells ( cells expressing >2*(max(expression)/3) ), it will be removed.\n")
  
  my.fcolors = dynamicColors
  
  MEDiss = 1-cor(MEs)
  MEDiss[!lower.tri(MEDiss)] = NA
  
  if (any(apply(MEDiss, 1, function(x) min(x, na.rm = T))<0.25)) {
    
    cat("\nInfo: Some clusters will be merged\n")
    
    MEtree = hclust(as.dist(MEDiss), method = "average")
    plot(MEtree)
    abline(h=0.25,col="red")
    
    my.merge = WGCNA::mergeCloseModules(datExpr, dynamicColors, cutHeight = 0.25, verbose = 3)
    
    raw.MEList = WGCNA::moduleEigengenes(raw.datExpr, colors = merge$colors)
    p.MEList = raw.MEList
    my.fcolors = merge$colors
    
  }
  
  my.oc=which(apply(p.MEList$averageExpr, 2, function(x) length(which(x > 2*(max(x)/3)))) < 3)
  
  if (length(my.oc) > 0) {
    my.oc = substr(colnames(p.MEList$averageExpr[my.oc]),3,nchar(colnames(p.MEList$averageExpr[my.oc])))
    cat(paste0("The following modules were only highly expressed in less than 3 cells :  ",
               paste(my.oc,collapse = ", " ), "\n","They will all be removed"))
    
    my.fcolors[which(my.fcolors %in% my.oc)] = NA
    
  }
  
  if (!identical(dynamicColors,my.fcolors)) {
    my.fMods = as.numeric(as.factor(my.fcolors))
    my.fcolors = WGCNA::labels2colors(my.fMods)
  }
  
  cat("This is the result of the filtering.\nThe color of the modules does not correspond between original and filtered\n")
  
  WGCNA::plotDendroAndColors(geneTree, cbind(dynamicColors,my.fcolors), c("Modules","Filtered"),
                             dendroLabels = NULL,
                             cex.dendroLabels = 0.6,
                             addGuide = TRUE,
                             main = "Gene dendrogram and module colors",
                             guideAll = F)
  
  if (!identical(dynamicColors,my.fcolors)) {
    if (any(is.na(my.fMods))) {
      datExpr = datExpr[,!is.na(my.fMods)]
      my.fMods = my.fMods[!is.na(my.fMods)]
    }
    dynamicMods = my.fMods
    dynamicColors = WGCNA::labels2colors(dynamicMods)
  
      MEList = WGCNA::moduleEigengenes(as.matrix(datExpr), colors = dynamicColors)
      MEs = MEList$eigengenes
      raw.datExpr = p.Wdata@assays$RNA@data[colnames(datExpr),]
      raw.datExpr = t(as.matrix(raw.datExpr))
      raw.MEList = WGCNA::moduleEigengenes(raw.datExpr, colors = dynamicColors)
      xy = lapply(as.list(names(table(dynamicColors))), function(x)
        return(colnames(datExpr)[which(dynamicColors==x)]))
      p.MEList = raw.MEList
      my.adjacency =WGCNA::adjacency(datExpr,type = "signed", power = my.power, corFnc = "bicor")
      TOM=WGCNA::TOMsimilarityFromExpr(datExpr,networkType = "signed", TOMType = "signed",
                                       power = my.power, corType = "bicor")
  
  }
  
  #Make it a distance
  dissTOM = 1-TOM
  
  # Call the hierarchical clustering function
  geneTree = hclust(as.dist(dissTOM), method = "average")
  
  WGCNA::plotDendroAndColors(geneTree, dynamicColors, "Filtered modules",
                             dendroLabels = NULL,
                             cex.dendroLabels = 0.6,
                             addGuide = TRUE,
                             main = "Gene dendrogram and module colors",
                             guideAll = F)
  }

