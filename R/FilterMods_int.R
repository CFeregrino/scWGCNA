#' Filtering modules internally
#'
#' @param dynamicColors the dynamic colors
#' @param p.Wdata the single-cell Seurat object
#' @param datExpr the expression data used to calculate 
#' @param geneTree the gene tree from the iteration we are filtering
#' @param my.power the power used to construct the trees
#'
#' @return a list with a new datExpr and dynamic colors
#' @noRd
#'

FilterMods_int = function(dynamicColors, p.Wdata, datExpr, geneTree, my.power){
  # Take the colors, in case they don't change
  my.fcolors = dynamicColors
  # Take the single-cell expression, in order to check for similar expression
  raw.datExpr = p.Wdata@assays$RNA@data[colnames(datExpr),]
  # Flip it like it's hot
  raw.datExpr = t(as.matrix(raw.datExpr))
  # Calculate eigengenes from the sc data, a distance matrix, and delete half of it.
  p.MEList = WGCNA::moduleEigengenes(raw.datExpr, colors = dynamicColors)
  MEDiss = 1-cor(p.MEList$eigengenes)
  MEDiss[!lower.tri(MEDiss)] = NA
  
  # Check if there are any modules with a distance <0.25 / correlation >0.75
  if (any(apply(MEDiss, 1, function(x) min(x, na.rm = T))<0.25)) {
    # If so, start the merging
    cat("\nInfo: Some clusters will be merged\n")
    # Calculate a tree of modules
    MEtree = hclust(as.dist(MEDiss), method = "average")
    # # Plotit
    # plot(MEtree)
    # abline(h=0.25,col="red")
    # par(mfrow=c(1,1))
    # Merge the modules, we only want the colors
    my.merge = WGCNA::mergeCloseModules(raw.datExpr, dynamicColors, cutHeight = 0.25, verbose = 3)$colors
    # Calculate the eigengenes again, using the new colors
    p.MEList = WGCNA::moduleEigengenes(raw.datExpr, colors = my.merge)
    # Modify the colors
    my.fcolors = my.merge
  }
  
  # Check which modules have a high expression > 2*(max(x)/3) in 1-3 cells
  my.oc=which(apply(p.MEList$averageExpr, 2, function(x) length(which(x > 2*(max(x)/3)))) < 3)
  if (length(my.oc) > 0) {
    # If we have any, get the names of those modules
    my.oc = substr(colnames(p.MEList$averageExpr[my.oc]),3,nchar(colnames(p.MEList$averageExpr[my.oc])))
    cat(paste0("The following modules were only highly expressed in less than 3 cells :  ",
               paste(my.oc,collapse = ", " ), "\n","They will all be removed\n"))
    # And NA them in the modified colors
    my.fcolors[which(my.fcolors %in% my.oc)] = NA
    
  }
  
  # In case we have modified the colors
  if (!identical(dynamicColors,my.fcolors)) {
    # Make new mods, and new colors (these have the same length, NA =  grey)
    my.fMods = as.numeric(as.factor(my.fcolors))
    my.fcolors = WGCNA::labels2colors(my.fMods)
    # Make a press release
    cat("\nThis is the result of the filtering.\nThe color of the modules does not correspond between original and filtered\n")
    # Calculate a new tree (remember, this has already been filtered in the main function)
    TOM=WGCNA::TOMsimilarityFromExpr(datExpr,networkType = "signed", TOMType = "signed", power = my.power, corType = "bicor")
    geneTree = hclust(as.dist(1-TOM), method = "average")
    # Show the new modules, relative to the old ones
    # WGCNA::plotDendroAndColors(geneTree, cbind(dynamicColors,my.fcolors), c("Modules","Filtered"),
    #                            dendroLabels = NULL,
    #                            cex.dendroLabels = 0.6,
    #                            addGuide = TRUE,
    #                            main = "Gene dendrogram and module colors",
    #                            guideAll = F)
    # par(mfrow=c(1,1))
    # Check if we removed any module, cause then we need to filter
    if (any(is.na(my.fMods))) {
      # If so, remove the genes that we NAed, from the expression matrix
      datExpr = datExpr[,!is.na(my.fMods)]
      # Remove them also from the Mods, and make new colors
      my.fMods = my.fMods[!is.na(my.fMods)]
      my.fcolors = WGCNA::labels2colors(my.fMods)
    }
    # Then, subsistute the colors for the modified colors: merged, filtered, or both.
    dynamicColors = my.fcolors
  }
  
  # Return a list, expression 1, colors 2
  return(list(datExpr, dynamicColors))
  }