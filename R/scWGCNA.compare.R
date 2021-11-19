#' Runs a comparative scWGCNA analysis
#' 
#' This function runs a WGCNA analysis adapted for single cells. Based on single-cell or pseudocell data.
#' @param scWGCNA.data Seurat object. The WGCNA data to use as reference for comparative analysis. Output from scWGCNA.report.
#' @param test.list List of seurat objects. The samples to test conservation in. If data is calculated on pseudocells, these should also be pseudocells
#' @param test.names Vector of strings. Contains the names for each of the samples in test.
#' @param ortho Data frame. As many columns as species used only listing 1-to-1 orthologous genes. The species of the reference MUST always be column no. 1
#' @param ortho.sp Vector numeric. Indicating to which species in ortho each test sample belongs to.
#' @param is.pseudocell Is it?
#' @return No inline output. It saves an html report, as well as an object with the resulting WGCNA comparison indices. They are both named using the project name and the date.
#' @export
#' @importFrom WGCNA bicor
#' 

scWGNA.compare = function(scWGCNA.data,
                          test.list,
                          test.names,
                          ortho = NULL,
                          ortho.sp = NULL,
                          is.pseudocell=T,
                          ...) {
  
  backup.options = options()
  options(stringsAsFactors = FALSE)
  
  # a list to keep our expression matrices
  datExpr = list()
  
  # Go to the WGCNA data, take all the module genes
  Expr = colnames(scWGCNA.data[["expression"]])
  
  if (any(unlist(lapply(test.list, function(s){
      length(which(Expr %in% rownames(s)))
    })) == 0) &
    is.null(ortho)) {
    
    return(message("It seems like at least one dataset from the test list doesn't share feature names with the reference data.\n
               Please make sure they do, or include an orthology table if you're working with different species."))
    
  }
  
  if (!is.null(ortho.sp)) {
    if (length(ortho.sp) != length(test.list)) {
      return(message("Please make sure that ortho.sp is of the same length as the amount of samples (including the reference)"))
      }
  }
  
  if(is.null(ortho)){
    
    ortho = data.frame(x=Expr, y=Expr)
    ortho.sp = c(1,rep(2,length(test.list)))
    
  }
  
  if(is.pseudocell){
    
    test.list = lapply(test.list, function(s){
      GetAssayData(s, slot = "counts", assay="RNA")
    })
    
  } else {test.list = lapply(test.list, function(s){
    GetAssayData(s)
  })
  }
  
  if (!is.null(ortho)) {
    
    # Keep only the ones present as 1-to-1 orthologues
    Expr = Expr[which(Expr %in% ortho[,1])]
    
    #Change the ortho rownames to the reference species, to be able to index it
    rownames(ortho) = ortho[,1]
    #A subset of the ortho, using only the genes in the modules
    Expr = ortho[Expr,,drop=F]
    
    #Each of the test sets
    for (t in 1:length(test.list)) {
      #Take only genes that are present in all of our expression matrices.
      Expr=Expr[which(Expr[,ortho.sp[t]] %in% rownames(test.list[[t]])),,drop=F]
    }
    
    # Retrieve gene id's from the reference species.
    Expr=rownames(Expr)
    
    #Each of the test sets
    for (t in 1:length(test.list)) {
      #Subset the ortho for the test speices, and genes we just got
      x = ortho[,ortho.sp[t]][which(ortho[,1] %in% Expr)]
      #Take the expression matrix
      datExpr[[t]] = t(as.matrix( test.list[[t]][x,] ))
      #Translate the names of the genes, to the reference species
      colnames(datExpr[[t]]) = ortho[,1][match(colnames(datExpr[[t]]), ortho[,ortho.sp[t]])]
    }  
    
  } else{
    
    #Each of the test sets
    for (t in 1:length(test.list)) {
      #Take only genes that are present in all of our expression matrices.
      Expr=Expr[which(Expr %in% rownames(test.list[[t]]))]
    }
    
  }
  
  
  #There might be genes present in our samples, but not expressed, or not variable. Take them all
  nonex=c()
  for (i in datExpr) {
    nonex = c(nonex,which(!WGCNA::goodGenes(i)))
  }
  nonex = unique(nonex)
  
  #Take then these genes out of the genes we will use to test.
  if (length(nonex)) {
    Expr = colnames(datExpr[[1]])[-nonex]
  } else {Expr = colnames(datExpr[[1]])}
  
  #Take only those genes in the expression matrix
  for (t in 1:length(datExpr)) {
    datExpr[[t]] = datExpr[[t]][,Expr]
  }
  
  # Here we asign the reference data
  refdat = scWGCNA.data[["expression"]][,Expr]
  
  # The reference colors
  refcol = scWGCNA.data$dynamicCols[Expr]
  
  # Objects with both the reference and the test datasets
  multiExpr = list(Reference = list(data = refdat))
  
  for (t in 1:length(datExpr)) {
    multiExpr[[t+1]] = list(data = datExpr[[t]])
  }
  
  names(multiExpr) = c("Reference", test.names)
  
  multiColor = list(Reference = refcol)
  
  my.dummy=utils::capture.output({
    mp = WGCNA::modulePreservation(multiExpr,
                                   multiColor,
                                   corFnc = "bicor",
                                   maxGoldModuleSize = 300,
                                   networkType = "signed",
                                   referenceNetworks = 1,
                                   nPermutations = 20,
                                   randomSeed = 42,
                                   quickCor = 0,
                                   verbose = 0,
                                   savePermutedStatistics = F,
                                   ...)
  })
  
  return(mp)
  
}