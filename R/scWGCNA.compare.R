#' Runs a comparative scWGCNA analysis
#' 
#' This function runs a WGCNA analysis adapted for single cells. Based on single-cell or pseudocell data.
#' @param scWGCNA.data The WGCNA data to use as reference for comparative analysis. Output from run.scWGCNA
#' @param test.list List of Seurat objects. The samples to test conservation in. If data is calculated on pseudocells, these should also be pseudocell objects
#' @param test.names Vector of strings. Contains the names for each of the samples in test.
#' @param ortho Data frame. If the test objects and the reference data have different feature names, a data frame should be provided. Only listing 1-to-1 orthologous genes. The species of the reference MUST always be column no. 1
#' @param ortho.sp Vector numeric. Indicating to which species in ortho each test sample belongs to.
#' @param is.pseudocell Logic. Is the analysis based on pseudocells?
#' @param ... Further parameters, passed to WGCNA::modulePreservation
#' @return It returns a new object, main result of WGCNA::modulePreservation
#' @export
#' @importFrom WGCNA bicor
#' @examples
#' # A pre-calculated WGCNA data list
#' class(MmLimbE155.scWGCNA)
#' Mm.scW = MmLimbE155.scWGCNA
#' 
#' # A Seurat object with chicken limb cells
#' my.small_GgLimbHH29
#' 
#' # We calculate pseudocells for this object
#' Gg.ps=calculate.pseudocells(my.small_GgLimbHH29, dims = 1:10)
#' 
#' # Our Seurat objects contain gene names equivalencies in the misc slot.
#' # We use them to biuld a table of orthologous genes
#' my.ortho = merge(my.small_MmLimbE155@@misc$gnames,my.small_GgLimbHH29@@misc$gnames, by = "Gene.name")
#' my.ortho=my.ortho[,2:3]
#' 
#' # Comparative tests
#' MmvGg = scWGNA.compare(Mm.scW, test.list=list(Gg.ps), test.names="Gg", ortho=my.ortho, ortho.sp=2)

scWGNA.compare = function(scWGCNA.data,
                          test.list,
                          test.names,
                          ortho = NULL,
                          ortho.sp = NULL,
                          is.pseudocell=T,
                          ...) {
  
  backup.options = options()
  options(stringsAsFactors = FALSE)
  
  # A list, to keep the genes that would be discarded in the process.
  my.genes = list()
  
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
  
  if(is.pseudocell){
    
    test.list = lapply(test.list, function(s){
      Seurat::GetAssayData(s, slot = "counts", assay="RNA")
    })
    
  } else {test.list = lapply(test.list, function(s){
    Seurat::GetAssayData(s)
  })
  }
  
  my.colfrac = data.frame(table(scWGCNA.data$dynamicCols))
  
  if (!is.null(ortho)) {
    
    #Keep track of which genes are lost
    my.genes[["non-orthologous"]] = data.frame(gene=Expr[which(!Expr %in% ortho[,1])],
               module=scWGCNA.data$dynamicCols[Expr[which(!Expr %in% ortho[,1])]])
    
    # Keep only the ones present as 1-to-1 orthologues
    Expr = Expr[which(Expr %in% ortho[,1])]
    
    # Keep the count on the final numbers, on the table
    my.colfrac$ortho = as.numeric(table(scWGCNA.data$dynamicCols[Expr]))
    
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
  
  #Keep track of which genes are lost
  my.genes[["non-expressed"]] = data.frame(gene=colnames(datExpr[[1]])[nonex],
                                           module=scWGCNA.data$dynamicCols[colnames(datExpr[[1]])[nonex]])
  
  #Take then these genes out of the genes we will use to test.
  if (length(nonex)) {
    Expr = colnames(datExpr[[1]])[-nonex]
  } else {Expr = colnames(datExpr[[1]])}
  
  my.colfrac$expr = as.numeric(table(scWGCNA.data$dynamicCols[Expr]))
  colnames(my.colfrac) = c("module","total","ortho","expr")
  
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
  
  mp$misc = list(
    modulefrac = my.colfrac,
    geneslost = my.genes,
    testnames = test.names
    )

  options(backup.options)
  
  return(mp)
  
}