#' Calculates module eigengenes, average expression, membership and p.values from different sc data.
#' 
#' This function calculates module eigengenes, average expression, membership and p.values, from a Seurat object, using as reference a set of gene modules.
#' @param modules scWGCNA.data or named vector. An scWGCNA.data object, as calculated by run.scWGCNA().
#' @param seurat a Seurat object, which has THE SAME gene names as the ones provided in the modules argument.
#' @return A list with different data
#' @importFrom stats aggregate
#' @export
#' @examples
#' #  We can calculate the eigengenes taking the expression of the single cells
#' eigengenes.sc = scW.eigen(modules = MmLimbE155.scWGCNA, seurat = my.small_MmLimbE155)
#' 

scW.eigen = function(modules, seurat){
  
  if (class(modules) == "list") {
    my.cols = modules[["dynamicCols"]]
  } else {my.cols = modules}
  
  if (is.null(names(my.cols))) {
    stop("Modules should be provided, either as a scWCGNA object calculated from run.scWGCNA.\n
         Alternatively, as a named vector of colors (or other module names), where the names correspond to the gene names")
  }
  
  my.s = Seurat::GetAssayData(seurat)
  
  my.s = t(as.matrix(my.s[names(my.cols),]))
  
  my.me = WGCNA::moduleEigengenes(my.s, colors = my.cols)
  
  my.gmm = as.data.frame(WGCNA::signedKME(my.s, my.me$eigengenes))
  my.mmp = as.data.frame(WGCNA::corPvalueStudent(as.matrix(my.gmm), nrow(my.s)))
  
  my.list = list(Module.Eigengenes = my.me, Gene.Module.Membership = my.gmm, Module.Membership.Pvalue = my.mmp)
  
  return(my.list)
  
}