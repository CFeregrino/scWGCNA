
#' Calculate variable genes
#' @param x an integer
#' @noRd

# If no variable genes are provided

calc.vargenes = function(se.sc, min.cells){
  
  nonex = which(apply(se.sc@assays$RNA@counts, 1, function(x) length(which(x >0))) < my.min.cell)
  # First get rid of non-expressed genes
  se.sc = subset(se.sc, features = rownames(se.sc@assays$RNA)[-nonex])
  
  # Find the variable genes
  se.sc=Seurat::FindVariableFeatures(se.sc, dispersion.cutoff = c(0.25,Inf),
                                     mean.cutoff = c(0,Inf), selection.method = "mvp", assay = "RNA")
  
  return(Seurat::VariableFeatures(object = se.sc, assay = "RNA"))

}








