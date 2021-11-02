#' Chooses a soft-tresholding power
#' 
#' This function examines the data to find a suitable soft-thresholding power. It also helps to know if there is a scale-free topology.
#' @param datexpr A matrix of gene expression. Thi SHOULD BE the expression data normally stores under @@assays$RNA@@data for single cell objects or ´@assays$RNA@counts´ for pseudocells. This has to contain ONLY the genes used for analysis.
#' @return A number representing the soft-thresholding power with the highest scale-free topology fit index. And a message about it
#' @importFrom graphics abline text
#' @export
#' @examples
#' ps.pbmc_small=choose.sfttp(SeuratObject::pbmc_small@@assays$RNA@@data)

choose.sfttp = function(datexpr){
  
  # Choose a set of soft-thresholding powers
  powers = c(c(1:10), seq(from = 12, to=30, by=2))
  
  # Call the network topology analysis function
  sft = WGCNA::pickSoftThreshold(datexpr, powerVector = powers, verbose = 0, networkType = "signed", corFnc = "bicor")
  
  # Plot of the scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,col="red")
  which((-sign(sft$fitIndices[,3])*sft$fitIndices[,2]) > 0.9)
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.85,col="red")
  
  # These are the scale-free topology indexes
  indexes = (-sign(sft$fitIndices[,3])*sft$fitIndices[,2])
  
  # If we don't have aany index above 0.75, we stop the script
  if ( !any( indexes > 0.75 ) ) {
    print("The scale-free topology index didn't reach 0.75 with any of the chosen powers, please consider changing the set of genes or cells provided")
    
  }
  
  if ( !any( indexes > 0.4 ) ) {
    print("The scale-free topology index didn't reach 0.4 with any of the chosen powers, the process will not continue, since it might fail")
    stop("The scale-free topology index didn't reach 0.4 with any of the chosen powers, the process will not continue, since it might fail")
  }
  
  # Take the smalles power that gives us an index over 0.9, or the highest index if we don't reach 0.9
  if ( any( indexes > 0.9 ) ) {
    my.power = sft$fitIndices$Power[min(which(indexes > 0.9))]
  } else { my.power = sft$fitIndices$Power[which.max(indexes)] }
  
  print(paste0("Or power is ", my.power))
  
  return(my.power)
}
