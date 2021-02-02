#' Calculates pseudocells from a Seurat object
#' 
#' This function calculates pseudocells from a Seurat object, based on pre-calculated cell clusters and dimentionality reduction.
#' @name calculate.pseudocells
#' @usage calculate.pseudocells(seurat, seeds=0.2, nn = 10, reduction = "pca", dims = 1:20, features =NULL, cells= NULL)
#' @param seurat The seurat object, with pre-computed PCA or other reductions, and the relevant clustering as IDs
#' @param seeds The proportion of cells to be used as seeds. Alternatively, a string with the name of the seeds to use. Numeric between 0.1 and 0.9 or string. Default 0.2
#' @param nn Number of nearest neighbors to compute and use for pseudocell aggregation. Default 10
#' @param reduction The name of the reduction to use. Should be present in the @@reductions slot of the seurat object. Default is "pca"
#' @param dims The relevant dimensions that will be used to compute nearest neighbors. Default 1:20
#' @param features The features to be used. Takes a string of feature names as present in the expression matrices. Defaults to NULL, which will use all the genes.
#' @param idents The clusters or identities to be used. Takes a string of identities as present in the @@active.ident slot. Defaults to NULL, which will use all the identities.
#' @return A seurat object of aggregated pseudocells. With average expression. The slot misc contains the pseudocells dataframe, with each original cell and its assigned pseudocell
#' @export
#' @examples
#' ps.pmbc_small=calculate.pseudocells(Seurat::pbmc_small, dims = 1:10)

library("Seurat")

calculate.pseudocells <- function(seurat, seeds=0.2, nn = 10, reduction = "pca", dims = 1:20, features =NULL, cells= NULL) {

  if (!is.null(cells)) {
    seurat = subset(seurat, idents=cells)
  }
  
  if (!is.null(features)) {
    seurat = subset(seurat, features=features)
  }
  
  seurat = Seurat::FindNeighbors(seurat,
                         reduction = reduction,
                         dims = dims,
                         compute.SNN = F,
                         k.param = 10)
  
  nn.matrix <- as.matrix(seurat@graphs[[paste0(Seurat::DefaultAssay(seurat),"_nn")]])
  
  my.seeds = list()
  seeds.count = c()
  
  if (is.numeric(seeds)) {
    
    message("Choosing seeds")
    
    for (i in 1:50) {
      
        seed.set = c()
        for (cluster in levels(Seurat::Idents(seurat)) ) {
          n.seeds = floor(table(Seurat::Idents(seurat))[[cluster]]/(1/seeds))
          seed.set=c(seed.set, sample(rownames(seurat@meta.data)[Seurat::Idents(seurat) == cluster], n.seeds) )
      }
      
      my.seeds[[i]] = seed.set
      rm(seed.set)
      
      seeds.count = c(seeds.count, length(which(colSums(nn.matrix[my.seeds[[i]],]) > 0)) )
      
    }
  
    seeds = my.seeds[[which.max(seeds.count)]]
    
  }
  
  nn.matrix <- nn.matrix[seeds, ]
  nn.matrix <- data.frame(nn.matrix[,colSums(nn.matrix) > 0 ])
  
  message(ncol(nn.matrix)," out of ", ncol(seurat)," Cells will be agreggated into ",nrow(nn.matrix),"Pseudocells")
  
  my.pseudocells = data.frame(pseudocell = rep("00",nrow(seurat@meta.data)), row.names = rownames(seurat@meta.data))
  
  my.pseudocells$pseudocell = as.character(my.pseudocells$pseudocell)
  
  my.pseudocells[seeds,1] = seeds
  
  nn.matrix = nn.matrix[,-which(colnames(nn.matrix) %in% seeds)]
  
  remaining.cells = colnames(nn.matrix)
  
  message("Assining pseudocells")
  
  while (length(remaining.cells) > 0) {
    
    for (s in rownames(nn.matrix)[order(rowSums(nn.matrix[,remaining.cells, drop = F]))]) {
      
      if ( sum(nn.matrix[s, remaining.cells]) > 0 & length(remaining.cells) > 0) {
        
        c = sample(which(nn.matrix[s, remaining.cells] == 1),1)
        
        my.pseudocells[remaining.cells[c],1] = s
        
        remaining.cells = remaining.cells[-c]
        
      }
      
      if (length(remaining.cells) == 0) {break}
      
    }
    
  }
  
  
  ps.seurat = subset(seurat, cells = rownames(my.pseudocells)[my.pseudocells$pseudocell != "00"])
  
  ps.seurat = Seurat::AddMetaData(object = ps.seurat, metadata = my.pseudocells, col.name = "pseudo.ident")
  
  ps.seurat = Seurat::SetIdent(ps.seurat, value = "pseudo.ident")
  
  message("Aggregating pseudocell expression")
  
  ps.seurat = Seurat::AverageExpression(object = ps.seurat,
                                return.seurat = T,
                                verbose = F)
  
  ps.seurat@meta.data$orig.cluster = Seurat::Idents(seurat)[match(rownames(ps.seurat@meta.data), rownames(seurat@meta.data))]
  
  ps.seurat@misc[["peudocells"]] = my.pseudocells
  
  return(ps.seurat)

}


