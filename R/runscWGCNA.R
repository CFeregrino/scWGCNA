#' Runs a semi-automatic, iterative scWGCNA analysis
#' 
#' This function runs our semi-automatic single-cell WGCNA analysis. It runs in an iterative way. Based on single-cell or pseudocell data.
#' @param p.cells Seurat object. The expression data used to run the co-expression analysis. Can be pseudocell or single-cell data but pseudocells are recommended.
#' @param s.cells Seurat object. The single cell data, if running on single cell data already, please repeat the argument.
#' @param idents Variable. Are certain clusters to be used? Please use group identities and not cell names.
#' @param features Variable. The features to be used for the analysis. Default is F, which makes the script calculate variable genes.
#' @param is.pseudocell Logical. Is the main data pseudocell data? Default is T
#' @param min.cells Numeric. The minimum cells in which genes need to be expressed, to be considered for variable genes calculation. Default is 10
#' @param less Logical. In case you are getting too many small modules, use this option as T
#' @param g.names Data frame. If you're using gene IDs and no symbols, you might wanna provide a list of gene names for plotting. Two columns: 1= ids present in expression matrix, 2= names to appear in plots. Rownames= same as 1st row
#' @return A list object with the resulting WGCNA data. 
#' @export
#' @importFrom WGCNA bicor
#' @examples
#' # Calculate pseudocells
#' ps.pbmc_small=calculate.pseudocells(SeuratObject::pbmc_small, dims = 1:10)
#' 
#' # Use pseudocells and single cells to calculate WGCNA
#' scWGNA.pbmc_small = run.scWGCNA(p.cells = ps.pbmc_small, s.cells = SeuratObject::pbmc_small)
#' 

run.scWGCNA = function(p.cells,
                       s.cells,
                       idents,
                       features,
                       is.pseudocell=T,
                       min.cells=10,
                       less=F,
                       g.names){
  
  backup.options = options()
  options(stringsAsFactors = FALSE)
  
  #Here we will keep track of the trees, in case we want to plot 'em later
  my.trees=list()
  
  # If we need to subset the data
  if (!missing(idents)) {
    s.Wdata = subset(s.cells, idents=idents)
  } else s.Wdata = s.cells
  
  # If we are using IDs or symbols
  if (missing(g.names)) {
    gnames= data.frame(x=rownames(p.cells),y=rownames(p.cells), row.names = rownames(p.cells))
  }  else {gnames = g.names; rownames(gnames) = gnames[1,]}
  
  nonex = which(apply(s.Wdata@assays$RNA@counts, 1, function(x) length(which(x >0))) < min.cells)
  
  # If no variable genes are provided
  if (missing(features)) {
    Expr = calc.vargenes(s.Wdata, min.cells)
    } else { Expr = features }

  if (is.pseudocell==T) {
    
    datExpr=p.cells@assays$RNA@counts[Expr,]
    datExpr = datExpr[which(Matrix::rowSums(datExpr)>0),]
    if (length(which(Matrix::rowSums(datExpr)==0))<1) {
      print(paste0("The following variable genes were not found expressed in the pseudocell object:  ",
                   names(which(Matrix::rowSums(datExpr)==0))))
    }
    Expr = rownames(datExpr)
    
  } else{datExpr=s.Wdata@assays$RNA@data[Expr,]}
  
  # Check the length
  print(paste0("We have ", length(Expr), " genes in the variable genes object"))
  
  # Check the size and transform
  datExpr = t(as.matrix(datExpr))
  
  # Choose a set of soft-thresholding powers
  my.power = choose.sfttp(datExpr)
  
  my.Clnumber = 20
  change = 0
  genesets=list()
  nonsig = 1
  
  while(nonsig != 0) {
    
    my.adjacency = WGCNA::adjacency(datExpr,type = "signed", power = my.power, corFnc = "bicor")
    
    # Turn adjacency into topological overlap (high overlap if they share the same "neighborhood")
    TOM=WGCNA::TOMsimilarityFromExpr(datExpr,networkType = "signed", TOMType = "signed",
                                     power = my.power, corType = "bicor")
    
    #Put the names in the tree
    colnames(TOM) <- gnames[colnames(datExpr),2]
    
    rownames(TOM) <- gnames[colnames(datExpr),2]
    
    #Make it a distance
    dissTOM = 1-TOM
    
    # Call the hierarchical clustering function
    geneTree = hclust(as.dist(dissTOM), method = "average")
    
    # Here I calculate the cutting height. Using the same formula and approach that the WGCNA package uses for the automatic function
    nMerge = length(geneTree$height) # The whole height of the tree
    refQuantile = 0.05 # What's the quantile that we want to exclude
    refMerge = round(nMerge * refQuantile) 
    refHeight = geneTree$height[refMerge]
    cutheight = signif(0.99 * (max(geneTree$height) - refHeight) + refHeight,4)
    
    # We construct THE TABLE that will help us make decisions
    # Min cluster sizes, from 7 to 30
    x=seq(7,30,1)
    # The height, up and down from the calculated height. We expect some "No module detected"
    y=seq(cutheight-0.0005,cutheight + 0.0005,0.0001)
    
    # The actual dataframe
    w=data.frame()
    # Populate, with i=min cluster size, j=cutting height, z=total number of clusters, z.1.=what's the first cluster? 0 is grey 1 is something else,
    # z.1.'=what's the size of the first cluster?
    for (i in x) {
      for (j in y) {
        sink("aux")
        z=table(dynamicTreeCut::cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = i, deepSplit = T, cutHeight = j, verbose = 0))
        sink(NULL)
        v=data.frame(i,j,dim(z),names(z[1]),unname(z[1]))
        w=rbind(w,v)
      }
    }
    
    # The height is then the one where we have the least number of genes in the first cluster
    my.height = w$j[which(w$unname.z.1..==min(w$unname.z.1..))]
    
    
    # Since different heights can give us the minimum grey size, we chose the computed height, if present, or the highest one.
    if (cutheight %in% my.height) {
      my.Clsize = w[which(w$j == cutheight),]
    } else { my.Clsize = w[which(w$j == max(my.height)),] }
    
    
    # This is to know, if we're looking for a minimum of cluster numbers
    # If we still have a lot of genes, we don't want to limit the number of clusters
    # if ( ((dim(datExpr)[2]) / length(Expr)) > 0.6 ) {change = 0}
    
    # If this is the first iteration after 0.6 of the genes are gone, we assign the number of clusters (and an extra for the grey in the case)
    if (change == 2) {
      my.Clnumber = length(table(dynamicColors)) + 1
    }
    # Count another iteration
    change = change + 1
    
    # If we don't have a gray cluster anymore, then we subset for those rows, and set a new height. ONLY if we get the same amount of clusters!
    if (any(w$names.z.1.. == 1)) { #any combination gives us no grey
      
      my.Clsize = w[which(w$names.z.1.. == 1),] # Take all combinations that gives us no grey
      
      if (any(my.Clsize$dim.z. >= (my.Clnumber -1) )) { # If there is any giving us the determined amount or more
        my.Clsize = my.Clsize[which(my.Clsize$dim.z. >= (my.Clnumber - 1) ),,drop=F] # Subset for those
      } else { my.Clsize = my.Clsize[which(my.Clsize$dim.z. == max(my.Clsize$dim.z.)),,drop=F] } # Or for the highest
      
      # Take the ones with the smallest number of clusters
      my.Clsize = my.Clsize[which(my.Clsize$dim.z. == min(my.Clsize$dim.z.)),,drop=F]
      # Take the one with the highest min cluster size
      my.Clsize = my.Clsize[which(my.Clsize$i == max(my.Clsize$i)),,drop=F]
      
      if (cutheight %in% my.Clsize$j) { # if original computed height is in,
        my.height = cutheight # take it
      } else { my.height = max(my.Clsize$j) } # Otherwise, the highest height
      
      my.Clsize = max(my.Clsize$i)
      
    }
    
    # Subset the table again, for those sizes that will gives the same number of clusters or more. IF NONE, use the highest number
    if (!any(w$names.z.1.. == 1)){
      if (any(my.Clsize$dim.z. >= my.Clnumber)) {
        my.Clsize = my.Clsize[which(my.Clsize$dim.z. >= my.Clnumber),,drop=F]
      } else {
        my.Clsize = my.Clsize[which(my.Clsize$dim.z. == max(my.Clsize$dim.z.)),,drop=F]}
      
      # Take the ones with the smallest number of clusters
      my.Clsize = my.Clsize[which(my.Clsize$dim.z. == min(my.Clsize$dim.z.)),,drop=F]
      # Take the one with the highest min cluster size
      my.Clsize = my.Clsize[which(my.Clsize$i == max(my.Clsize$i)),,drop=F]
      
      if (cutheight %in% my.Clsize$j) { # if original computed height is in,
        my.height = cutheight # take it
      } else { my.height = max(my.Clsize$j) } # Otherwise, the highest height
      
      my.Clsize = max(my.Clsize$i)
    }
    
    # If we still have more than 60% of the genes, we just use the min size of 15 regardles
    
    if ( change < 3 ) {
      my.Clsize = 15
      my.height = cutheight
      # if (my.less == T) {####
      #   my.Clsize = w[which(w$j == my.height),]####
      #   my.Clsize = my.Clsize$i[which.min(my.Clsize$dim.z.)]####
      # } ####
    }
    
    print(paste("my.height: ",my.height," .... my.Clsize: ", my.Clsize))
    
    dynamicMods = dynamicTreeCut::cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = my.Clsize, deepSplit = T, cutHeight = my.height)
    
    #dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, minClusterSize = minModuleSize)
    
    table(dynamicMods)
    
    # Convert numeric lables into colors
    dynamicColors = WGCNA::labels2colors(dynamicMods)
    table(dynamicColors)
    
    #Save, to keep track of the trees
    my.trees[[length(my.trees)+1]] = list(geneTree, dynamicColors)

    # Calculate eigengenes
    MEList = WGCNA::moduleEigengenes(as.matrix(datExpr), colors = dynamicColors)
    
    MEs = MEList$eigengenes
    
    # Calculate the module membership
    geneModuleMembership = as.data.frame(WGCNA::signedKME(datExpr, MEs))
    MMPvalue = as.data.frame(WGCNA::corPvalueStudent(as.matrix(geneModuleMembership), nrow(datExpr)))
    
    # We're gonna make a list, where we keep the genes that are significantly associated with each module
    x=c()
    xy=list()
    
    # We also need a vector with all the dynamic colors
    dcols = 1:length(levels(as.factor(dynamicColors)))
    
    # Getting rid of the grey module
    grey.genes = length(which(dynamicColors == "grey"))
    if (any(levels(as.factor(dynamicColors)) == "grey")) {
      dcols = dcols[-which(levels(as.factor(dynamicColors)) == "grey")]
    }
    
    # Run the loop to get the genes
    for (i in dcols) {
      modGenes = rownames(MMPvalue)[which(dynamicColors==levels(as.factor(dynamicColors))[i] & MMPvalue[,i]<0.01)]
      x=c(x,modGenes)
      xy[[i]]=modGenes
      #print(paste0(levels(as.factor(dynamicColors))[i]," ",length(modGenes),
      #" of ", length(which(dynamicColors==levels(as.factor(dynamicColors))[i]))))
      #print(gnames[modGenes,2])
    }
    
    # Make a new list, where we keep ALL the gens thar are left from the iteration, that will be used to make the new object. To keep track
    genesets[[length(genesets)+1]] = colnames(datExpr)
    
    # Give me a message saying how many genes are gone this time
    cat( paste0( grey.genes, " genes not assigned to any module.", '\n',
                 length(which(!(colnames(datExpr)%in%x))) - grey.genes, " genes excluded due to significance."))
    # Save this also, cause if it's 0 then we stop the whole thing
    nonsig = length(which(!(colnames(datExpr)%in%x)))
    
    # If it ain't 0, subset the dynamic colors and the expression data
    if (length(which(!(colnames(datExpr)%in%x))) != 0) {
      dynamicColors=dynamicColors[-which(!(colnames(datExpr)%in%x))]
      datExpr=datExpr[,-(which(!(colnames(datExpr)%in%x)))]
    }
    
    if (change == 2 & less == T) {
      
      cat("\n\nIMPORTANT NOTE!!!\nYou have run this analysis witht the option less=TRUE. This means that the analysis will try to reduce the number of modules detected, based on their expression pattern. If modules have very similar expression profile (distance < 0.25), they will be merged. Moreover, if a module seems to be highly expressed in only 1-3 cells ( cells expressing >2*(max(expression)/3) ), it will be removed.\n")
      
      my.filtered = FilterMods_int(dynamicColors=dynamicColors, s.Wdata=s.Wdata, datExpr=datExpr,
                                   geneTree=geneTree, my.power = my.power)
      par(mfrow=c(1,1))
      datExpr = my.filtered[[1]]
      dynamicColors = my.filtered[[2]]
      
    }
    
  }
  
  #We prepare the outputs
  
  my.memberships=list()
  
  my.cols = as.factor(dynamicColors)
  
  for (m in 1:length(levels(my.cols))) {
    
    my.mem = cbind(geneModuleMembership[my.cols==levels(my.cols)[m],m,drop=F],
                   MMPvalue[my.cols==levels(my.cols)[m],m,drop=F])
    
    colnames(my.mem) = c("Membership", "p.val")
    
    my.memberships[[paste0(m,"_",levels(my.cols)[m])]] = my.mem
      
  }
  
  s.MEList = MEList
  
  raw.datExpr = s.Wdata@assays$RNA@data[colnames(datExpr),]
  
  raw.datExpr = t(as.matrix(raw.datExpr))
  
  raw.MEList = WGCNA::moduleEigengenes(raw.datExpr, colors = dynamicColors)
  
  s.MEList = raw.MEList
  
  names(dynamicColors) = colnames(datExpr)
  
  WGCNA_data = list()
  WGCNA_data[["modules"]] = my.memberships
  WGCNA_data[["expression"]] = datExpr
  WGCNA_data[["dynamicMods"]] = dynamicMods
  WGCNA_data[["dynamicCols"]] = dynamicColors
  WGCNA_data[["MEList"]] = MEList
  WGCNA_data[["MEs"]] = MEs
  WGCNA_data[["module.genes"]] = xy
  WGCNA_data[["sc.MEList"]] = s.MEList
  WGCNA_data[["genesets.history"]]= genesets
  WGCNA_data[["moduletrees.history"]] = my.trees
  WGCNA_data[["TOM"]]= TOM
  WGCNA_data[["adjacency"]]= my.adjacency
  
  return(WGCNA_data)
  
  }





