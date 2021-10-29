

scWGCNA.networks = function(scWGCNA.data){
  
  geneModuleMembership = as.data.frame(WGCNA::signedKME(scWGCNA.data[["datExpr"]], scWGCNA.data[["MEs"]]))
  my.cols = levels(as.factor(scWGCNA.data[["dynamicCols"]]))
  
  
  
  scWGCNA.data[["networks"]] = lapply(list(my.cols), function(col){
    mod=is.finite(match(scWGCNA.data[["dynamicCols"]], col))
    mynetwork = reshape::melt(WGCNA.data[["TOM"]][mod,mod])
    colnames(mynetwork) = c("fromNode", "toNode", "weight")
    
    # We get rid of all the edges that have a very small weight. The cutoff set to keep at least ONE edge on the nodes
    x = max( c( min(aggregate(mynetwork$weight, by = list(mynetwork$fromNode), max)$x),
                min(aggregate(mynetwork$weight, by = list(mynetwork$toNode), max)$x) ) )
    mynetwork=mynetwork[-which(mynetwork$weight < x),]
    
    # We rescale the weights so that we have them from 0 to 1
    mynetwork$weight01 = GGally::rescale01(mynetwork[,3])
    # And then multiply them by 2, to give the heavy edges a 2 thickness, we ad 0.2 to not have completelly invisible edges
    mynetwork[,3] = (GGally::rescale01(mynetwork[,3]) * 2) + 0.2
    # Convert this edgelist into a network object
    mynet = network::network(mynetwork, matrix.type="edgelist", ignore.eval=F)
    
    mynodes = geneModuleMembership[levels(mynetwork$fromNode),which(my.cols == col), drop=F]
    
    colnames(mynodes) <- "membership"
    
    GGally::rescale01(mynodes)
    
    # We multiply by 30 to give a good range of sizes, plus 1 to avoid innexistant nodes
    mynodes$membership= (mynodes$membership*30)+1
    
    return(network::set.vertex.attribute(mynet, "membership", 
                                         mynodes[network::network.vertex.names(mynet),
                                                 "membership"]))
    })
  
  return(scWGCNA.data)
  
}
