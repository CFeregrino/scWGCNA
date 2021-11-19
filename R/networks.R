
#' Calculates the graphic networks of scWGCNA data
#' 
#' This function calculates the graphic networks of each module in an scWGCNA list object. Returns an updated object.
#' @param scWGCNA.data scWGCNA.data. An scWGCNA.data object, as calculated by run.scWGCNA().
#' @return An scWGCNA.data object, updated with the slot "networks"
#' @importFrom stats aggregate
#' @export
#' @examples
#' # calculate networks
#' MmLimbE155.scWGCNA = scWGCNA.networks(MmLimbE155.scWGCNA)
#' 

scWGCNA.networks = function(scWGCNA.data){
  
  geneModuleMembership = as.data.frame(WGCNA::signedKME(scWGCNA.data[["expression"]], scWGCNA.data[["MEs"]]))
  my.cols = levels(as.factor(scWGCNA.data[["dynamicCols"]]))
  
  
  
  scWGCNA.data[["networks"]] = lapply(as.list(my.cols), function(col){
    mod=is.finite(match(scWGCNA.data[["dynamicCols"]], col))
    mynetwork = suppressWarnings(reshape::melt(scWGCNA.data[["TOM"]][mod,mod]))
    mynetwork = mynetwork[!mynetwork$X1 == mynetwork$X2,]
    colnames(mynetwork) = c("fromNode", "toNode", "weight")
    mynetwork$fromNode = as.character(mynetwork$fromNode)
    mynetwork$toNode = as.character(mynetwork$toNode)
    
    # We get rid of all the edges that have a very small weight. The cutoff set to keep at least ONE edge on the nodes
    x = max( c( min(aggregate(mynetwork$weight, by = list(mynetwork$fromNode), max)$x),
                min(aggregate(mynetwork$weight, by = list(mynetwork$toNode), max)$x) ) )
    mynetwork=mynetwork[-which(mynetwork$weight < x),]
    
    # We rescale the weights so that we have them from 0 to 1
    mynetwork$weight01 = GGally::rescale01(mynetwork[,3]) + 0.1
    # And then multiply them by 2, to give the heavy edges a 2 thickness, we ad 0.2 to not have completelly invisible edges
    mynetwork$weight02 = (GGally::rescale01(mynetwork[,3]) * 2) + 0.2
    # Convert this edgelist into a network object
    mynet = network::network(mynetwork, matrix.type="edgelist", ignore.eval=F)
    
    mynodes = geneModuleMembership[levels(as.factor(mynetwork$fromNode)),which(my.cols == col), drop=F]
    
    colnames(mynodes) <- "membership"
    
    mynodes$membership01 = GGally::rescale01(mynodes$membership)
    
    # We multiply by 30 to give a good range of sizes, plus 1 to avoid innexistant nodes
    mynodes$membership01= (mynodes$membership01*30)+1
    
    mynet = network::set.vertex.attribute(mynet, "membership", 
                                          mynodes[network::network.vertex.names(mynet),
                                                  "membership"])
    
    mynet = network::set.vertex.attribute(mynet, "membership01", 
                                          mynodes[network::network.vertex.names(mynet),
                                                  "membership01"])
    
    mynet
    })
  
  return(scWGCNA.data)
  
}
