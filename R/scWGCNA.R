#' Runs a semi-automatic, iterative scWGCNA analysis
#' 
#' This function runs our semi-automatic single-cell WGCNA analysis. It runs in an iterative way. Based on single-cell or pseudocell data.
#' @name scWGCNA.report
#' @usage scWGNA.report(data,sc.data,gene.names,project.name,sp=Mm,cells=F,features=F,reduction="tsne",dir="./",is.pseudocell=T,GO=T,interactive=T)
#' @param data Seurat object. The expression data used to run the co-expression analysis. Can be pseudocell or single-cell data but pseudocells are recommended.
#' @param sc.data Seurat object. The single cell data, if running on single cell data already, please repeat the argument.
#' @param gene.names Data frame. Two columns: 1= ids present in expression matrix, 2= names to appear in plots. Can be the same, but the two columns are necessary
#' @param project.name String. The name of the project, will be used to name the saved objects and the report.
#' @param sp String. The species used for the analysis in string form. Only two letters in the format Sp. Defaults to "Mm". Check that package "org.Sp.eg.db" exists
#' @param cells Variable. Are certain clusters to be used? Please use group identities and not cell names
#' @param features Variable. The features to be used for the analysis. Default is F, which makes the script calculate variable genes
#' @param reduction String. The reduction to use for plots. It must be present in the @@reductions slot of sc.data. Default is "tsne"
#' @param dir String. The directory where the output will be saved. It defaults to the working directory
#' @param is.pseudocell Logical. Is the main data pseudocell data? Default is T
#' @param GO Logical. Should GO term enrichment analyses be performed? Default is T
#' @param interactive Logical. Turns off the initial question to continue with the analyses.
#' @return No inline output. It saves an html report, as well as a list object with the resulting WGCNA data. They are both named using the project name and the date. It also creates a folder containing the network files per module.
#' @export
#' @importFrom WGCNA bicor
#' @examples
#' # Calculate pseudocells
#' ps.pbmc_small=calculate.pseudocells(Seurat::pbmc_small, dims = 1:10)
#' 
#' gnames= data.frame(a=as.character(rownames(ps.pbmc_small))); gnames[,2]=gnames[,1]
#' 
#' # Use pseudocells and single cells to calculate WGCNA
#' #scWGNA.report(data = ps.pbmc_small, sc.data = Seurat::pbmc_small, gene.names = gnames, project.name = "test", sp = "Hs")
#' 

scWGNA.report = function(data,sc.data,gene.names, project.name, sp="Mm", cells=F, features=F, reduction="tsne", dir="./", is.pseudocell=T,GO=T,interactive=T) {
  
  if (interactive) {
    my.question <- utils::menu(c("Y", "N"), title = "WARNING: This funciton might take a long time, depending on the data, and will make changes in your computer (save different files).
                             Would you like to continue? (Y/N)")
    if(my.question == 2){
      return(cat(""))
    }
  }
  
  #If running GO analyses
  if (GO==T) {
    #Check that the species package of biomart is installed
     if (!paste0("org.",sp,".eg.db") %in% rownames(utils::installed.packages())) {
      stop(paste0("org.",sp,".eg.db")," is not installed. This is necessary to carry out GO Term analyses. Please install and try again or run the analyses with GO=F")
    }
  }
  
  # Knit the report
  rmarkdown::render(
    system.file("WGCNA_scripts", "scWGCNA_report.Rmd", package = "scWGCNA"), params = list(
      #Parse the parametes
      data = data,
      sc.data = sc.data,
      gene.names = gene.names,
      project.name = project.name,
      sp = sp,
      cells = cells,
      features = features,
      reduction = reduction,
      dir = dir,
      is.pseudocell = is.pseudocell,
      GO=GO
    ),
    output_file = paste0(dir,"WGCNA_report_", project.name, "_", format(Sys.Date(), "%d%m%y"), ".html"),
    output_format = "html_document"
  )
}
