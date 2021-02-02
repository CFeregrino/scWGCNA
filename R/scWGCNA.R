#' Runs a semi-automatic, iterative scWGCNA analysis
#' 
#' This function runs our semi-automatic single-cell WGCNA analysis. It runs in an iterative way. Based on single-cell or pseudocell data.
#' @name scWGCNA.report
#' @usage scWGNA.report(data,sc.data,gene.names,project.name,sp=Mm,cells=F,features=F,reduction="tsne",dir="./",is.pseudocell=T)
#' @param data The expression data to run the analysis in. In a Seurat object. Can be pseudocell or single-cell data. We recommend pseudocells
#' @param sc.data The single cell data, if running on single cell data already, please repear the argument
#' @param gene.names Translation of the identifiers found in the expression data. Two columns: 1= ids present in expression matrix, 2= names to appear in plots. Can be the same
#' @param project.name The name of the project, will be used to save objects and the report
#' @param sp The species used for the analysis in string form. Only two letters in the format Sp. Defaults to "Mm"
#' @param cells Are certain clusters to be used? Please use identities and not cell names
#' @param features The features to be used for the analysis. Default is F, which makes the script calculate variable genes
#' @param reduction The reduction to use for plots. It must be present in the @@reductions slot of sc.data. Default is "tsne"
#' @param dir The directory to save the output to. It defaults to the working directory
#' @param is.pseudocell is the main data pseudocell data? Default is T
#' @return Saves a report as well as different objects.
#' @export
#' 

scWGNA.report = function(data,sc.data,gene.names, project.name, sp="Mm", cells=F, features=F, reduction="tsne", dir="./", is.pseudocell=T) {
  
  if (!paste0("org.",sp,".eg.db") %in% rownames(utils::installed.packages())) {
         stop(paste0("org.",sp,".eg.db")," is not installed. This is necessary to carry out GO Term analyses. Please install and try again")
     }
  
  rmarkdown::render(
    system.file("WGCNA_scripts", "scWGCNA_report", package = "scWGCNA"), params = list(
      data = data,
      sc.data = sc.data,
      gene.names = gene.names,
      project.name = project.name,
      sp = sp,
      cells = cells,
      features = features,
      reduction = reduction,
      dir = dir,
      is.pseudocell = is.pseudocell
    ),
    output_file = paste0(dir,"WGCNA_report_", project.name, "_", format(Sys.Date(), "%d.%m.%y"), ".pdf"),
    output_format = "pdf_document"
  )
}
