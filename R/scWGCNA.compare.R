#' Runs a comparative scWGCNA analysis
#' 
#' This function runs a WGCNA analysis adapted for single cells. Based on single-cell or pseudocell data.
#' @name scWGCNA.compare.report
#' @usage scWGNA.compare.report(data,test,test.names,project.name,ortho,ortho.sp,groups,dir="./",interactive=T)
#' @param data Seurat object. The WGCNA data to use as reference for comparative analysis. Output from scWGCNA.report.
#' @param test List of seurat objects. The samples to test conservation in. If data is calculated on pseudocells, these should also be pseudocells
#' @param test.names Vector of strings. Contains the names for each of the samples in test.
#' @param project.name String. The name of the project, will be used to name the saved objects and the report.
#' @param ortho Data frame. As many columns as species used only listing 1-to-1 orthologous genes. The species of the reference MUST always be column no. 1
#' @param ortho.sp Vector numeric. Indicating to which species in ortho each test sample belongs to.
#' @param groups List of vectors. Indicating if the plots should be grouped. Otherwise, all samples will be plotted together. Suggested for groups of >5 samples. Default is F
#' @param dir String. The directory where the output will be saved. It defaults to the working directory.
#' @param interactive Logical. Turns off the initial question to continue with the analyses.
#' @return No inline output. It saves an html report, as well as an object with the resulting WGCNA comparison indices. They are both named using the project name and the date.
#' @export
#' @importFrom WGCNA bicor
#' 

scWGNA.compare.report = function(data,test,test.names,project.name,ortho,ortho.sp,groups=F,dir="./",interactive=T) {
  
  if (interactive) {
    my.question <- utils::menu(c("Y", "N"), title = "WARNING: This funciton might take a long time, depending on the data, and will make changes in your computer (save different files).
                             Would you like to continue? (Y/N)")
    if(my.question == 2){
      return(cat(""))
    }
  }
  
  # Knit the report
  rmarkdown::render(
    system.file("WGCNA_scripts", "scWGCNA.comparative_report.Rmd", package = "scWGCNA"), params = list(
      #Parse the parametes
      data = data,
      test = test,
      test.names = test.names,
      project.name = project.name,
      ortho = ortho,
      ortho.sp = ortho.sp,
      groups = groups
    ),
    output_file = paste0(dir,"WGCNA_comparative_report_", project.name, "_", format(Sys.Date(), "%d%m%y"), ".html"),
    output_format = "html_document"
  )
}
