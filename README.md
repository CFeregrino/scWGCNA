# scWGCNA
scWGCNA is an adaptation of WGCNA to work with single-cell datasets.
The functionality is presented in  <a href="https://doi.org/10.1101/2021.02.09.430383" target="_blank">Feregrino & Tschopp 2021</a>

Our package is still in beta testing and not fully functional, but the scripts in the "int" folder can be accessed and used to perform the WGCNA analyses.

scWGCNA works with Seurat objects and produces integrated HTML reports of WGCNA analyses.

## Installation

To install scWGCNA run in R:
```
devtools::install_github("cferegrino/scWGCNA", ref="main")
library(scWGCNA)

```

## Applications
* Calculate pseudocells from a Seurat object with pre.calculated PCA (or other reductions) and cell clusters
```
calculate.pseudocells(seurat, seeds=0.2, nn = 10, reduction = "pca")
```
* Perform WCGNA analyses and output an integrated HTML report
```
scWGNA.report(data, sc.data, gene.names, project.name, sp="Mm")
```
* Perform comparative WCGNA analyses and output an integrated HTML report
```
scWGNA.compare.report(data, test, test.names, project.name, ortho, ortho.sp)
```

## Output examples
These are the HTML outputs you can expect from the functions.
The data used in our publication produced [this HTML report output](https://htmlpreview.github.io/?https://github.com/CFeregrino/scWGCNA/blob/main/HTMLexamples/WGCNA_report_E15test3_080421.html) from the scWGNA.report function, and [this HTML report output](https://htmlpreview.github.io/?https://github.com/CFeregrino/scWGCNA/blob/main/HTMLexamples/WGNA_comparative_E15.nb.html) from the scWGNA.compare.report function.

## References

<a href="https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559" target="_blank">Langfelder, P., Horvath, S. (2008) WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 9, 559</a>

<a href="https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1001057" target="_blank">Langfelder P, Luo R, Oldham MC, Horvath S (2011) Is My Network Module Preserved and Reproducible?. PLOS Computational Biology 7(1): e1001057</a>

<a href="https://doi.org/10.1101/2021.02.09.430383" target="_blank">Feregrino, C, Tschopp, P (2021) Assessing evolutionary and developmental transcriptome dynamics in homologous cell types. bioRxiv 2021.02.09.430383</a>
