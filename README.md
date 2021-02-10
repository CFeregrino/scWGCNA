# scWGCNA
scWGCNA is an adaptation of WGCNA to work with single-cell datasets.
The functionality is presented in  <a href="https://doi.org/10.1101/2021.02.09.430383" target="_blank">Feregrino & Tschopp 2021</a>

Our package is still in beta testing and not fully functional, but the scripts in the "int" folder can be accessed and used to perform the WGCNA analyses.

scWGCNA works with Seurat objects and produces integrated html reports of WGCNA analyses.

## Instalation

To instal scWGCNA run
```
devtools::install_github("cferegrino/scWGCNA", ref="main")
```

## Applications
* Calculate pseudocells from a Seurat object with pre.calculated PCA (or other reductions) and cell clusters
```
calculate.psedocells(seurat, seeds=0.2, nn = 10, reduction = "pca")
```
* Perform WCGNA analyses and output an integrated HTML report
```
scWGNA.report(data,sc.data,gene.names,project.name,sp="Mm")
```
* Perform comparative WCGNA analyses and output an integrated HTML report
```
scWGNA.compare.report(data,test,test.names,project.name,ortho,ortho.sp)
```
## References

<a href="https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559" target="_blank">Langfelder, P., Horvath, S. (2008) WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 9, 559</a>

<a href="https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1001057" target="_blank">Langfelder P, Luo R, Oldham MC, Horvath S (2011) Is My Network Module Preserved and Reproducible?. PLOS Computational Biology 7(1): e1001057</a>

<a href="https://doi.org/10.1101/2021.02.09.430383" target="_blank">Feregrino, C, Tschopp, P (2021) Assessing evolutionary and developmental transcriptome dynamics in homologous cell types. bioRxiv 2021.02.09.430383</a>
