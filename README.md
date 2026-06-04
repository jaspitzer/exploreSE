
<!-- README.md is generated from README.Rmd. Please edit that file -->

# exploreSE

<!-- badges: start -->

<!-- badges: end -->

The goal of exploreSE is to create and maintian a shiny app for
transcriptomic results presentation. It heavily relies on the
[SummarizedExperiment](https://www.bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html)
and
[DeeDeeExperiment](https://bioconductor.org/packages//release/bioc/html/DeeDeeExperiment.html)
packages on bioconductor.

## Installation

You can install the development version of exploreSE from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("jaspitzer/exploreSE")
```

## Example

Let’s say you’ve loaded in your data and created two separate models:

``` r
library(exploreSE)
library(airway)
#> Loading required package: SummarizedExperiment
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: 'MatrixGenerics'
#> The following objects are masked from 'package:matrixStats':
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: 'generics'
#> The following objects are masked from 'package:base':
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
#>     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
#>     get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
#>     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
#>     rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: 'S4Vectors'
#> The following object is masked from 'package:utils':
#> 
#>     findMatches
#> The following objects are masked from 'package:base':
#> 
#>     expand.grid, I, unname
#> Loading required package: IRanges
#> Loading required package: Seqinfo
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: 'Biobase'
#> The following object is masked from 'package:MatrixGenerics':
#> 
#>     rowMedians
#> The following objects are masked from 'package:matrixStats':
#> 
#>     anyMissing, rowMedians
library(DeeDeeExperiment)
#> Loading required package: SingleCellExperiment
library(DESeq2)
data(airway)
```

``` r
airway <- DESeqDataSet(airway, design = ~dex)
airway <- DESeq(airway)
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing
baseline <- results(airway)

airway <- DESeqDataSet(airway, design = ~dex + cell)
airway <- DESeq(airway)
#> using pre-existing size factors
#> estimating dispersions
#> found already estimated dispersions, replacing these
#> gene-wise dispersion estimates
#> found already estimated gene-wise dispersions, removing these
#> mean-dispersion relationship
#> final dispersion estimates
#> found already estimated dispersions, removing these
#> fitting model and testing
cell_controlled <- results(airway)
```

These results can be added to the DeeDeeExperiment

``` r
airway <- addDEA(DeeDeeExperiment(airway), baseline)
airway <- addDEA(airway, cell_controlled)
```

``` r
airway
#> class: DeeDeeExperiment 
#> dim: 63677 8 
#> metadata(2): '' version
#> assays(4): counts mu H cooks
#> rownames(63677): ENSG00000000003 ENSG00000000005 ... ENSG00000273492
#>   ENSG00000273493
#> rowData names(62): gene_id gene_name ... cell_controlled_pvalue
#>   cell_controlled_padj
#> colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
#> colData names(10): SampleName cell ... BioSample sizeFactor
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> dea(2): baseline, cell_controlled 
#> fea(0):
```

``` r
airway <- get.gos(obj = airway, NAME = "baseline", gene_type = "ENSEMBL")
#> 
#> 
#> Found 3424 gene sets in `enrichResult` object, of which 132 are significant.
#> Converting for usage within the DeeDeeExperiment framework...
#> ✔ Renamed FEA entries: "up_go" to "baseline_up_go"
#> Found 3589 gene sets in `enrichResult` object, of which 127 are significant.
#> Converting for usage within the DeeDeeExperiment framework...
#> ✔ Renamed FEA entries: "dn_go" to "baseline_down_go"
airway <- get.gos(obj = airway, NAME = "cell_controlled", gene_type = "ENSEMBL")
#> Found 2248 gene sets in `enrichResult` object, of which 23 are significant.
#> Converting for usage within the DeeDeeExperiment framework...
#> ✔ Renamed FEA entries: "up_go" to "cell_controlled_up_go"
#> Found 3205 gene sets in `enrichResult` object, of which 124 are significant.
#> Converting for usage within the DeeDeeExperiment framework...
#> ✔ Renamed FEA entries: "dn_go" to "cell_controlled_down_go"
```

``` r
airway
#> class: DeeDeeExperiment 
#> dim: 63677 8 
#> metadata(2): '' version
#> assays(4): counts mu H cooks
#> rownames(63677): ENSG00000000003 ENSG00000000005 ... ENSG00000273492
#>   ENSG00000273493
#> rowData names(62): gene_id gene_name ... cell_controlled_pvalue
#>   cell_controlled_padj
#> colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
#> colData names(10): SampleName cell ... BioSample sizeFactor
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> dea(2): baseline, cell_controlled 
#> fea(4): baseline_up_go, baseline_down_go, cell_controlled_up_go, cell_controlled_down_go
```

With the enrichments added, we can now call our shiny app for
exploration! (Note: this is a work in progress, I still need to figure
out how to properly show this in the vignette)

``` r
app <- exploreSE()
if(interactive()){
  shiny::runApp(app, port = 1234)
}
```
