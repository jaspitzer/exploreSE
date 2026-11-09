# exploreSE 0.99.2

* Preparing for initial submission

# exploreSE 0.99.3
* fixed build warning related to the use of shiny related functions

# exploreSE 0.99.4
* fixed an error in the expression plots related to geom_violin() and the text aesthetic
* fixed an error with point coloring based on cutoffs in the volcano plot
* fixed various typos
* added installation instructions to the vignette
* removed unused code

# exploreSE 0.99.5
* fixed R dependency in the description file


# exploreSE 0.99.6
* adjusted some things in the README and added a table below the FCFC plot 
* fixed a bug where the results stored in the metadata as no data frames where not converted

# exploreSE 0.99.7
* removed messages and warnings from the readme and the vignette to increase legibility
* added a a=b  line in the FC FC plot for visibility
* fixed an issue in the test suite caused by filtering not significant hits from the fc-fc plots
* fixed a bug where the gene selection did not update based on gene type
* added an overview table of the different DE comparisons where a given gene is considered DE in the bottom of the gene expression panel 

# exploreSE 0.99.8
* put the data loading interface into a conditional, so that it only appears when no data is loaded

# exploreSE 0.99.9
* various bug fixes
* missed a version bump in both description and news

# exploreSE 0.99.10
* fixed a bug with the conditional data loading 

# exploreSE 0.99.11
* fixed a bug with the PCA when demo data was loaded first. 

# exploreSE 0.99.12
* adjusted the rendering of the DE table for better legibility
* excluded list columns as possible gene identifier columns

# exploreSE 0.99.13
* adjusted alpha for the FC FC lpot to be inline with the rest of the visualisations
* added a toggle for log scale on the expression plots
* added to toggle to hide the text for gene driving enrichments in the enrichment plots

# exploreSE 0.99.14
* added a context tab to overlay specified context over the respective comparison
* bug fixes with regards to demo data and gene identifiers
* minor refactor for parsing of input lists

# exploreSE 0.99.15
* fixed spelling and style in the vignette and README

# exploreSE 0.99.16
* bug fix in get.gos that ran code twice by accident

# exploreSE 0.99.17
* added example contexts
* updates to vignette and readme to reflect the example context

# exploreSE 0.99.18
* edited the example context, updates sources 
* added a landing page with thumbnails and introduction 
* added cross-identifier matching for the context lists
* accordingly, updated to the vignette and readme

# exploreSE 0.99.19
* author update