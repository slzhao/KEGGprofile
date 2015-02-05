KEGGprofile
============

# Introduction #
KEGGprofile combined the KEGG pathway map with expression
profiles of genes in that pathway and facilitated more detailed analysis
about the specific function changes inner pathway or temporal correlations
in different genes and samples. KEGGprofile also supports to display compound and gene expression data at the same time.

# Download and install #
You can install KEGGprofile package in R from [Bioconductor](http://bioconductor.org/packages/release/bioc/html/KEGGprofile.html) by following R codes:

    source("http://bioconductor.org/biocLite.R")
    biocLite("KEGGprofile")

# Document #
KEGGprofile has a vignette in [Bioconductor](http://bioconductor.org/packages/release/bioc/vignettes/KEGGprofile/inst/doc/KEGGprofile.pdf) to demonstrate its application on high-throughput expression data.

# Web interface #
KEGGprofile has a web interface at [CQS website](https://cqs.mc.vanderbilt.edu/shiny/KEGGprofile/). It allows users without R background to perform KEGG database based analysis and visualization easily.