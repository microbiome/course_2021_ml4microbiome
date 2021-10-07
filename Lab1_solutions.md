\#command to load that Rds file ( )

Tutorial Aims
=============

In this tutorial we will learn feature selection, dimension reduction,
clustering, visualization, analysis and interpretation using  
miaverse (mia = MIcrobiome Analysis) and miaViz - (Microbiome Analysis
Plotting and Visualization) to explore patterns in human gut microbiome
datasets.

The example in this tutorial are mainly based on this online book  
[Orchestrating Microbiome Analysis](https://microbiome.github.io/OMA/)

In this tutorial we use data from [The Human Gut Microbiome Atlas
(HGMA)](https://www.microbiomeatlas.org/)

Lab 1: feature selection & dimension reduction
==============================================

Load the data
-------------

Let us load the readily processed data (TreeSummarizedExperiment
object).

You can download the R data file tse.Rds from
[here](https://github.com/microbiome/course_2021_ml4microbiome/blob/master/tse.Rds).

    tse <- readRDS("tse.Rds")

If you like to see how the data was created, you can download the
[zip](data_zipfile.zip) and run the data preparation
[script](tse_script.R).

Load packages in R
------------------

    library(mia)

    ## Loading required package: SummarizedExperiment

    ## Loading required package: MatrixGenerics

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'MatrixGenerics'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars

    ## Loading required package: GenomicRanges

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## Loading required package: GenomeInfoDb

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attaching package: 'Biobase'

    ## The following object is masked from 'package:MatrixGenerics':
    ## 
    ##     rowMedians

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     anyMissing, rowMedians

    ## Loading required package: SingleCellExperiment

    ## Loading required package: TreeSummarizedExperiment

    ## Loading required package: Biostrings

    ## Loading required package: XVector

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

    library(miaViz)

    ## Loading required package: ggplot2

    ## Loading required package: ggraph

Tasks
-----

1.  Aggregate the data to Phylum level

<!-- -->

    tse_phylum <- agglomerateByRank(tse, rank = "phylum")

1.  Subset the data at the Species level with only taking bacterial
    Species that are greater than 10% of prevalence in the total sample

<!-- -->

    altExp(tse,"species") <- agglomerateByRank(tse, "species")
    getPrevalentTaxa(tse, detection = 0, prevalence = 10/100)

    ##   [1] "msp_0003"  "msp_0005"  "msp_0007"  "msp_0008"  "msp_0009"  "msp_0010" 
    ##   [7] "msp_0011"  "msp_0012"  "msp_0013"  "msp_0014"  "msp_0015"  "msp_0016" 
    ##  [13] "msp_0017"  "msp_0018"  "msp_0019"  "msp_0020"  "msp_0021"  "msp_0022" 
    ##  [19] "msp_0023"  "msp_0024"  "msp_0025"  "msp_0026"  "msp_0027"  "msp_0028" 
    ##  [25] "msp_0029"  "msp_0030"  "msp_0031"  "msp_0032"  "msp_0033"  "msp_0034" 
    ##  [31] "msp_0035"  "msp_0036"  "msp_0037"  "msp_0038"  "msp_0039"  "msp_0040" 
    ##  [37] "msp_0041"  "msp_0042"  "msp_0043"  "msp_0044"  "msp_0045"  "msp_0046" 
    ##  [43] "msp_0047"  "msp_0048"  "msp_0050"  "msp_0052"  "msp_0053"  "msp_0054" 
    ##  [49] "msp_0055"  "msp_0056"  "msp_0057"  "msp_0058"  "msp_0059"  "msp_0060" 
    ##  [55] "msp_0062"  "msp_0063"  "msp_0065"  "msp_0066"  "msp_0068"  "msp_0069" 
    ##  [61] "msp_0071"  "msp_0072"  "msp_0073"  "msp_0074"  "msp_0075"  "msp_0076" 
    ##  [67] "msp_0077"  "msp_0079"  "msp_0081"  "msp_0083"  "msp_0086"  "msp_0087" 
    ##  [73] "msp_0089"  "msp_0090"  "msp_0092"  "msp_0093"  "msp_0095"  "msp_0103" 
    ##  [79] "msp_0105"  "msp_0106"  "msp_0107"  "msp_0110"  "msp_0111"  "msp_0117" 
    ##  [85] "msp_0121"  "msp_0124"  "msp_0125"  "msp_0126"  "msp_0128"  "msp_0129" 
    ##  [91] "msp_0130"  "msp_0131"  "msp_0132"  "msp_0133"  "msp_0134"  "msp_0138" 
    ##  [97] "msp_0139"  "msp_0141"  "msp_0142"  "msp_0144"  "msp_0145"  "msp_0148c"
    ## [103] "msp_0151"  "msp_0152"  "msp_0153"  "msp_0154"  "msp_0158"  "msp_0159" 
    ## [109] "msp_0160"  "msp_0161"  "msp_0162"  "msp_0163"  "msp_0164"  "msp_0166" 
    ## [115] "msp_0168"  "msp_0169"  "msp_0170"  "msp_0172"  "msp_0173"  "msp_0175" 
    ## [121] "msp_0178"  "msp_0183"  "msp_0186"  "msp_0188"  "msp_0189"  "msp_0192" 
    ## [127] "msp_0194"  "msp_0195"  "msp_0196"  "msp_0198"  "msp_0199"  "msp_0202" 
    ## [133] "msp_0203"  "msp_0205"  "msp_0206"  "msp_0207"  "msp_0208"  "msp_0212" 
    ## [139] "msp_0213"  "msp_0215"  "msp_0220"  "msp_0224"  "msp_0225"  "msp_0226" 
    ## [145] "msp_0227"  "msp_0230"  "msp_0232"  "msp_0236"  "msp_0237"  "msp_0238" 
    ## [151] "msp_0239"  "msp_0243"  "msp_0244"  "msp_0246"  "msp_0249"  "msp_0250" 
    ## [157] "msp_0258"  "msp_0259"  "msp_0262"  "msp_0263"  "msp_0265"  "msp_0269" 
    ## [163] "msp_0270"  "msp_0271"  "msp_0274"  "msp_0277"  "msp_0280"  "msp_0281" 
    ## [169] "msp_0285"  "msp_0287"  "msp_0288"  "msp_0290"  "msp_0291"  "msp_0296" 
    ## [175] "msp_0298"  "msp_0299"  "msp_0301"  "msp_0302"  "msp_0303"  "msp_0305" 
    ## [181] "msp_0306"  "msp_0307"  "msp_0308"  "msp_0312"  "msp_0313"  "msp_0314" 
    ## [187] "msp_0315"  "msp_0317"  "msp_0318"  "msp_0319"  "msp_0324"  "msp_0325" 
    ## [193] "msp_0331"  "msp_0334"  "msp_0335"  "msp_0337"  "msp_0340"  "msp_0342" 
    ## [199] "msp_0345"  "msp_0346"  "msp_0348"  "msp_0352"  "msp_0353"  "msp_0356" 
    ## [205] "msp_0357"  "msp_0360"  "msp_0361"  "msp_0362"  "msp_0364"  "msp_0369" 
    ## [211] "msp_0373"  "msp_0374"  "msp_0375"  "msp_0378"  "msp_0380"  "msp_0381" 
    ## [217] "msp_0388"  "msp_0389"  "msp_0392"  "msp_0396"  "msp_0398"  "msp_0399" 
    ## [223] "msp_0404"  "msp_0407"  "msp_0412"  "msp_0418"  "msp_0419"  "msp_0422" 
    ## [229] "msp_0424"  "msp_0425"  "msp_0430"  "msp_0431"  "msp_0432"  "msp_0436" 
    ## [235] "msp_0440"  "msp_0442"  "msp_0447"  "msp_0449"  "msp_0452"  "msp_0454" 
    ## [241] "msp_0456"  "msp_0457"  "msp_0461"  "msp_0462"  "msp_0463"  "msp_0464" 
    ## [247] "msp_0467"  "msp_0468"  "msp_0471"  "msp_0473c" "msp_0477"  "msp_0478" 
    ## [253] "msp_0480"  "msp_0484"  "msp_0489"  "msp_0490"  "msp_0500"  "msp_0506" 
    ## [259] "msp_0507"  "msp_0509"  "msp_0510"  "msp_0512"  "msp_0521"  "msp_0522" 
    ## [265] "msp_0526"  "msp_0530"  "msp_0534"  "msp_0538"  "msp_0542"  "msp_0546" 
    ## [271] "msp_0551"  "msp_0553"  "msp_0554"  "msp_0555"  "msp_0558"  "msp_0561" 
    ## [277] "msp_0563"  "msp_0565"  "msp_0570"  "msp_0572"  "msp_0573"  "msp_0576" 
    ## [283] "msp_0578"  "msp_0581"  "msp_0584"  "msp_0585"  "msp_0586"  "msp_0591" 
    ## [289] "msp_0592"  "msp_0594"  "msp_0595"  "msp_0596"  "msp_0597"  "msp_0604" 
    ## [295] "msp_0613"  "msp_0615"  "msp_0619"  "msp_0621"  "msp_0622"  "msp_0623" 
    ## [301] "msp_0624"  "msp_0630"  "msp_0639"  "msp_0643"  "msp_0648"  "msp_0650" 
    ## [307] "msp_0652"  "msp_0654"  "msp_0655"  "msp_0665"  "msp_0676"  "msp_0679" 
    ## [313] "msp_0685"  "msp_0692"  "msp_0693"  "msp_0697"  "msp_0699"  "msp_0705" 
    ## [319] "msp_0707"  "msp_0709"  "msp_0713"  "msp_0716"  "msp_0722"  "msp_0723" 
    ## [325] "msp_0733"  "msp_0735"  "msp_0740"  "msp_0742"  "msp_0747"  "msp_0756" 
    ## [331] "msp_0757"  "msp_0760"  "msp_0761"  "msp_0763"  "msp_0766"  "msp_0769" 
    ## [337] "msp_0773"  "msp_0776"  "msp_0777"  "msp_0780"  "msp_0781"  "msp_0794" 
    ## [343] "msp_0800"  "msp_0804"  "msp_0808"  "msp_0811"  "msp_0818"  "msp_0820" 
    ## [349] "msp_0821"  "msp_0833"  "msp_0841"  "msp_0850"  "msp_0854"  "msp_0858" 
    ## [355] "msp_0860"  "msp_0861"  "msp_0874"  "msp_0879"  "msp_0881"  "msp_0883" 
    ## [361] "msp_0884"  "msp_0885"  "msp_0893"  "msp_0898"  "msp_0903"  "msp_0906" 
    ## [367] "msp_0908"  "msp_0917"  "msp_0930"  "msp_0931"  "msp_0932"  "msp_0942" 
    ## [373] "msp_0945"  "msp_0947"  "msp_0953"  "msp_0959"  "msp_0970"  "msp_0976" 
    ## [379] "msp_0977"  "msp_0979"  "msp_0986"  "msp_1009"  "msp_1012"  "msp_1021" 
    ## [385] "msp_1032"  "msp_1033"  "msp_1038"  "msp_1062"  "msp_1093"  "msp_1114" 
    ## [391] "msp_1143"  "msp_1173c" "msp_1177"  "msp_1219"  "msp_1236"  "msp_1244" 
    ## [397] "msp_1270"  "msp_1279"  "msp_1291"  "msp_1302"  "msp_1308"  "msp_1323" 
    ## [403] "msp_1339"  "msp_1342"  "msp_1347c" "msp_1349"  "msp_1356"  "msp_1381" 
    ## [409] "msp_1403"  "msp_1461"  "msp_1492"  "msp_1533"  "msp_1556"  "msp_1564" 
    ## [415] "msp_1622"  "msp_1641"  "msp_1643c" "msp_1724"

    prev <- getPrevalentTaxa(tse, detection = 0, prevalence = 10/100)

1.  Subset a tse object which consist only two phyla i.e
    “Actinobacteria”, “Cyanobacteria”

<!-- -->

    tse_subset_by_feature <- tse[rowData(tse)$phylum %in% c("Actinobacteria", "Cyanobacteria")]
    tse_subset_by_feature

    ## class: TreeSummarizedExperiment 
    ## dim: 28 5166 
    ## metadata(0):
    ## assays(1): counts
    ## rownames(28): msp_0130 msp_0328 ... msp_2003 msp_2013
    ## rowData names(6): phylum class ... genus species
    ## colnames(5166): SRR585765 SRR585697 ... SRR5580397 SRR5580399
    ## colData names(11): dataset.ID metadata.ID ... GeneRichness enteroType
    ## reducedDimNames(0):
    ## mainExpName: NULL
    ## altExpNames(1): species
    ## rowLinks: a LinkDataFrame (28 rows)
    ## rowTree: 1 phylo tree(s) (1880 leaves)
    ## colLinks: NULL
    ## colTree: NULL

1.  Calculate relative abundances, and store the table to assays

<!-- -->

    tse <- relAbundanceCounts(tse)

5.Perform a centered log-ratio transformation (clr)
(i.e. mia::transformSamples)

    tse <- transformSamples(x = tse, abund_values = "counts", method = "clr", 
                            pseudocount = 1, name = "CLR")
