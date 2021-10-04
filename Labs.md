Unsupervised learning: tutorial aims
=============

This tutorial walks us through some key concepts in unsupervised
learning using a set of tools from R/Bioconductor, called miaverse.

The assignments explore a recently released human gut microbiome data
collection from [The Human Gut Microbiome Atlas
(HGMA)](https://www.microbiomeatlas.org/).

We cover aspects of feature selection, dimension reduction,
clustering, visualization, analysis and interpretation, following the
online book [Orchestrating Microbiome
Analysis](https://microbiome.github.io/OMA/) (beta version).



Lab 1
=====

## Feature selection & dimension reduction

### Preparation

Prepare the data. We have provided a script that you can use to load
the data in a TreeSummarizedExperiment (tse) format in R for further
analysis.

You can familiarize with this R data format based on examples in
[online
tutorial](https://microbiome.github.io/course_2021_miaverse/microbiome-data-exploration.html). We
will go through some basics after you have successfully loaded the
data.


### Tasks


1.  Aggregate data to your desired taxonomic Phylum level (mia::agglomerateByRank)

2.  Subset the data at the Species level with only taking bacterial
    Species that are greater than 10% of prevalence in the total sample

3.  Subset a tse object which consist only two phyla i.e
    Actinobacteria“,”Cyanobacteria

4.  Calculate relative abundances, and store the table to assays

6.  Perform a centred log-ratio transformation (clr)
    (i.e. mia::transformSamples)



Lab 2
=====

unsupervised learning: clustering & visualization
-------------------------------------------------

1.  Visualize beta diversity using principal coordinate analysis
    (PCoA);based on the Bray-Curtis dissimilarities

2.  Visualize beta diversity using principal coordinates analysis
    (PCoA); with Aitchison distance (clr transformation+ Euclidean
    distance)

3.  Cluster the samples using Dirichlet- multinomial mixture model

4.  Visualize the clusters in the PCoA plot

Lab 3
=====

unsupervised learning:Analysis and interpretation
-------------------------------------------------

1.  What taxa are driving the axis? Calculate the Spearman correlation
    between PC1 and the relative abundance of the bacteria and visualize
    the results in a bar plot.

2.  Visualize the dominant taxonomic group of each sample by colour on
    PCoA plot

3.  Visualize gender by colour on PCoA plot
