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



Lab 1: feature selection & dimension reduction
=====


## Preparation

**Load data**.

You can download the R data file tse.Rds from
[here](https://github.com/microbiome/course_2021_ml4microbiome/blob/master/tse.Rds). The
read it to R with:

```{r}
tse <- readRDS("tse.Rds")
```

If you like to see how the data was created, you can download the
[zip](data_zipfile.zip) and run the data preparation
[script](tse_script.R).


## Load necessary packages in R

```{r}
library(mia)
library(miaViz)
```


**Familiarize** with this R data format based on examples in [online
tutorial](https://microbiome.github.io/course_2021_miaverse/microbiome-data-exploration.html). We
will go through some basics after you have successfully loaded the
data.

Use the following resources to solve the tasks:

- Online book [Orchestrating Microbiome Analysis](https://microbiome.github.io/OMA/)
- [mia R package function help pages](https://microbiome.github.io/mia/reference/index.html)
- [miaViz R package function help pages](https://microbiome.github.io/miaViz/reference/index.html)


## Tasks (hints of useful functions in parentheses)

Feature selection will determine which aspects of the data we will focus on.

1.  **Aggregate** data to your desired taxonomic Phylum level ([mia::agglomerateByRank](https://microbiome.github.io/mia/reference/agglomerate-methods.html))

2.  **Subset** the data at your desired taxonomic level with only
    taking taxonomic groups that are greater than 10% of prevalence in
    the total sample. ([mia::subsetByPrevalentTaxa](https://microbiome.github.io/mia/))
    

3.  **Subset** the data to selected two phyla, for instance:
    "Actinobacteria“ and ”Cyanobacteria"

4.  **Transform** the data with relative abundances and centered
    log-ratio (clr) (i.e. [mia::transformSamples](https://microbiome.github.io/mia/reference/transformCounts.html))



Lab 2: visualization & clustering
=====

Visualization and clustering help to discover patterns in the data in a reproducible way.

## Tasks

1. Visualize beta diversity using principal coordinate analysis
   (PCoA); based on the Bray-Curtis dissimilarities

2. Test how results are affected by feature selection and
   dissimilarity measure (for instance, Aitchison distance i.e. clr
   transformation + Euclidean distance)

3. Cluster the samples using Dirichlet-multinomial mixture model (DMM)

4.  Visualize the clusters in the PCoA plot


Lab 3: analysis and interpretation
=====


Finally, the identified patterns require interpretation.

## Tasks

1.  What taxa are driving the axis? Calculate the Spearman correlation
    between PC1 and the relative abundance of the bacteria and visualize
    the results in a bar plot.

2.  Visualize the dominant taxonomic group of each sample by colour on
    PCoA plot

3.  Visualize gender by colour on PCoA plot

4.  Experiment with the other data exploration techniques in [OMA](https://microbiome.github.io/OMA/)