# deconvo_method-paper

Contains analysis for the deconvolution methods paper

## Details about this repository

The path to this repository on JHPCE is:

This repo contains analysis scripts, outputs, figures, etc. for the deconvolution methods paper. The subdirs contained here include:

* `scripts` : Main location where scripts and scripts to generate content live.

* `figures_and_tables`: Main location of figures and tables, with same layout as `scripts` subdir

* `outputs`: Objects (i.e. not figures or tables) which are produced by scripts from `scripts`. Has same layout as `scripts` subdir.

* `source`: Scripts, vignettes, etc. which are central to the planned deconvolution resource. This will have its own internal structure and be developed in parallel to analyses for the paper.

## Data files

`SingleCellExperiment` objects containing single nucleus RNA-seq data.

* `data/sce`  

** `data/sce/SCE_DLPFC...` snRNAseq Dorsolateral prefrontal cortex data from Tran et al 2021. Downloaded from ![link](https://github.com/LieberInstitute/10xPilot_snRNAseq-human)

** `data/sce/sce...` snRNAseq Dorsolateral prefrontal cortex data, new to this study.

## Manuscript figures

[Figure 1](./figures/lute_diagram.jpg) : Diagram of `lute` deconvolution framework.

* Figure 2
[Figure 2A](./notebooks/02_pseudobulk/) : Scatterplot of true versus predicted proportions from pseudobulk experiment, showing results in neuron from k2 labels tests, using DLPFC Cohort1.
[Figure 2B](./notebooks/02_pseudobulk/) : Scatterplot of true versus predicted proportions from pseudobulk experiment, showing results in neuron from k2 labels tests, using DLPFC Cohort2.
[Figure 2C](./notebooks/02_pseudobulk/) : Scatterplot of true versus predicted proportions from pseudobulk experiment, showing results in Excit, Inhib, and glial from k3 label tests, using DLPFC Cohort1.
[Figure 2D](./notebooks/02_pseudobulk/) : Scatterplot of true versus predicted proportions from pseudobulk experiment, showing results in Excit, Inhib, and glial from k3 label tests, DLPFC Cohort2.
[Figure 2E](./notebooks/02_pseudobulk/) : Scatterplot of true versus predicted proportions from pseudobulk experiment, showing results in Excit, Inhib, and glial from k4 label tests, using DLPFC Cohort1.
[Figure 2F](./notebooks/02_pseudobulk/) : Scatterplot of true versus predicted proportions from pseudobulk experiment, showing results in Excit, Inhib, and glial from k4 label tests, DLPFC Cohort2.

* Figure 3
[Figure 3A](./figures/03_shuffle/) Shuffle $S_{pseudobulk}$
[Figure 3B-C](./notebooks/03_shuffle/) Shuffle $S_{pseudobulk}$
[Figure 3D-E](./notebooks/03_shuffle/) Shuffle $S_{deconvolution}$
[Figure 3F-G](./notebooks/03_shuffle/) Shuffle $Z_{pseudobulk}$

* Figure 4
[Figure 4A](./scripts/04_experiment) 
[Figure 4B](./scripts/04_experiment) 
[Figure 4C](./scripts/04_experiment)

* Supplemental Figure S1
[Figure S1A](./scripts/02_pseudobulk/) : Paneled scatterplots of adjusted and unadjusted results of true versus predicted proportions from pseudobulk experiment, DLPFC cohort 1.
[Figure S1B](./scripts/02_pseudobulk/) : Jittered points overlaying cyan quantile boxplots of absolute errors for adjusted and unadjusted results, DLPFC cohort 1.
[Figure S1C](./scripts/02_pseudobulk/) : Paneled scatterplots of adjusted and unadjusted results of true versus predicted proportions from pseudobulk experiment, DLPFC cohort 2.
[Figure S1D](./scripts/02_pseudobulk/) : Jittered points overlaying cyan quantile boxplots of absolute errors for adjusted and unadjusted results, DLPFC cohort 1 k2 tests.
[Figure S1E](./scripts/02_pseudobulk/) : Jittered points overlaying cyan quantile boxplots of absolute errors for adjusted and unadjusted results, DLPFC cohort 1, k3 tests.
[Figure S1F](./scripts/02_pseudobulk/) : Jittered points overlaying cyan quantile boxplots of absolute errors for adjusted and unadjusted results, DLPFC cohort 1, k3 tests.
[Figure S1E](./scripts/02_pseudobulk/) : Jittered points overlaying cyan quantile boxplots of absolute errors for adjusted and unadjusted results, DLPFC cohort 1, k4 tests.
[Figure S1F](./scripts/02_pseudobulk/) : Jittered points overlaying cyan quantile boxplots of absolute errors for adjusted and unadjusted results, DLPFC cohort 2, k4 tests.

* Supplemental Figure 2
[Supplemental Figure 2A](./scripts/06_bulk/) : 
[Supplemental Figure 2B](./scripts/06_bulk/) : 
[Supplemental Figure 2C](./scripts/06_bulk/) : 

* Supplemental Figure 3


## File locations
This is an integrative, multi-assay project including individual-matched data generated from human dorsolateral prefrontal cortex (DLPFC) using the RNAscope, 10X Chromium, and Illumina HiSeq platforms. 

The main directory path on JHPCE is located at:

`/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/`

Note that integrative analyses and content for the deconvolution methods paper are stored in the subdirectory `Human_DLPFC_Deconvolution`.

### 10X Chromium single-nucleus RNA-seq datasets

The snRNA-seq datasets for this project currently include 9 donors and 10 samples. For the deconvolution method paper, we access these data as a `SingleCellExperiment` object. The latest (final?) version is a file called `sce_DLPFC.Rdata` located at the subdirectory:

`DLPFC_snRNAseq/processed-data/sce/`

### Bulk RNA-seq datasets

The bulk RNA-seq datasets include 6 samples groups comprised of mRNA isolated from either nucleus, cytoplasm, or bulk, and prepped with either polyA selection or rRNA depletion strategies. Data come from 19 samples, meaning 113 total data points. 

These bulk RNA-seq data are currently being QC'd, but they will be accessible as a `SummarizedExperiment` object located at the subdirectory:

`Human_DLPFC_Deconvolution/processed-data/01_SPEAQeasy/`

### RNAScope/HALO image datasets

Image data is generated for RNAscope slides by analysis with HALO and outputting the analysis results as `.csv` tables. These tables are read from the file tree located at the subdirectory:

`Human_DLPFC_Deconvolution/raw-data/HALO/Deconvolution_HALO_analysis/`

Note, this contains one subdirectory for each of the RNAscope experiments, called `Circle` and `Star`, respectively, These experiments are distinct in that they each comprise of an analysis of an independent, albeit adjacent, tissue and each includes a different set of molecular markers (see table below).

## Variables and data dictionaries

The deconvolution project makes use of a number of metadata attributes and variables in the results files mentioned about. This section describes the key terms and definitions of these attributes and variables for the deconvolution method paper.

### Cell type labels

Cell type labels for `snRNA-seq` datasets are determined from the variable `cellType_broad_hc`. This can be accessed from the `sce` object in various ways such as `colData(sce)$cellType_broad_hc` or `sce[["cellType_broad_hc"]]`. Note that the deconvolution methods paper focuses on just 6 cell types of interest, and these are identified from among the cell type labels in the `cellType_broad_hc` variable.

Cell type labels for RNAscope experiments are obtained from the image analysis outputs produced by the HALO software. In brief, outputs each contain a series of columns corresponding to the cell type makers. Since each row in these outputs corresponds to an individual detected nucleus, we simply look at which marker is positive for that nucleus to determine its cell type. Cell type proportions and abundances are then calculated from these outputs.

Cell type labels aren't available for the bulk RNA-seq and other datasets produced for this project.

### Data dictionaries for marker tables, results, etc.

### RNAscope labels and marker types

The molecular markers for the `Circle` and `Star` RNAscope experiments are as follows:

```
# rnascope markers/expts
0 cellType marker  Combo	Type		LongName
1 Endo     CLDN5   Circle	Ab 			Claudin_5 
2 Astro    GFAP    Circle	Ab 			GFAP
3 Inhib    GAD1    Circle	RNA_probe	GAD1
4 Excit    SLC17A7 Star		RNA_probe	SLC17A7
5 Micro    TMEM119 Star		Ab 			TMEM119
6 Oligo    OLIG2   Star		Ab			OLIG2
```

## Terminology

Throughout the deconvolution project we try to use standard terms to refer to key project entities. Here are some of the
key terms to be aware of when using these project files and understanding analysis outputs.

### The deconvolution equation

We take deconvolution to be the prediction of cell type amounts in a mixed sample by leveraging data from a non-mixed sample. For the deconvolution method paper, we focus on predicting either cell abundances or proportions for each of 6 cell types in bulk RNA-seq data by leveraging a single-nucleus RNA-seq reference dataset. The deconvolution equation looks like:

$Y = \pi * Z$

Where $Y$ is a matrix of mixed signals having dimensions $[G,J]$, $\pi$ is a vector of cell type amounts (proportions or abundances), and $Z$ is a matrix of non-mixed signals having dimensions $[G,K]$. The dimensions of these objects include the $G$ marker gene set, $J$ samples with mixed signal, and $K$ cell types.

### Strict deconvolution

Strict deconvolution refers to solving for $\pi$ given a set of matrices $[Y,Z]$. In other words, if we assume we have a set of matrices $Y$ and $Z$, we then solve for $\pi$ based on these.

### Deconvolution preparation

This simply refers to steps taken to prepare $Y$ and $Z$ prior to performing deconvolution. This may include differential weighting of the marker genes, bias corrections, data transformations and rescaling, etc.

### Marker genes

The $G$ marker genes represent markers of cell types which we can use in deconvolution. Note, this is the "row" dimension in the $Y$ and $Z$ matrices.

### Samples

The $J$ samples are the samples having mixed signal, for which we predict cell type amounts. Note, this is the "column" dimension in the matrix $Y$.

### Cell types

The $K$ cell types represent the cell types for which we deconvolute signal and obtain amount estimates. Note, this is the "column" dimension in matrix $Y$ and the length of the vector represented by $\pi$.

## Funding sources

* 1R01MH123183-01 ([grantome_link](https://grantome.com/grant/NIH/R01-MH123183-01))