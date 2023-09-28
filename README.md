# deconvo_method-paper

Supplementary materials to prepare data and run analyses for the paper "".

## Paper citation

## Repository structure

Describes the contents of this repo.

### Folder tree

File tree structure returned from `tree deconvo_method-paper -d`:

```
deconvo_method-paper
├── cohort1
│   ├── data
│   ├── env
│   │   ├── 02_pseudobulk
│   │   ├── 03_shuffle
│   │   ├── 04_experiment
│   │   ├── 05_bulk
│   │   ├── 06_estimate
│   │   ├── 07_summary
│   │   ├── 08_adjustment
│   │   └── 09_fast
│   ├── figures
│   │   ├── 02_pseudobulk
│   │   ├── 03_shuffle
│   │   ├── 04_experiment
│   │   ├── 05_bulk
│   │   ├── 06_estimate
│   │   ├── 08_adjustment
│   │   └── 09_fast
│   ├── notebooks
│   │   ├── 00_preprocess
│   │   ├── 01_mae
│   │   ├── 02_pseudobulk
│   │   ├── 03_shuffle
│   │   ├── 04_experiment
│   │   ├── 05_bulk
│   │   ├── 06_estimate
│   │   ├── 07_summary
│   │   ├── 08_adjustment
│   │   ├── 09_fast
│   │   └── bioc2023
│   ├── outputs
│   │   ├── 00_preprocess
│   │   ├── 01_mae
│   │   ├── 02_pseudobulk
│   │   ├── 03_shuffle
│   │   ├── 04_experiment
│   │   ├── 05_bulk
│   │   ├── 06_estimate
│   │   ├── 07_summary
│   │   ├── 08_adjustment
│   │   └── 09_fast
│   ├── scripts
│   │   ├── 00_preprocess
│   │   ├── 01_mae
│   │   ├── 02_pseudobulk
│   │   │   └── old
│   │   ├── 03_shuffle
│   │   ├── 04_experiment
│   │   ├── 05_bulk
│   │   │   └── abtest
│   │   ├── 06_estimate
│   │   ├── 07_summary
│   │   ├── 08_adjustment
│   │   └── 09_fast
│   └── source
├── cohort2
│   ├── env
│   ├── figures
│   ├── outputs
│   └── scripts
├── methods
│   ├── N
│   ├── New folder
│   └── New folder (2)
└── sosina_et_al_2021
```

### Folder guide

* `cohort1` : Supplementary materials for analyses of Cohort1.

* `cohort2` : Supplementary materials for analyses of Cohort2.

* `sosina_et_al_2021` : Supplementary materials for analyses of Sosina et al 2021.

* `cohort1/data` : External data used by supplement.

* `cohort1/env` : Environment (`.RData`) saves from `save.image()`.

* `cohort1/figures` : Figure saves from `cohort1/notebooks`.

* `cohort1/notebooks` : Notebooks (`.Rmd`) which load saved images and save figures.

* `cohort1/outputs` : Objects (`.rda`) saved by scripts.

* `cohort1/scripts` : Scripts (`.Rmd`) performing data preparation and analyses.

* `cohort1/source` : Scripts (`.R`) shared across analysis sections.

This repo contains analysis scripts, outputs, figures, etc. for the deconvolution methods paper. The subdirs contained here include:

* `scripts` : Main location where scripts and scripts to generate content live.

* `figures_and_tables`: Main location of figures and tables, with same layout as `scripts` subdir

* `outputs`: Objects (i.e. not figures or tables) which are produced by scripts from `scripts`. Has same layout as `scripts` subdir.

* `source`: Scripts, vignettes, etc. which are central to the planned deconvolution resource. This will have its own internal structure and be developed in parallel to analyses for the paper.

## Figure links

Contains the links to files saving figures for the manuscript.

### Main figures

**Figure 1**

[Figure 1](./figures/lute_diagram.jpg)

**Figure 2**

[Figure 2A](./notebooks/02_pseudobulk/)

[Figure 2B](./notebooks/02_pseudobulk/)

[Figure 2C](./notebooks/02_pseudobulk/)

[Figure 2D](./notebooks/02_pseudobulk/)

[Figure 2E](./notebooks/02_pseudobulk/)

[Figure 2F](./notebooks/02_pseudobulk/)

**Figure 3**

[Figure 3A](./figures/03_shuffle/)

[Figure 3B-C](./notebooks/03_shuffle/)

[Figure 3D-E](./notebooks/03_shuffle/)

[Figure 3F-G](./notebooks/03_shuffle/)

**Figure 4**

[Figure 4A](./scripts/04_experiment) 

[Figure 4B](./scripts/04_experiment) 

[Figure 4C](./scripts/04_experiment)

### Supplementary Figures

**Supplemental Figure 1**

[Figure S1A](./scripts/02_pseudobulk/)

[Figure S1B](./scripts/02_pseudobulk/)

[Figure S1C](./scripts/02_pseudobulk/)

[Figure S1D](./scripts/02_pseudobulk/)

[Figure S1E](./scripts/02_pseudobulk/)

[Figure S1F](./scripts/02_pseudobulk/)

[Figure S1E](./scripts/02_pseudobulk/)

[Figure S1F](./scripts/02_pseudobulk/)

**Supplemental Figure 2**

[Supplemental Figure 2A](./notebooks/06_estimate/01_run/Rmd)

[Supplemental Figure 2B](./notebooks/09_fast/01_run_1sample.Rmd)

**Supplemental Figure 3**

[Supplemental Figure 3A](./notebooks/05_bulk/01_plot_manual_save.Rmd)

[Supplemental Figure 3B](./notebooks/05_bulk/01_plot_manual_save.Rmd) 

**Supplemental Figure 4**

[Supplemental Figure 4A](./notebooks/05_bulk/01_plot_manual_save.Rmd)

[Supplemental Figure 4B](./notebooks/05_bulk/01_plot_manual_save.Rmd)

[Supplemental Figure 4C](./notebooks/05_bulk/02_crossvalidate_save.Rmd)

[Supplemental Figure 4D](./notebooks/05_bulk/02_crossvalidate_save.Rmd)

[Supplemental Figure 4E](./notebooks/05_bulk/02_crossvalidate_save.Rmd)

**Supplemental Figure 5**

[Supplemental Figure 5](./notebooks/09_fast/03_run_sopt_realbulk_all.Rmd)


## JHPCE

The main directory path on JHPCE is located at `/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/deconvo_lute-paper`.

Further details about sample processing and QC are available at `/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/Human_DLPFC_Deconvolution`.

## Datasets

### snRNA-seq

The snRNA-seq datasets for this project currently include 9 donors and 10 samples. For the deconvolution method paper, we access these data as a `SingleCellExperiment` object. The latest (final?) version is a file called `sce_DLPFC.Rdata` located at the subdirectory:

`DLPFC_snRNAseq/processed-data/sce/`

Cell type labels for `snRNA-seq` datasets are determined from the variable `cellType_broad_hc`. This can be accessed from the `sce` object in various ways such as `colData(sce)$cellType_broad_hc` or `sce[["cellType_broad_hc"]]`. Note that the deconvolution methods paper focuses on just 6 cell types of interest, and these are identified from among the cell type labels in the `cellType_broad_hc` variable.

### Real bulk RNA-seq samples

The bulk RNA-seq datasets include 6 samples groups comprised of mRNA isolated from either nucleus, cytoplasm, or bulk, and prepped with either polyA selection or rRNA depletion strategies. Data come from 19 samples, meaning 113 total data points. 

These bulk RNA-seq data are currently being QC'd, but they will be accessible as a `SummarizedExperiment` object located at the subdirectory:

`Human_DLPFC_Deconvolution/processed-data/01_SPEAQeasy/`

### RNAscope outputs

Image data is generated for RNAscope slides by analysis with HALO and outputting the analysis results as `.csv` tables. These tables are read from the file tree located at the subdirectory:

`Human_DLPFC_Deconvolution/raw-data/HALO/Deconvolution_HALO_analysis/`

Note, this contains one subdirectory for each of the RNAscope experiments, called `Circle` and `Star`, respectively, These experiments are distinct in that they each comprise of an analysis of an independent, albeit adjacent, tissue and each includes a different set of molecular markers (see table below).

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

Details about deconvolution terminology and formulations are available in the [`lute` User's Guide](https://github.com/metamaden/lute/blob/main/vignettes/lute_users_guide.Rmd).

## Funding sources

* 1R01MH123183-01 ([grantome_link](https://grantome.com/grant/NIH/R01-MH123183-01))

## See also

* [lute](https://github.com/metamaden/lute/tree/main) : Framework for bulk transcriptomics deconvolution using cell size scale factor normalization.

