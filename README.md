[![DOI](https://zenodo.org/badge/543720810.svg)](https://zenodo.org/doi/10.5281/zenodo.10568795)

# deconvo_lute-paper

Supplementary materials to prepare data and run analyses for the paper "lute: estimating the cell composition of heterogeneous tissue with varying cell sizes using gene expression". The `lute` (1) R package is available from Bioconductor.

## Citation

Sean K. Maden, Louise A. Huuki-Myers, Sang Ho Kwon, Leonardo Collado-Torres, Kristen R. Maynard, Stephanie C. Hicks. "lute: estimating the cell composition of heterogeneous tissue with varying cell sizes using gene expression". [link]()

## Reproducing analyses

Scripts are provided in clearly-labeled folder. Script and folder names include numbers indicating their intended run order. Environment (`.RData`) and output (`.rda`) files may be generated by navigating one level above the `scripts` folder and entering the following from a shell prompt:

```
Rscript -e "source('PATH_TO_SCRIPT')"
```

or the following from an R session or script:

```
source('PATH_TO_SCRIPT')
```

Where `PATH_TO_SCRIPT` is the valid script path (e.g. `./scripts/02_pseudobulk/01_k2.R`)

Figures and notebook builds are generated by building the `.Rmd` files in the `notebooks` folder.

### Datasets

Folders `cohort1`, `cohort2`, and `cohort3` include bulk and single-nucleus transcriptomics data from 3 sources including dorsolateral prefrontal cortex (DLPFC) from (1) and (2), and peripheral blood mononuclear cell data from (3).

### Folders

The file tree structure of this repo, returned from `tree deconvo_method-paper -d`, is:

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
│   │   ├── 07_adjustment
│   │   └── 08_sizes
│   ├── figures
│   │   ├── 02_pseudobulk
│   │   ├── 03_shuffle
│   │   ├── 04_experiment
│   │   ├── 05_bulk
│   │   ├── 06_estimate
│   │   └── 08_adjustment
│   ├── notebooks
│   │   ├── 00_preprocess
│   │   ├── 01_mae
│   │   ├── 02_pseudobulk
│   │   ├── 03_shuffle
│   │   ├── 04_experiment
│   │   ├── 05_bulk
│   │   ├── 06_estimate
│   │   ├── 07_adjustment
│   │   ├── 08_sizes
│   │   └── bioc2023
│   ├── outputs
│   │   ├── 00_preprocess
│   │   ├── 01_mae
│   │   ├── 02_pseudobulk
│   │   ├── 03_shuffle
│   │   ├── 04_experiment
│   │   ├── 05_bulk
│   │   ├── 06_estimate
│   │   ├── 07_adjustment
│   │   └── 08_sizes
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
│   │   ├── 07_adjustment
│   │   └── 08_sizes
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

The folder contents are described as follows:

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

### Compute environment

Files to repeat analyses have been provided as `.yml` files in the `./deconvo_lute-paper/yml` folder. These may be used to generate new conda environments.

## Dataset descriptions

### Single-nucleus RNA-seq (snRNA-seq) data

Single-nucleus RNA-seq (snRNA-seq) data from DLPFC were generated and preprocessed as described in Huuki-Meyers et al 2023. `SingleCellExperiment` files of these data are available from ![LieberInstitute/DLPFC_snRNAseq](https://github.com/LieberInstitute/DLPFC_snRNAseq). Cell type labels for `snRNA-seq` datasets are determined from the variable `cellType_broad_hc`. This can be accessed from the `sce` object in various ways such as `colData(sce)$cellType_broad_hc`. Note that the deconvolution methods paper focuses on just 6 cell types of interest, and these are identified from among the cell type labels in the `cellType_broad_hc` variable.

### Bulk RNA-seq data

The bulk RNA-seq datasets include 6 samples groups comprised of mRNA isolated from either nucleus, cytoplasm, or bulk, and prepped with either polyA selection or rRNA depletion strategies. Data come from 19 samples, meaning 113 total data points. Bulk RNA-seq data are available accessible as `SummarizedExperiment` files located at the project folder path:

`Human_DLPFC_Deconvolution/processed-data/01_SPEAQeasy/`

### RNAscope data

Image data is generated for RNAscope slides by analysis with HALO and outputting the analysis results as `.csv` tables. These tables are read from the file tree located at the project folder path:

`Human_DLPFC_Deconvolution/raw-data/HALO/Deconvolution_HALO_analysis/`

Note, this contains one subdirectory for each of the RNAscope experiments, called `Circle` and `Star`, respectively. The processed data with confidence annotations files are located at the project folder path:

`Human_DLPFC_Deconvolution/processed-data/03_HALO/halo_all.Rdata`

## Funding sources

* 1R01MH123183-01 ([grantome_link](https://grantome.com/grant/NIH/R01-MH123183-01))

## Works cited

1. Maden SK, Hicks S (2023). lute: Framework for cell size scale factor normalized bulk transcriptomics deconvolution experiments. doi:10.18129/B9.bioc.lute, R package version 0.99.24, https://bioconductor.org/packages/lute.

2. Huuki-Myers, Louise, Abby Spangler, Nick Eagles, Kelsey D. Montgomery, Sang Ho Kwon, Boyi Guo, Melissa Grant-Peters, et al. 2023. “Integrated Single Cell and Unsupervised Spatial Transcriptomic Analysis Defines Molecular Anatomy of the Human Dorsolateral Prefrontal Cortex.” bioRxiv. https://doi.org/10.1101/2023.02.15.528722.

3. Tran, Matthew N., Kristen R. Maynard, Abby Spangler, Louise A. Huuki, Kelsey D. Montgomery, Vijay Sadashivaiah, Madhavi Tippani, et al. 2021. “Single-Nucleus Transcriptome Analysis Reveals Cell-Type-Specific Molecular Signatures across Reward Circuitry in the Human Brain.” Neuron 109 (19): 3088-3103.e5. https://doi.org/10.1016/j.neuron.2021.09.001.

4. Monaco, Gianni, Bernett Lee, Weili Xu, Seri Mustafah, You Yi Hwang, Christophe Carré, Nicolas Burdin, et al. 2019. “RNA-Seq Signatures Normalized by mRNA Abundance Allow Absolute Deconvolution of Human Immune Cell Types.” Cell Reports 26 (6): 1627-1640.e7. https://doi.org/10.1016/j.celrep.2019.01.041.
