[![DOI](https://zenodo.org/badge/543720810.svg)](https://zenodo.org/doi/10.5281/zenodo.10568795)

# deconvo_lute-paper

Supplementary materials to prepare data and run analyses for the paper "lute: estimating the cell composition of heterogeneous tissue with varying cell sizes using gene expression" (1). The `lute` R/Bioconductor package (2) is publicly available for download.

## Citation

[Maden, Sean, Louise A. Huuki-Myers, Sang Ho Kwon, Leonardo Collado-Torres, Kristen R. Maynard, Stephanie C. Hicks. "lute: estimating the cell composition of heterogeneous tissue with varying cell sizes using gene expression" (manuscript in preparation)]()

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

Folders `huuki-meyers-et-al-2023_dlpfc-cohort1` (3,4), `tran-et-al-2021_dlpfc-cohort2` (5), include data and analyses from neurotypical subject dorsolateral prefrontal cortex (DLPFC). Folder `monaco-et-al-2019_pbmc-cohort1` (6) includes data and analyses of healthy subject peripheral blood mononuclear cell (PBMC) data.

### Folders and directory structure

In addition to the dataset folders, `software_manuscript` includes additional files for the `lute` software vignettes and reporting RMSEs. `dann-et-al-2023_env` and `sosina-et-al-2021_env` include environment information from studies (7,8).

The file tree structure of this repo returned from `tree -d`:

```
├conda_env
├── dann-et-al-2023_env
│   └── yml
├── huuki-meyers-et-al-2023_dlpfc-cohort1
│   ├── data
│   │   └── 01_mae
│   ├── env
│   │   ├── 02_pseudobulk
│   │   ├── 02_summarize_mae
│   │   ├── 03_shuffle
│   │   ├── 04_experiment
│   │   └── 07_adjustment
│   ├── figures
│   │   ├── 02_pseudobulk
│   │   ├── 03_shuffle
│   │   └── 07_adjustment
│   ├── notebooks
│   │   ├── 00_preprocess
│   │   ├── 01_mae
│   │   ├── 02_pseudobulk
│   │   ├── 02_summarize_mae
│   │   ├── 03_shuffle
│   │   └── 07_adjustment
│   ├── outputs
│   │   ├── 00_preprocess
│   │   ├── 01_mae
│   │   ├── 02_pseudobulk
│   │   ├── 02_summarize_mae
│   │   ├── 03_shuffle
│   │   └── 07_adjustment
│   ├── scripts
│   │   ├── 00_preprocess
│   │   ├── 01_mae
│   │   ├── 02_pseudobulk
│   │   ├── 02_summarize_mae
│   │   ├── 03_shuffle
│   │   └── 07_adjustment
│   ├── source
│   └── unit_tests
│       ├── 02_mae_summaries
│       ├── 03_shuffle
│       └── 05_bulk
├── monaco-et-al-2019_pbmc-cohort1
│   ├── data
│   │   ├── cellScaleFactors
│   │   └── monaco_et_al_2019
│   │       ├── geo
│   │       ├── manuscript
│   │       │   └── original
│   │       └── shared
│   ├── env
│   │   └── 01_pseudobulk
│   ├── figures
│   │   └── 01_pseudobulk
│   ├── notebooks
│   │   └── 01_pseudobulk
│   │       └── 04_simulation_files
│   │           └── figure-latex
│   ├── outputs
│   │   └── 01_pseudobulk
│   ├── scripts
│   │   └── 01_pseudobulk
│   ├── source
│   └── yml
├── software_manuscript
│   ├── data
│   ├── env
│   │   └── 01_multipanel
│   ├── figures
│   │   ├── 01_multipanel
│   │   └── 02_diagram
│   ├── notebooks
│   │   ├── 01_multipanel
│   │   └── 02_aggregateRMSE
│   ├── outputs
│   │   └── 02_aggregateRMSE
│   ├── scripts
│   │   ├── 01_multipanel
│   │   ├── 02_diagram
│   │   └── 03_rmseTable
│   └── src
├── sosina-et-al-2021_env
│   ├── notebooks
│   ├── scripts
│   └── yml
└── tran-et-al-2021_dlpfc-cohort2
    ├── data
    ├── env
    │   └── 01_pseudobulk
    ├── figures
    │   └── 01_pseudobulk
    ├── notebooks
    │   └── 01_pseudobulk
    ├── outputs
    │   ├── 01_pseudobulk
    │   └── 02_summaries
    ├── scripts
    │   ├── 01_pseudobulk
    │   └── 02_summaries
    └── yml
```

The folder contents are described as follows:

* `huuki-meyers-et-al-2023_dlpfc-cohort1` : Materials for analyses of principal DLPFC dataset.
    * `~/data` : External data used by supplement.
    * `~/env` : Environment (`.RData`) saves from `save.image()`.
    * `~/figures` : Figure saves from `~/notebooks`.
    * `~/notebooks` : Notebooks (`.Rmd`) which load saved images and save figures.
    * `~/outputs` : Objects (`.rda`) saved by scripts.
    * `~/scripts` : Scripts (`.Rmd`) performing data preparation and analyses.
    * `~/source` : Scripts (`.R`) shared across analysis sections.
    * `~/unit_tests` : Scripts (`.R`) testing sample and sample source filters across analysis sections.

* `tran-et-al-2021_dlpfc-cohort2` : Materials for analyses of validation DLPFC dataset.

* `monaco-et-al-2019_pbmc-cohort1` : Materials for analyses of peripheral blood mononuclear cells (PBMC) dataset.

* `sosina-et-al-2021_env` : Supplementary materials for analyses of previous work.

* `dann-et-al-2023_env` : Supplementary materials for preprocessing samples for running with `lute`.

### Manuscript Figures

Paths to manuscript figure sources.

* Figure 1 : (see manuscript)

* Figure 2A : `huuki-meyers-et-al-2023_dlpfc-cohort1/notebooks/02_pseudobulk/01_k2.Rmd`

* Figure 2B : `cohort3/notebooks/01_pseudobulk/04_simulation.Rmd"`

* Figure 3A : `huuki-meyers-et-al-2023_dlpfc-cohort1/notebooks/07_adjustment/
03_run_sopt_realbulk_all.Rmd`

* Figure 3B : `huuki-meyers-et-al-2023_dlpfc-cohort1/notebooks/07_adjustment/03_run_sopt_realbulk_all.Rmd`

* Figure S1A : `huuki-meyers-et-al-2023_dlpfc-cohort1/notebooks/02_pseudobulk/01_k2.Rmd`

* Figure S1B : `tran-et-al-2021_dlpfc-cohort2/notebooks/01_pseudobulk/01_k2_mrb.Rmd`

* Figure S2A : `huuki-meyers-et-al-2023_dlpfc-cohort1/notebooks/02_pseudobulk/02_k3.Rmd`

* Figure S2B : `huuki-meyers-et-al-2023_dlpfc-cohort1/notebooks/02_pseudobulk/02_k3.Rmd`

* Figure S2C : `tran-et-al-2021_dlpfc-cohort2/notebooks/01_pseudobulk/02_k3_mrb.Rmd`

* Figure S2D : `huuki-meyers-et-al-2023_dlpfc-cohort1/notebooks/01_pseudobulk/02_k3_mrb.Rmd`

* Figure S3A : `huuki-meyers-et-al-2023_dlpfc-cohort1/notebooks/03_shuffle/figS3A.Rmd`

* Figure S3B : `huuki-meyers-et-al-2023_dlpfc-cohort1/notebooks/03_shuffle/figS3B.Rmd`

* Figure S4 : `monaco-et-al-2019_pbmc-cohort1/notebooks/01_pseudobulk/04_simulations.Rmd"`

* Figure S5 : `huuki-meyers-et-al-2023_dlpfc-cohort1/notebooks/07_adjustment/
03_run_sopt_realbulk_all.Rmd`

* Figure S6 : `huuki-meyers-et-al-2023_dlpfc-cohort1/notebooks/07_adjustment/03_run_sopt_realbulk_all.Rmd`

### Manuscript Tables

Paths to manuscript table sources.

* Table S1 : `huuki-meyers-et-al-2023_dlpfc-cohort1/notebooks/02_summarize_mae/TABLE_S1.csv`

* Table S2 : `huuki-meyers-et-al-2023_dlpfc-cohort1/outputs/02_summarize_mae/TABLE_S2.csv`

* Table S3 : `software/outputs/02_aggregateRMSE/TABLE_S3.csv`

* Table S4 : `tran-et-al-2021_dlpfc-cohort2/outputs/02_summaries/TABLE_S4.csv`

* Table S5 : `tran-et-al-2021_dlpfc-cohort2/outputs/02_summaries/cohort2_nuclei_summaries.csv`

* Table S6 : `huuki-meyers-et-al-2023_dlpfc-cohort1/outputs/08_sizes/tables_cell_sizes.csv`

* Table S7 : `huuki-meyers-et-al-2023_dlpfc-cohort1/outputs/03_shuffle/TABLE_S7.csv`

* Table S8 : `huuki-meyers-et-al-2023_dlpfc-cohort1/outputs/02_summarize_mae/TABLE_S8.csv`

### Compute environment

Files to repeat analyses have been provided as `.yml` files in the `./deconvo_lute-paper/yml` folder. These may be used to generate new conda environments. Additionally, top-level folders with `.*_env` in their name contain additional compute environment files for manuscript analyses.

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

	1. Maden, Sean, Louise A. Huuki-Myers, Sang Ho Kwon, Leonardo Collado-Torres, Kristen R. Maynard, Stephanie C. Hicks. "lute: estimating the cell composition of heterogeneous tissue with varying cell sizes using gene expression" (manuscript in preparation)

    2. Maden, Sean K., Stephanie C. Hicks (2024). lute: Framework for cell size scale factor normalized bulk transcriptomics deconvolution experiments. doi:10.18129/B9.bioc.lute, R package version 0.99.25, DOI: 10.18129/B9.bioc.lute

    3. Huuki-Myers, Louise, Abby Spangler, Nick Eagles, Kelsey D. Montgomery, Sang Ho Kwon, Boyi Guo, Melissa Grant-Peters, Heena R. Divecha, Madhavi Tippani, Chaichontat Sriworarat, Annie B. Nguyen, Prashanthi Ravichandran, Matthew N. Tran, Arta Seyedian, PsychENCODE consortium, Thomas M. Hyde, Joel E. Kleinman, Alexis Battle, Stephanie C. Page, Mina Ryten, Stephanie C. Hicks, Keri Martinowich, Leonardo Collado-Torres, and Kristen R. Maynard. 2023. “Integrated Single Cell and Unsupervised Spatial Transcriptomic Analysis Defines Molecular Anatomy of the Human Dorsolateral Prefrontal Cortex.” bioRxiv. https://doi.org/10.1101/2023.02.15.528722.

    4. Huuki-Myers, Louise A., Kelsey D. Montgomery, Sang Ho Kwon, Stephanie Cinqueman, Sean K. Maden, Nicholas J. Eagles, Joel E. Kleinman, Thomas M. Hyde, Stephanie C. Hicks, Kristen R. Maynard, Leonardo Collado-Torres. "Benchmark of cellular deconvolution methods using a multi-assay reference dataset from postmortem human prefrontal cortex". BioRxiv. 2024 Feb 12. doi: https://doi.org/10.1101/2024.02.09.579665

    5. Tran, Matthew N., Kristen R. Maynard, Abby Spangler, Louise A. Huuki, Kelsey D. Montgomery, Vijay Sadashivaiah, Madhavi Tippani, et al. 2021. “Single-Nucleus Transcriptome Analysis Reveals Cell-Type-Specific Molecular Signatures across Reward Circuitry in the Human Brain.” Neuron 109 (19): 3088-3103.e5. https://doi.org/10.1016/j.neuron.2021.09.001.

    6. Monaco, Gianni, Bernett Lee, Weili Xu, Seri Mustafah, You Yi Hwang, Christophe Carré, Nicolas Burdin, Lucian Visan, Michele Ceccarelli, Michael Poidinger, Alfred Zippelius, João Pedro de Magalhães, Anis Larbi. 2019. "RNA-Seq Signatures Normalized by mRNA Abundance Allow Absolute Deconvolution of Human Immune Cell Types." Cell Reports 26 (6): 1627-1640.e7. https://doi.org/10.1016/j.celrep.2019.01.041.

    7. Sosina, Olukayode A., Matthew N. Tran, Kristen R. Maynard, Ran Tao, Margaret A. Taub, Keri Martinowich, Stephen A. Semick, Bryan C. Quach, Daniel R. Weinberger, Thomas Hyde, Dana B. Hancock, Joel E. Kleinman, Jeffrey T. Leek, Andrew E. Jaffe. Strategies for cellular deconvolution in human brain RNA sequencing data. F1000Res. 2021 Aug 4;10:750. doi.org/10.12688/f1000research.50858.1

    8. Dann, Emma, Ana-Maria Cujba, Amanda J. Oliver, Kerstin B. Meyer, Sarah A. Teichmann, and John C. Marioni. 2023. “Precise Identification of Cell States Altered in Disease Using Healthy Single-Cell References.” Nature Genetics 55 (11): 1998–2008. https://doi.org/10.1038/s41588-023-01523-7.
