# deconvo_method-paper/cohort3

Experiment and results using previously published peripheral blood mononuclear cells (PBMCs) profiled from 17 normal/healthy young adult donors (1). Data were accessed from either the main publication, publication supplement, GEO, or provided at author request.

## Folders

Descriptions of directories contained in `cohort3`.

### Folder tree

The file tree for the `cohort3` directory looks like:

```
cohort3
├── README.md
├── data
│   └── monaco_et_al_2019
│       ├── geo
│       ├── manuscript
│       └── shared
├── env
│   └── 01_pseudobulk
├── figures
├── notebooks
│   └── 01_pseudobulk
├── outputs
│   └── 01_pseudobulk
├── scripts
│   └── 01_pseudobulk
├── source
└── yml
```

Folder descriptions

* `data` : Primary publication data files, supplement, and metadata, including RNA-seq expression and flow cytometry cell abundances.

* `env` : Environments and objects resulting from running `script` experiments.

* `figures` : Figures stored from running `script` experiments.

* `notebooks` : Notebook built from stored `env` files, outputs from running `script` experiments.

* `outputs` : Ouput data files from running `scripts`.

* `scripts` : Scripts to perform experiments and analyses. Results saved to `files`, `env`, and `outputs`. Visualizations built from `notebook` files.

* `source` : Common functions shared across scripts, including to read experiment files from `data` and perform pseudobulk simulation experiment runs.

* `yml` : Resources to make the virtual environment for executing scripts and building notebooks.

## Citation

1. Monaco, Gianni, Bernett Lee, Weili Xu, Seri Mustafah, You Yi Hwang, Christophe Carré, Nicolas Burdin, et al. 2019. “RNA-Seq Signatures Normalized by mRNA Abundance Allow Absolute Deconvolution of Human Immune Cell Types.” Cell Reports 26 (6): 1627-1640.e7. [https://doi.org/10.1016/j.celrep.2019.01.041](https://doi.org/10.1016/j.celrep.2019.01.041).
