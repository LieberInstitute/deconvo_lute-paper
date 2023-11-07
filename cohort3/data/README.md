# data

External data for cohort analysis. Data comes from Monaco et al 2019 Cell Reports online materials.

# Files

## Tree

Results from running `tree`:

```
README.md
├── cellScaleFactors
│   └── csf_table.rda
└── monaco_et_al_2019
    ├── geo
    │   └── GSE107011_Processed_data_TPM.txt
    ├── manuscript
    │   ├── abisseq_rnaseq_cell_types_references_k17.csv
    │   ├── flow_cytometry_true_cell_type_proportions.csv
    │   └── original
    └── shared
        └── mRNAyield_PBMC.txt
```

## Contents

* `csf_table.rda` : cell scale factors table from `cellScaleFactors` package.

* `monaco_et_al_2019` : datasets published and affiliated with the Monaco et al 2019 Cell Reports study.

* `monaco_et_al_2019/geo/GSE107011_Processed_data_TPM.txt` : 
expression table from GEO record [GSE107011](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107011).

* `abisseq_rnaseq_cell_types_references_k17.csv` : ABIS-seq RNA-seq cell type reference table from study supplemental tables (Table S5).

* `flow_cytometry_true_cell_type_proportions.csv` : Flow Cytometry true cell type proportions from study supplemental tables (Table S6).

* `original` : Original tables published from the study.

* `mRNAyield_PBMC.txt` : Flow cytometry cell abundances from FCS, standard file format for flow cytometry.
