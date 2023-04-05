load("~/GitHub/deconvo_method-paper/outputs/09_manuscript/list-scef_markers-k2-k3-k4_mrb-dlpfc.rda")
names(lscef)
# [1] "k2" "k3" "k4"
assays(lscef[["k2"]])[["counts"]] %>% kappa()
# [1] 51.10824

# try subsampling 31 markers
set.seed(0)
assays(lscef[["k2"]])[["counts"]][sample(40, 31),] %>% kappa()
# [1] 22.99347
assays(lscef[["k2"]])[["counts"]][sample(40, 31),] %>% kappa()
# [1] 54.70324
assays(lscef[["k2"]])[["counts"]][sample(40, 31),] %>% kappa()
# [1] 44.65265
assays(lscef[["k2"]])[["counts"]][sample(40, 31),] %>% kappa()
# [1] 30.36579

kappa(assays(lscef[["k3"]])[["counts"]])
# [1] 58.63972
kappa(assays(lscef[["k4"]])[["counts"]])
# [1] 186.883
names(assays(lscef[[1]]))
# [1] "counts"     "counts_adj"
kappa(assays(lscef[["k2"]])[["counts_adj"]])
# [1] 52.12258
kappa(assays(lscef[["k3"]])[["counts_adj"]])
# [1] 62.93934
kappa(assays(lscef[["k4"]])[["counts_adj"]])
# [1] 199.0362

# kappa values for signature matrices
lscef[["k2"]] %>% signature_matrix_from_sce() %>% kappa()
# [1] 1.170085
lscef[["k3"]] %>% signature_matrix_from_sce() %>% kappa()
# [1] 1.316566
lscef[["k4"]] %>% signature_matrix_from_sce() %>% kappa()
# [1] 1.666715
lscef[["k2"]] %>% signature_matrix_from_sce(assay.name = "counts_adj") %>% kappa()
# [1] 1.208293
lscef[["k3"]] %>% signature_matrix_from_sce(assay.name = "counts_adj") %>% kappa()
# [1] 1.36493
lscef[["k4"]] %>% signature_matrix_from_sce(assay.name = "counts_adj") %>% kappa()
# [1] 1.679446
