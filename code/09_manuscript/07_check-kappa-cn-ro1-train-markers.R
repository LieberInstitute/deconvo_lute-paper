load("~/GitHub/deconvo_method-paper/outputs/09_manuscript/list-scef_markers-k2-k3-k4_mrb-dlpfc.rda")
names(lscef)
# [1] "k2" "k3" "k4"
kappa(assays(lscef[["k2"]])[["counts"]])
# [1] 51.10824
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