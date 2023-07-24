snrnaseq.filename <- "list-scef_markers-k2-k3-k4_ro1-dlpfc.rda"
snrnaseq.path <- file.path(base.path, snrnaseq.filename)
snrnaseq <- get(load(snrnaseq.path))
snrnaseq <- get(load())
sn.qc.k2 <- scuttle::perCellQCMetrics(lscef$k2)