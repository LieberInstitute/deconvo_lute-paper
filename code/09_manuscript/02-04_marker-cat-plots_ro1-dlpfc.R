
libv <- c("ggplot2", "ggforce", "gridExtra")
sapply(libv, library, character.only = TRUE)

set.seed(0) # seed for random colors

#----------
# load data
#----------
# ro1 dlpfc sce, markers
fname1 <- "list-scef_markers-k2-k3-k4_ro1-dlpfc.rda"
lscef1 <- get(load(fname1))

# mrb dlpfc sce, markers
fname2 <- "list-scef_markers-k2-k3-k4_mrb-dlpfc.rda"
lscef2 <- get(load(fname2))