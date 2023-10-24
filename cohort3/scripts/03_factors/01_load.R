#!/usr/bin/env R

# Author: Sean Maden
#
# Load cell scale factors.
#

#-----
# load
#-----
monaco.data <- get(load("./data/02_abisseq/csf_table.rda"))
load("./env/02_abisseq/01_abisseq_script.RData")

#---------------
# combine/append
#---------------
# map cell labels
monaco.data <- monaco.data[grepl("Monaco", monaco.data$citation.s.),]
vector.ref.size <- monaco.data$scale.factor.value
names(vector.ref.size) <- monaco.data$cell_type

# map values
df.size <- df.tall[,c(10,11)]
df.size$s.cell.size.ref <- "NA"
df.size[grepl("T\\.CD4\\..*", df.size$cell.type),]$s.cell.size.ref <-
  vector.ref.size["T cells CD4"]
df.size[grepl("T\\.CD8\\..*", df.size$cell.type),]$s.cell.size.ref <-
  vector.ref.size["T cells CD8"]
df.size[grepl("B\\..*", df.size$cell.type),]$s.cell.size.ref <-
  vector.ref.size["B cells"]
df.size[grepl("Monocytes", df.size$cell.type),]$s.cell.size.ref <-
  vector.ref.size["Monocytes"]
df.size[grepl("Neutrophils\\..*", df.size$cell.type),]$s.cell.size.ref <-
  vector.ref.size["Neutrophils"]
df.size[grepl("NK", df.size$cell.type),]$s.cell.size.ref <-
  vector.ref.size["NK cells"]
df.size[grepl("Plasma.*", df.size$cell.type),]$s.cell.size.ref <-
  vector.ref.size["Plasma cells"]
df.size[grepl(".*DCs", df.size$cell.type),]$s.cell.size.ref <-
  vector.ref.size["Dendritic cells"]
df.size$s.cell.size.ref <- as.numeric(df.size$s.cell.size.ref)

# correlation test
cor.test <- cor.test(df.size$s.cell.size.log, df.size$s.cell.size.ref)

#-----
# save
#-----
save.image("./03_factors/01_load_script.RData")
