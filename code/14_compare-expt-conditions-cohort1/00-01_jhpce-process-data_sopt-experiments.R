

# Author: Sean Maden 
# 
# Run experiment script on JHPCE. Get de novo (i.e. not batch-robust) markers for given K.
#
#

libv <- c("lute")
sapply(libv, library, character.only = T)

#-----
# load
#-----
# get assay data
# load sce data
sce.path <- ""
sce <- get(load(sce.path))
# load bulk data
# y.unadj <- mae[["bulk.rnaseq"]]
y.path <- ""
y.unadj <- get(load(y.path))
rownames(y.unadj) <- rowData(y.unadj)[,rd.geneid.var.y]
y.unadj <- y.unadj[rownames(y.unadj) %in% rownames(sce),]
y.unadj <- y.unadj[order(match(rownames(y.unadj), rownames(sce))),]
y.unadj <- logNormCounts(y.unadj)
y.unadj$sample.id <- y.unadj$batch.id2
y.train <- y.unadj[,!y.unadj$BrNum %in% validation.sample.id]
y.validate <- y.unadj[,y.unadj$BrNum %in% validation.sample.id]
# load rnascope data
# df.rnascope.kdata
df.rn.path <- ""
df.rn <- get(load(df.rn.path))

# source scripts
# source functions
# source deconvo buddies info

#---------------------------
# prep experiment parameters
#---------------------------
assay.name.lute <- "logcounts"
group.id.variable <- group.name.experiment <- "sample.id"

colData(sce)[[group.id.variable]] <- sce[["BrNum"]]
colData(y.unadj)[[group.id.variable]] <- y.unadj[["Sample"]]
sample.id.vector <- unique(sce[[group.id.variable]])

# assay parameters
assay.name.lute <- "logcounts"
group.name.experiment <- "sample.id"
validation.sample.id <- c("Br6432", "Br6522", "Br8667")
rd.geneid.var.y <- "Symbol"
rd.geneid.var.sce <- "gene_name"
y.group.variable.name <- "batch.id2"

#------------------
# 1. K2 Experiments
#------------------
# Prep k2 parameters
# get list.df.true, rnascope data
celltype.variable <- k.variable.name <- "k2"
unique.types <- c("glial", "neuron")
sample.id.vector <- unique(y.unadj$group.id)
list.df.true <- df.true.list(df.rn, sample.id.vector, celltype.variable, unique.types)

# 1C. Experiment: 
# * P_true : RNAscope proportions
# * G : De novo K2 markers (N = 20 genes/celltype)
# * Z : snRNAseq matched or unmatched

k2.1c.result.list <- kmatch_experiment(k.variable.name = k.variable.name,
                                        sce = sce, 
                                        sample.id.vector = sample.id.vector, 
                                        list.df.true = list.df.true, 
                                        y.eset = y.unadj, y.train = y.train, 
                                        y.validate = y.validate, 
                                        assay.name = assay.name.lute, 
                                        celltype.variable = k.variable.name, 
                                        group.name = group.name.experiment)

# 1D. Experiment: 
# * P_true : snRNAseq proportions
# * G : De novo K2 markers (N = 20 genes/celltype)
# * Z : snRNAseq matched or unmatched

k2.1d.result.list <- kmatch_experiment(k.variable.name = k.variable.name,
                                       sce = sce, 
                                       sample.id.vector = sample.id.vector, 
                                       list.df.true = list.df.true, 
                                       y.eset = y.unadj, y.train = y.train, 
                                       y.validate = y.validate, 
                                       assay.name = assay.name.lute, 
                                       celltype.variable = k.variable.name, 
                                       group.name = group.name.experiment)

# save results
k2.1c.res.path <- ""
k2.1d.res.path <- ""
save(k2.1c.result.list, file = k2.1c.res.path)
save(k2.1d.result.list, file = k2.1c.res.path)

#------------------
# 1. K3 Experiments
#------------------
