#---------------
# get cell sizes
#---------------
# get cell size factor series
# load cell size scale factors
df.csf <- get_csf_reference()
df.csf <- df.csf[df.csf$tissue=="brain",]
df.csf.area <- df.csf[grepl("mRNA", df.csf$scale.factor.type),]
df.csf.mrna <- df.csf[df.csf$scale.factor.type=="cell area",]
s.set.osm.area <- c("glial" = df.csf.area[df.csf.area$cell_type=="glial",]$scale.factor.value,
                    "neuron" = df.csf.area[df.csf.area$cell_type=="neuron",]$scale.factor.value)
s.set.osm.mrna <- c("glial" = df.csf.mrna[df.csf.mrna$cell_type=="glial",]$scale.factor.value,
                    "neuron" = df.csf.mrna[df.csf.mrna$cell_type=="neuron",]$scale.factor.value)
# coerce to list for experiment
#list.s.pred <- list(s.set.null = c("glial" = 1, "neuron" = 1),
#                    s.set1 = c("glial" = 0.5, "neuron" = 10),
#                    s.set2 = c("glial" = 1, "neuron" = 10),
#                    s.set3 = c("glial" = 2, "neuron" = 10),
#                    s.set4 = c("glial" = 3, "neuron" = 10),
#                    s.set5 = c("glial" = 4, "neuron" = 10))

list.s.pred <- list(s.set4 = c("glial" = 3, "neuron" = 10))

# notes:
# difference scale seems to matter more than magnitude scale (e.g. similar c(3,10), c(30,100), etc.)
# ratio difference matters (e.g. different c(3,10), c(5,10), etc.)