lm <- list.markers.final
markers.neuron <- lm$neuron[,1]
markers.glial <- lm$glial[,1]
length(unique(markers.neuron, markers.glial)) # 7703
length(intersect(markers.neuron, markers.glial)) # 3693
# 3693/7703 # 0.4794236

# get type-overlapping markers
type.overlapping.markers <- intersect(markers.neuron, markers.glial)

# get overlap with prior markers
length(intersect(type.overlapping.markers, rownames(lscef$k2))) # [1] 36
length(intersect(type.overlapping.markers, rownames(lscef$k3))) # [1] 35
length(intersect(type.overlapping.markers, rownames(lscef$k4))) # [1] 61