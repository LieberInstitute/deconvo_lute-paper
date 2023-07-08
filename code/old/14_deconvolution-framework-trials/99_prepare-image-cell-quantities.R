#!/usr/bin/env R

# Author: Sean Maden
#
# Prepare cell quantities from halo image outputs.

source("deconvo_method-paper/code/14_deconvolution-framework-trials/00_parameters.R")
sapply(libv, library, character.only = T)
image.table <- get(load(halo.output.path))

# get k2 labels
image.table$k2.type <- ifelse(cell.type.vector %in% c("Inhib", "Excit"), "neuron",
                              ifelse(cell.type.vector %in% 
                                       c("Oligo", "OPC", "Astro", "Micro"), 
                                     "glial", "other"))
# filter other
image.table <- image.table[!image.table$k2.type=="other",]

# circle-star comparisons
image.count1 <- as.data.frame(table(image.table$SAMPLE_ID, image.table$k2.type))
image.count1$sample <- gsub("_.*", "", image.count1$Var1)
by.list <- list(category = paste0(image.count1$sample, ";", image.count1$Var2))
image.count1.agg1 <- aggregate(image.count1[,3], by = by.list, FUN = "mean")

image.count1.agg2 <- aggregate(image.count1[,3], by = by.list, FUN = "sum")

image.count2 <- as.data.frame(table(image.table$Sample, image.table$cell_type))


# get cell amounts on brnum
image.brnum.vector <- unique(image.table[,"BrNum"])
image.cells <- do.call(rbind, lapply(image.brnum.vector, function(brnum){
  filter.brnum <- image.table[,"BrNum"]==brnum
  counts <- image.table[filter.brnum,"cell_type"] %>% table() 
  cell.types <- names(counts)
  counts <- counts %>% as.numeric()
  proportions <- counts %>% prop.table() %>% as.numeric()
  names(counts) <- names(proportions) <- cell.types
  names(counts) <- paste0(names(counts), ".count")
  names(proportions) <- paste0(names(proportions), ".proportion")
  c(counts, proportions, brnum)
})) %>% as.data.frame()
colnames(image.cells)[15] <- "BrNum"
for(index in seq(14)){
  image.cells[,index] <- image.cells[,index] %>% as.numeric()}
# get k2 values
column.names <- colnames(image.cells)
glial.id.vector <- c("Astro", "Oligo", "OPC", "Micro")
neuron.id.vector <- c("Inhib", "Excit")
which.glial.counts <- grepl(paste0(glial.id.vector, ".count", collapse = "|"), column.names)
which.neuron.counts <- grepl(paste0(neuron.id.vector, ".count", collapse = "|"), column.names)
image.cells$glial.count <- apply(image.cells[,which.glial.counts], 1, sum)
image.cells$neuron.count <- apply(image.cells[,which.neuron.counts], 1, sum)
image.cells$total.k2 <- apply(image.cells[,c("glial.count", "neuron.count")], 1, sum)
image.cells$glial.proportion <- image.cells$glial.count/image.cells$total.k2
image.cells$neuron.proportion <- image.cells$neuron.count/image.cells$total.k2
# save
save(image.cells, file = image.cells.path)
# compare k2 estimates
cell.type.vector <- image.table$cell_type
image.table$k2.type <- ifelse(cell.type.vector %in% c("Inhib", "Excit"), "neuron",
                              ifelse(cell.type.vector %in% 
                                       c("Oligo", "OPC", "Astro", "Micro"), 
                                     "glial", "other"))
# remove "other"
image.table.filter <- image.table[!image.table$k2.type == "other",]
brnum.id.vector <- image.table.filter$BrNum
unique.brnum <- unique(brnum.id.vector)
sample.id.vector <- paste0(image.table.filter$BrNum, "_", image.table.filter$Position)
unique.sample.id <- unique(sample.id.vector)
# get prop.table
prop.table <- do.call(rbind, lapply(unique.brnum, function(brnum.iter){
  filter.brnum <- brnum.id.vector == brnum.iter
  image.brnum <- image.table.filter[filter.brnum,]
  prop.brnum <- table(image.brnum$k2.type) %>% prop.table() %>% as.data.frame()
  filter.sample.id <- grepl(paste0(brnum.iter, "_.*"), unique.sample.id)
  unique.sample.id.iter <- unique.sample.id[filter.sample.id]
  prop.sample <- do.call(rbind, lapply(unique.sample.id.iter, function(sample.id){
    filter.image.sample <- sample.id.vector == sample.id
    image.sample <- image.table.filter[filter.image.sample,]
    table(image.sample$k2.type) %>% prop.table() %>% t()
  })) %>% as.data.frame()
  prop.sample$sample.id <- unique.sample.id.iter
  prop.sample$brnum <- brnum.iter
  prop.sample$glial.prop.brnum <- prop.brnum[prop.brnum[,1]=="glial",2]
  prop.sample$neuron.prop.brnum <- prop.brnum[prop.brnum[,1]=="neuron",2]
  prop.sample
}))

# append variables for plots
prop.table$all <- "all"
prop.table$position <- gsub(".*_", "", prop.table$sample.id)
prop.table$neuron.difference <- prop.table$neuron - prop.table$neuron.prop.brnum
prop.table$neuron.abs.difference <- abs(prop.table$neuron.difference)

# scatterplot
#new.plot <- ggplot(prop.table, aes(x = neuron.prop.brnum, y = neuron)) + 
#  geom_point(alpha = 0.4) + geom_abline(slope = 1, intercept = 0, alpha = 0.5) +
#  theme_bw() + xlab("BrNum") + ylab("BrNum_Position") + ggtitle("Proportion neuron")
#new.plot + facet_wrap(~sample.id)
#new.plot + facet_wrap(~brnum)

# barplot -- brnum by difference
yaxis.label <- "Proportion difference\n(Slide - Donor)"
xaxis.label <- "Slide (BrNum_Position)"
title <- "Cell type: neuron"
new.plot <- ggplot(prop.table, aes(x = sample.id, y = neuron.difference)) + 
  geom_bar(stat = "identity") +
  theme_bw() + 
  xlab(xaxis.label) + 
  ylab(yaxis.label) + 
  ggtitle(title) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot.barplot.diff.byslide.path <- here(save.path, "ggplot-barplot_neuron-prop-diff_by-slide.jpg")
jpeg(plot.barplot.diff.byslide.path, width = 6, height = 3.5, units = "in", res = 400)
new.plot
dev.off()

# barplot -- brnum by abs difference
yaxis.label <- "Proportion absolute difference\nabs(Slide - Donor)"
xaxis.label <- "Slide (BrNum_Position)"
title <- "Cell type: neuron"
new.plot <- ggplot(prop.table, aes(x = sample.id, y = neuron.abs.difference)) + 
  geom_bar(stat = "identity") +
  theme_bw() + 
  xlab(xaxis.label) + 
  ylab(yaxis.label) + 
  ggtitle(title) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
jpeg(plot.barplot.absdiff.byslide.path, width = 6, height = 3.5, units = "in", res = 400)
new.plot
dev.off()

# set jitterbox params
ylab.diff <- "Proportion difference\n(Slide - Donor)"
ylab.absdiff <- "Absolute proportion difference\nabs(Slide - Donor)"
title <- "Cell type: neuron"

# jitterbox -- group: position
xaxis.label <- "Position"
# difference
new.plot <- ggplot(prop.table, aes(x = position, y = neuron.difference)) + 
  geom_jitter(alpha = 0.5) + geom_boxplot(alpha = 0, color = "cyan") +
  theme_bw() + 
  xlab(xaxis.label) + 
  ylab(ylab.diff) + 
  ggtitle(title) +
  geom_hline(yintercept = 0) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
jpeg(plot.jitterbox.diff.byposition.path, width = 3, height = 2.5, units = "in", res = 400)
new.plot
dev.off()
# absolute difference
new.plot <- ggplot(prop.table, aes(x = position, y = neuron.abs.difference)) + 
  geom_jitter(alpha = 0.5) + geom_boxplot(alpha = 0, color = "cyan") +
  theme_bw() + 
  xlab(xaxis.label) + 
  ylab(ylab.absdiff) + 
  ggtitle(title) +
  geom_hline(yintercept = 0) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
jpeg(plot.jitterbox.absdiff.byposition.path, width = 3, height = 2.5, units = "in", res = 400)
new.plot
dev.off()

# jitterbox -- group: donor
xaxis.label <- "Donor (BrNum)"
# difference
new.plot <- ggplot(prop.table, aes(x = brnum, y = neuron.difference)) + 
  geom_jitter(alpha = 0.5) + geom_boxplot(alpha = 0, color = "cyan") +
  theme_bw() + 
  xlab(xaxis.label) + 
  ylab(ylab.diff) + 
  ggtitle(title) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, alpha = 0.5)
jpeg(plot.jitterbox.diff.bydonor.path, width = 5, height = 2.5, units = "in", res = 400)
new.plot
dev.off()
# absolute difference
new.plot <- ggplot(prop.table, aes(x = brnum, y = neuron.abs.difference)) + 
  geom_jitter(alpha = 0.5) + geom_boxplot(alpha = 0, color = "cyan") +
  theme_bw() + 
  xlab(xaxis.label) + 
  ylab(ylab.absdiff) + 
  ggtitle(title) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, alpha = 0.5)
jpeg(plot.jitterbox.absdiff.bydonor.path, width = 5, height = 2.5, units = "in", res = 400)
new.plot
dev.off()

# jitterbox -- group: all
xaxis.label <- "All"
# difference
new.plot <- ggplot(prop.table, aes(x = all, y = neuron.difference)) + 
  geom_jitter(alpha = 0.5) + geom_boxplot(alpha = 0, color = "cyan") +
  theme_bw() + 
  xlab(xaxis.label) + 
  ylab(ylab.diff) + 
  ggtitle(title) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, alpha = 0.5)
jpeg(plot.jitterbox.diff.all.path, width = 2.5, height = 2.5, units = "in", res = 400)
new.plot
dev.off()
# absdiff
new.plot <- ggplot(prop.table, aes(x = all, y = neuron.abs.difference)) + 
  geom_jitter(alpha = 0.5) + geom_boxplot(alpha = 0, color = "cyan") +
  theme_bw() + 
  xlab(xaxis.label) + 
  ylab(ylab.absdiff) + 
  ggtitle(title) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, alpha = 0.5)
jpeg(plot.jitterbox.absdiff.all.path, width = 2.5, height = 2.5, units = "in", res = 400)
new.plot
dev.off()
