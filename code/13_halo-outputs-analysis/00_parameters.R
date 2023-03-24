#!/usr/bin/env R

# Author: Sean Maden
#
# Parameters for HALO image analyses.
#

#-------
# header
#------- 

libv <- c("here", "nlme", "ggplot2", "gridExtra", "dplyr", "ggforce")
sapply(libv, library, character.only = TRUE)

# set the save directory
save.path <- here("deconvo_method-paper", "outputs", "13_halo-outputs-analysis")

# set the halo data path
halo.output.file.name <- "halo_all.Rdata"
halo.output.path <- here("Human_DLPFC_Deconvolution", "processed-data", "03_HALO", 
                         halo.output.file.name)

# cell labels
labels <- c("Endo" = "CLDN5", "Astro" = "GFAP", "Inhib" = "GAD1", "Excit" = "SLC17A7", 
            "Micro" = "TMEM119", "Oligo" = "OLIG2")

# helper functions
normalization1 <- function(variable){log10(variable)}

summary_term_group <- function(table, term, cell.type.label){
  aggregate(table[,summary.variable.label],
            by = list(cell_type = table[,cell.type.label]), 
            FUN = term)
}

summary_term_list <- function(table, group.id.variable, cell.type.variable,
                              summary.terms = c("variance", "mean", "max", "min")){
  lapply(summary.terms, function(term){
    unique.groups <- table[,group.id.variable] %>% unique()
    list.aggregate <- list(group = table[,group.id.variable])
    aggregate(table, by = list.aggregate, FUN = "summary_term_group")
    do.call(rbind, lapply(unique.groups, function(group){
      group.filter <- table[,group.id.variable]
      outputs.filtered <- table[group.filter,]
      list.argument <- list(cell_type = halo.outputs.table[,cell.type.label])
      summary_term_group(table = halo.outputs.table, term = term, cell.type.label = cell.type.label)
      colnames(summary.table)[2] <- term; summary.table$term <- term
      summary.table
    }))
  })
}

#--------
# scripts
#--------

# 01
sample.id.label <- levels.variable <- "Sample"
cell.area.variable <- "Nucleus_Area"
gene.marker.label <- "AKT3_Copies"
normalized.area.variable <- "log10_nucleus-area"
normalized.marker.variable <- "log10_akt3-copies"
#
output.updated.filename <- "halo-outputs_updated.Rdata"
output.updated.path <- here(save.path, output.updated.filename)
#
boxplot.area.jpg.name <- "ggboxplot_nucleus-area.jpg"
boxplot.area.jpg.path <- here(save.path, boxplot.area.jpg.name)
boxplot.marker.jpg.name <- "ggboxplot_marker-akt3-counts.jpg"
boxplot.marker.jpg.path <- here(save.path, boxplot.marker.jpg.name)

# 02
halo.marker.yaxis.label <- "AKT3 counts"
#
boxplot.marker.filename <-"ggboxplot_akt3-counts.jpg"
boxplot.log.marker.filename <-"ggboxplot-logscale_akt3-counts.jpg"
histogram.area.filename <- "histogram_nucleus-area.jpg"
histogram.marker.filename <- "histogram_akt3-marker.jpg"
#
boxplot.marker.path <- here(save.path, boxplot.marker.filename)
boxplot.log.marker.path <- here(save.path, boxplot.log.marker.filename)
histogram.area.path <- here(save.path, histogram.area.filename)
histogram.marker.path <- here(save.path, histogram.marker.filename)

# 03, summary statistics
boxplot.composite.area.summary.name <- "ggboxplot_slide-summary-statistics_norm-area.jpg"
boxplot.composite.area.summary.path <- here(save.path, boxplot.composite.area.summary.name)
anova.dependent.variable <- "Nucleus_Area"

# 04, shapiro test normality
shapiro.downsample.amount <- 5000
summary.terms = c("var", "mean", "max", "min")
cell.type.label <- "cell_type"

# 05, anova tests
anova.dependent.variable <- normalized.area.variable
model.results.list.path <- "list-results_linear-model_anova.rda"

anova.results.list.path <- here(save.path, "anova-results-list.rda")
explained.variance.title.string <- "Nucleus area (log10-transformed)"
# barplot with residuals
barplot.residuals.jpg.filepath <- "ggbar-perc-var_nalog10-complex_ro1-dlpfc.jpg"
barplot.residuals.jpg.path <- file.path(save.path, barplot.residuals.jpg.filepath)
# barplot with residuals, log10 scale
barplot.residuals.log10.jpg.filename <- "ggbar-perc-var-log10_nalog10-complex_ro1-dlpfc.jpg"
barplot.residuals.log10.jpg.path <- file.path(save.path, barplot.residuals.log10.jpg.filename)
# barplot, no residuals
barplot.noresiduals.jpg.filepath <- "ggbar-perc-expl-var_nalog10-complex_ro1-dlpfc.jpg"
barplot.noresiduals.jpg.path <- file.path(save.path, barplot.noresiduals.jpg.filepath)

# 99 quantile scale summaries
