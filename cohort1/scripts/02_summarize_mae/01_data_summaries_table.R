#!/usr/bin/env R

# Author: Sean Maden
#
# Make data summaries table for manuscript.
#
#
#
#
#

library(dplyr)

# dataset names
path.data.pre.filter <- "./outputs/01_mae/mae_allsamples_append.rda"
path.data.post.filter <- "./outputs/01_mae/mae_analysis_append.rda"

list.sample.id <- get(load("./outputs/00_preprocess/list_snrnaseq_sampleid.rda"))
sample.id.train <- list.sample.id[["train"]]
sample.id.validate <- list.sample.id[["validation"]]









#----------------------------
# PRE FILTER DATA SETS
#----------------------------
# load
mae <- get(load(path.data.pre.filter))

# sample source summaries
# PLATFORM_TYPE X SAMPLE_SOURCE_ID_COUNT
df.samples.map.all <- mae@sampleMap %>% 
  as.data.frame() %>%
  distinct(assay, primary) %>%
  group_by(assay) %>%
  summarize(n())

df.samples.map.train <- mae@sampleMap %>% 
  as.data.frame() %>%
  filter(primary %in% sample.id.train) %>%
  distinct(assay, primary) %>%
  group_by(assay) %>%
  summarize(n())

df.samples.map.validate <- mae@sampleMap %>% 
  as.data.frame() %>%
  filter(primary %in% sample.id.validate) %>%
  distinct(assay, primary) %>%
  group_by(assay) %>%
  summarize(n())

df.samples.map.all$data.type <- "all"
df.samples.map.train$data.type <- "train"
df.samples.map.validate$data.type <- "validate"

new.table <- rbind(df.samples.map.all,
                   rbind(
                     df.samples.map.train,
                     df.samples.map.validate
                   ))
colnames(new.table) <- c("platform", "unique_sample_source_count", "dataset_type")
new.table <- new.table[,c(3,1,2)]

new.table$platform.name <- 
  new.table$assay.type <- "NA"
platform.names.iter <- c(
  "snrnaseq.k2.all", "snrnaseq.k3.all", "snrnaseq.k4.all",
  "bulk.pb.k2", "bulk.pb.k3", "bulk.pb.k4"
)
new.table[new.table$platform %in% platform.names.iter,]$assay.type <- "snRNAseq"
new.table[new.table$platform %in% platform.names.iter,]$platform.name <- "10X Chromium"
platform.names.iter <- c(
  "bulk.rnaseq", "bulk.rpkm.rnaseq"
)
new.table[new.table$platform %in% platform.names.iter,]$assay.type <- "bulk RNAseq"
new.table[new.table$platform %in% platform.names.iter,]$platform.name <- "Illumina HiSeq"
platform.names.iter <- c(
  "cell.sizes", "sce.img"
)
new.table[new.table$platform %in% platform.names.iter,]$assay.type <- "fluorescent in situ hybridization"
new.table[new.table$platform %in% platform.names.iter,]$platform.name <- "RNAscope"

# save
save(new.table, file = "./outputs/01_mae/table2_platforms_prefilter.rda")
write.csv(new.table, file = "./outputs/01_mae/table2_platforms_prefilter.csv", row.names = FALSE)

new.table.prefilter <- new.table














#----------------------------
# POST FILTER DATA SETS
#----------------------------
# load
mae <- get(load(path.data.post.filter))


# sample source summaries
# PLATFORM_TYPE X SAMPLE_SOURCE_ID_COUNT
df.samples.map.all <- mae@sampleMap %>% 
  as.data.frame() %>%
  distinct(assay, primary) %>%
  group_by(assay) %>%
  summarize(n())

df.samples.map.train <- mae@sampleMap %>% 
  as.data.frame() %>%
  filter(primary %in% sample.id.train) %>%
  distinct(assay, primary) %>%
  group_by(assay) %>%
  summarize(n())

df.samples.map.validate <- mae@sampleMap %>% 
  as.data.frame() %>%
  filter(primary %in% sample.id.validate) %>%
  distinct(assay, primary) %>%
  group_by(assay) %>%
  summarize(n())

df.samples.map.all$data.type <- "all"
df.samples.map.train$data.type <- "train"
df.samples.map.validate$data.type <- "validate"

new.table <- rbind(df.samples.map.all,
      rbind(
        df.samples.map.train,
        df.samples.map.validate
      ))
colnames(new.table) <- c("platform", "unique_sample_source_count", "dataset_type")
new.table <- new.table[,c(3,1,2)]

new.table$platform.name <- 
  new.table$assay.type <- "NA"
platform.names.iter <- c(
  "snrnaseq.k2.all", "snrnaseq.k3.all", "snrnaseq.k4.all",
  "bulk.pb.k2", "bulk.pb.k3", "bulk.pb.k4"
)
new.table[new.table$platform %in% platform.names.iter,]$assay.type <- "snRNAseq"
new.table[new.table$platform %in% platform.names.iter,]$platform.name <- "10X Chromium"
platform.names.iter <- c(
  "bulk.rnaseq", "bulk.rpkm.rnaseq"
)
new.table[new.table$platform %in% platform.names.iter,]$assay.type <- "bulk RNAseq"
new.table[new.table$platform %in% platform.names.iter,]$platform.name <- "Illumina HiSeq"
platform.names.iter <- c(
  "cell.sizes", "sce.img"
)
new.table[new.table$platform %in% platform.names.iter,]$assay.type <- "fluorescent in situ hybridization"
new.table[new.table$platform %in% platform.names.iter,]$platform.name <- "RNAscope"

new.table.postfilter <- new.table

# save
save(new.table, file = "./outputs/01_mae/table2_platforms.rda")
write.csv(new.table, file = "./outputs/01_mae/table2_platforms.csv", row.names = FALSE)










#---------------------
# save the environment
#---------------------

save.image(file = './env/02_summarize_mae/05_data_summaries_script.RData')