#!/usr/bin/env R

# Author: Sean Maden
#
# Functions for MAE summaries.
#


summaries_df_list <- function(mae, filter.type.label){
  # summaries_df_list
  #
  # mae MultiAssayExperiment
  # filter.type.label label for the filter used, if any
  #
  #
  #
  #
  
  dfmap <- mae@sampleMap
  dfmap$region <- gsub(".*_", "", dfmap$primary)
  dfmap$brnum <- gsub("_.*", "", dfmap$primary)
  
  variables.vector <- c("primary", "region", "brnum")
  
  # get summaries by data set type
  df.all <- do.call(cbind, 
                    lapply(variables.vector, 
                           function(variable.iter){
                             dfmap$variable.iter <- dfmap[,variable.iter]
                             df.iter <- dfmap %>% 
                               as.data.frame() %>%
                               distinct(assay, variable.iter) %>%
                               group_by(assay) %>%
                               summarize(n()) %>%
                               as.data.frame()
                             colnames(df.iter) <- c("assay", variable.iter)
                             df.iter
                           })) %>% as.data.frame()
  df.all <- df.all[,c(1,2,4,6)]
  df.train <- do.call(cbind, 
                      lapply(variables.vector, 
                             function(variable.iter){
                               dfmap$variable.iter <- dfmap[,variable.iter]
                               df.iter <- dfmap %>% 
                                 as.data.frame() %>%
                                 filter(primary %in% sample.id.train) %>%
                                 distinct(assay, variable.iter) %>%
                                 group_by(assay) %>%
                                 summarize(n()) %>%
                                 as.data.frame()
                               colnames(df.iter) <- c("assay", variable.iter)
                               df.iter
                             })) %>% as.data.frame()
  df.train <- df.train[,c(1,2,4,6)]
  df.validate <- do.call(cbind, 
                         lapply(variables.vector, 
                                function(variable.iter){
                                  dfmap$variable.iter <- dfmap[,variable.iter]
                                  df.iter <- dfmap %>% 
                                    as.data.frame() %>%
                                    filter(primary %in% sample.id.validate) %>%
                                    distinct(assay, variable.iter) %>%
                                    group_by(assay) %>%
                                    summarize(n()) %>%
                                    as.data.frame()
                                  colnames(df.iter) <- c("assay", variable.iter)
                                  df.iter
                                })) %>% as.data.frame()
  df.validate <- df.validate[,c(1,2,4,6)]
  
  # make wide data.frame
  df.all.wide <- df.all; df.train.wide <- df.train; df.validate.wide <- df.validate
  colnames(df.all.wide)[2:4] <- paste0(colnames(df.all.wide)[2:4],".all")
  colnames(df.train.wide)[2:4] <- paste0(colnames(df.train.wide)[2:4],".train")
  colnames(df.validate.wide)[2:4] <- paste0(colnames(df.validate.wide)[2:4],".validate")
  
  df.wide <- cbind(df.all.wide, 
                   cbind(df.train.wide[,c(2:4)], df.validate.wide[,c(2:4)])) %>% 
    as.data.frame()
  
  # make df tall tables
  df.all$data.type <- "all"
  df.train$data.type <- "train"
  df.validate$data.type <- "validate"
  df.tall <- rbind(df.all,
                   rbind(df.train, df.validate)) %>% as.data.frame()
  
  # append formatted assay name and platform name
  assays.keep <- c("snrnaseq.k2.all", "cell.sizes", "bulk.rnaseq")
  df.tall <- df.tall[df.tall[,1] %in% assays.keep,]
  df.wide <- df.wide[df.wide[,1] %in% assays.keep,]
  df.tall$platform.name <- df.tall$assay.type <- 
    df.wide$platform.name <- df.wide$assay.type <- "NA"
  
  assay.name.iter <- assays.keep[1]
  df.wide[df.wide$assay==assay.name.iter,]$assay.type <- "snRNAseq"
  df.wide[df.wide$assay==assay.name.iter,]$platform.name <- "10X Chromium"
  df.tall[df.tall$assay==assay.name.iter,]$assay.type <- "snRNAseq"
  df.tall[df.tall$assay==assay.name.iter,]$platform.name <- "10X Chromium"
  
  assay.name.iter <- assays.keep[2]
  df.wide[df.wide$assay==assay.name.iter,]$assay.type <- "fluorescent in situ hybridization"
  df.wide[df.wide$assay==assay.name.iter,]$platform.name <- "RNAscope"
  df.tall[df.tall$assay==assay.name.iter,]$assay.type <- "fluorescent in situ hybridization"
  df.tall[df.tall$assay==assay.name.iter,]$platform.name <- "RNAscope"
  
  assay.name.iter <- assays.keep[3]
  df.wide[df.wide$assay==assay.name.iter,]$assay.type <- "bulk RNAseq"
  df.wide[df.wide$assay==assay.name.iter,]$platform.name <- "Illumina HiSeq"
  df.tall[df.tall$assay==assay.name.iter,]$assay.type <- "bulk RNAseq"
  df.tall[df.tall$assay==assay.name.iter,]$platform.name <- "Illumina HiSeq"
  
  df.tall$filter.type <- filter.type.label
  df.wide$filter.type <- filter.type.label
  
  return(list(wide = df.wide, tall = df.tall))
}