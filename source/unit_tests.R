######################
# running lpb expt on large set of samples

library(here)
proj.dpath <- "deconvo_method-paper"
source.dpath <- file.path("deconvo_method-paper",
                          "source")
source.fnv <- c("pb_methods", 
                "z_transformations")
for(fni in source.fnv){
  source(
    file.path(here(), 
              source.dpath, 
              paste0(fni, ".R")))}




expt1 <- get_pb_experiment(datv = rep(1,4))
expt2 <- get_pb_experiment(datv = seq(1000, 4*100))



#####################

# plotting issues -- plot wrapper functions dont save plots or return working objects.

library(here)
proj.dpath <- "deconvo_method-paper"
source.dpath <- file.path("deconvo_method-paper",
                          "source")
source.fnv <- c("pb_methods", 
                "z_transformations")
for(fni in source.fnv){
  source(
    file.path(here(), 
              source.dpath, 
              paste0(fni, ".R")))}





get_facet <- function(df){
  plot <- ggplot(df, aes(x = pi_est, y = pi_diff)) + 
    facet_wrap(~cell_type)
  return(plot)
}
lfun <- function(){
  dftall <- get_exe_dftall()
  plot <- get_facet(dftall)
  return(list(plot = plot))
}


dftall <- get_exe_dftall()
plot <- get_facet(dftall) # works

lp <- lfun()
lp$plot # works



















