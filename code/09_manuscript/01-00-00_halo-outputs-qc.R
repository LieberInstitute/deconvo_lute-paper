#!/usr/bin/env R

#
#
#
#

libv <- c("ggplot2", "gridExtra")

#-------------
# manage paths
#-------------
# get save dpath
code.dname <- "09_manuscript"
proj.dname <- "deconvo_method-paper"
save.dpath <- file.path(proj.dname, "outputs", code.dname)

#---------------
# load halo data
#---------------
fname <- "halo_all.Rdata"
path <- file.path("Human_DLPFC_Deconvolution", "processed-data", "03_HALO", fname)
dfh <- get(load(path))

#--------------------
# get transformations
#--------------------
# nuclear area
# do log10 transformation
dfh$nuc.area.log10 <- log10(dfh$Nucleus_Area)

# akt3 
# do quantile transform by subject
samplev <- unique(dfh$Sample)
datv <- unlist(lapply(samplev, function(si){
  message("working on sample: ", si)
  aktv <- dfh[dfh$Sample==si, ]$AKT3_Copies
  qv <- quantile(aktv, seq(0, 1, 1e-3))
  qlabv <- as.numeric(gsub("%", "", names(qv)))
  qv2 <- qv[1:(length(qv)-1)]
  newv <- rep("NA", length(aktv))
  for(ii in seq(length(qv2))){
    which.aktv <- aktv >= qv[ii] & aktv < qv[ii+1]
    newv[which.aktv] <- qlabv[ii+1]
  }
  return(newv)
}))
datv <- as.numeric(datv)
datv[is.na(datv)] <- median(datv, na.rm = T) # replace NAs
plot(density(datv))

dfh$akt3.copies.qscale <- as.numeric(datv)


# qnorm method 2
# library(preprocessCore)

dfh$akt3.qnorm2 <- as.numeric(normalize.quantiles(as.matrix(dfh$AKT3_Copies, ncol = 1)))
quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

#-------------------
# plot distributions
#-------------------
dfp1 <- data.frame(untransformed = dfh$Nucleus_Area,
                  transformed = dfh$nuc.area.log10,
                  variable = "Nucleus_Area")
dfp2 <- data.frame(untransformed = dfh$AKT3_Copies,
                   transformed = dfh$akt3.copies.qscale,
                   variable = "AKT3_Copies")
dfp <- rbind(dfp1, dfp2)
for(i in seq(2)){dfp[,i] <- as.numeric(dfp[,i])}

p1 <- ggplot(dfp, aes(x = untransformed)) + 
  geom_density() + theme_bw() +
  ggtitle("Untransformed") + 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  facet_wrap(~variable)
p2 <- ggplot(dfp, aes(x = transformed)) + 
  geom_density() + theme_bw() +
  ggtitle("Transformed") + 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  facet_wrap(~variable)

fname <- "ggdens-composite_cellsize-halo_ro1-dlpfc.jpg"
jpeg(file.path(save.dpath, fname), width = 5, height = 3.5,
     units = "in", res = 400)
grid.arrange(p1, p2, nrow = 2, left = "Density", bottom = "Value")
dev.off()

#----------------------
# check key assumptions
#----------------------
# equality of variance between groups
varv.na <- tapply(dfh$Nucleus_Area, INDEX = dfh$SAMPLE_ID, FUN=var)
varv.akt3 <- tapply(dfh$AKT3_Copies, INDEX = dfh$SAMPLE_ID, FUN=var)

# normality
plot(hist(dfh$Nucleus_Area))
plot(hist(dfh$AKT3_Copies))
plot(hist(1/(max(dfh$AKT3_Copies+1) - dfh$AKT3_Copies)))
plot(hist(log(dfh$AKT3_Copies)))

datv <- dfh$AKT3_Copies
shapiro.test(sample(datv, 5000))

datv <- 1/(dfh$AKT3_Copies+1e-10)
shapiro.test(sample(datv, 5000))

datv <- sqrt(dfh$AKT3_Copies)
plot(hist(datv))
shapiro.test(sample(datv, 5000))

#-------------------------------------------------
# fit some models and check information criterions
#-------------------------------------------------
library(nlme)

dep.var <- "Nucleus_Area"
model.str <- paste0(dep.var, " ~ cell_type + Slide + SAMPLE_ID")

mod1 <- eval(parse(text = paste0("lm(", model.str, ", data = dfh)")))
mod2 <- eval(parse(text = paste0("gls(", model.str, ", data = dfh, method = 'REML')")))
mod3 <- eval(parse(text = paste0("glm(", model.str, ", data = dfh)")))

anova(mod)
anova(mod, lm)

#---------------
# anova analysis
#---------------
# nucleus area
dep.var <- "nuc.area.log10"
# basic model
model.str <- paste0(dep.var, " ~ cell_type + Slide + SAMPLE_ID")
avi1 <- eval(parse(text = paste0("aov(", model.str, ", data = dfh)")))
# complex model
model.str <- paste0(dep.var, " ~ cell_type + Slide + Combo + Position + BrNum")
avi2 <- eval(parse(text = paste0("aov(", model.str, ", data = dfh)")))

# model residuals
res <- residuals(avi)
qqnorm(res)
qqline(res)

# akt3
dep.var <- "AKT3_Copies"
# basic model
model.str <- paste0(dep.var, " ~ cell_type + Slide + Combo + Position + BrNum")
avi3 <- eval(parse(text = paste0("aov(", model.str, ", data = dfh)")))
# complex model
model.str <- paste0(dep.var, " ~ cell_type + Slide + Combo + Position + BrNum")
avi4 <- eval(parse(text = paste0("aov(", model.str, ", data = dfh)")))

#----------------------------
# explained variance analysis
#----------------------------
library(ggplot2)

title.str <- "Nucleus area (log10-transformed)"
avi <- avi2
dfp <- as.data.frame(summary(avi)[[1]])
dfp <- data.frame(variable = gsub(" ", "", rownames(dfp)), ssq = dfp[,2])
dfp$perc.var <- 100*dfp$ssq/sum(dfp$ssq)
dfp$label <- paste0(round(dfp$perc.var, 1), "%")

level.order <- rev(order(dfp$perc.var))
dfp$variable <- factor(dfp$variable, levels = dfp$variable[level.order])

ggbar <- ggplot(dfp, aes(x = variable, y = perc.var, fill = variable)) +
  geom_bar(stat = "identity") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Total percent variance") + ggtitle(title.str) +
  geom_text(aes(x = variable, y = perc.var, label = label)) +
  theme(legend.position = "none")
fname <- "ggbar-perc-var_nalog10-complex_ro1-dlpfc.jpg"
jpeg(file.path(save.dpath, fname), width = 3.5, height = 2.5, 
     units = "in", res = 400)
ggbar; dev.off()

ggbar <- ggplot(dfp, aes(x = variable, y = perc.var, fill = variable)) +
  geom_bar(stat = "identity") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  ylab("Percent variance\n(log10-scaled)") + 
  scale_y_log10() + ggtitle(title.str) +
  geom_text(aes(x = variable, y = perc.var, label = label))
fname <- "ggbar-perc-var-log10_nalog10-complex_ro1-dlpfc.jpg"
jpeg(file.path(save.dpath, fname), width = 3.5, height = 2.5, 
     units = "in", res = 400)
ggbar; dev.off()

dfpe <- dfp[!dfp$variable=="Residuals",]
dfpe$perc.var <- 100*dfpe$ssq/sum(dfpe$ssq)
dfpe$label <- paste0(round(dfpe$perc.var, 1), "%")
ggbar <- ggplot(dfpe, aes(x = variable, y = perc.var, fill = variable)) +
  geom_bar(stat = "identity") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  ylab("Total explained variance") + ggtitle(title.str) +
  geom_text(aes(x = variable, y = perc.var, label = label))
fname <- "ggbar-perc-expl-var_nalog10-complex_ro1-dlpfc.jpg"
jpeg(file.path(save.dpath, fname), width = 3.5, height = 2.5, 
     units = "in", res = 400)
ggbar; dev.off()


