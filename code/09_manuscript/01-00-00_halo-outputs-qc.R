#!/usr/bin/env R

#
#
#
#

# get save dpath
code.dname <- "09_manuscript"
proj.dname <- "deconvo_method-paper"
save.dpath <- file.path(proj.dname, "outputs", code.dname)

# load halo data
fname <- "halo_all.Rdata"
path <- file.path("Human_DLPFC_Deconvolution", "processed-data", "03_HALO", fname)
dfh <- get(load(path))

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
dep.var <- "Nucleus_Area"
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

dfp <- as.data.frame(summary(avi)[[1]])
dfp <- data.frame(variable = gsub(" ", "", rownames(dfp)), ssq = dfp[,2])
dfp$perc.var <- 100*dfp$ssq/sum(dfp$ssq)

level.order <- rev(order(dfp$perc.var))
dfp$variable <- factor(dfp$variable, levels = dfp$variable[level.order])

ggplot(dfp, aes(x = variable, y = perc.var, fill = variable)) +
  geom_bar(stat = "identity") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Total percent variance")

dfpe <- dfp[!dfp$variable=="Residuals",]
dfpe$perc.var <- 100*dfpe$ssq/sum(dfpe$ssq)
dfpe$label <- paste0(round(dfpe$perc.var, 2), "%")
ggplot(dfpe, aes(x = variable, y = perc.var, fill = variable)) +
  geom_bar(stat = "identity") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  ylab("Total explained variance") +
  geom_text(aes(x = variable, y = perc.var, label = label))

ggplot(dfpe, aes(x = variable, y = perc.var)) +
  geom_bar(stat = "identity", position = "stack") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  ylab("Total explained variance")


