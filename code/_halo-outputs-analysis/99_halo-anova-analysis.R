source("deconvo_method-paper/code/_halo-outputs-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
halo.outputs.table <- get(load(halo.output.path))
halo.outputs.table <- halo.outputs.table %>% as.data.frame()

# check key assumptions
# equality of variance between groups
cell.area.vector <- halo.outputs.table[,cell.area.variable]
gene.marker.vector <- halo.outputs.table[,gene.marker.label]
sample.id.vector <- halo.outputs.table[,sample.id.label]
variance.nucleus.area <- tapply(cell.area.vector, INDEX = sample.id.vector, FUN = var)
variance.gene.marker <- tapply(gene.marker.vector, INDEX = sample.id.vector, FUN = var)

# check normality
# plot histograms

halo.outputs.table[,cell.area.variable] %>% hist() %>% plot()

halo.outputs.table[,gene.marker.label] %>% hist() %>% plot()

vector <- max(halo.outputs.table[,gene.marker.label] + 1) - halo.outputs.table[,gene.marker.label]

transformed.halo.outputs <- 1/vector

halo.outputs.table[,gene.marker.label] %>% log() %>% hist() %>% plot()

# shapiro-wilk test of normality

shapiro.test.transform1 <- gene.marker.vector %>% sample(5000) %>% shapiro.test()

inverse.transform <- 1/(halo.outputs.table[,gene.marker.label]+1e-10)
shapiro.test.transform2 <- inverse.transform %>% sample(5000) %>% shapiro.test()

sqrt.transform <- sqrt(gene.marker.vector)
shapiro.test.transform2 <- sqrt.transform %>% sample(5000) %>% shapiro.test()

# get anova models

dependent.variable <- "Nucleus_Area"

model.string <- paste0(dependent.variable, " ~ cell_type + Slide + SAMPLE_ID")

model1 <- paste0("lm(", model.string, ", data = halo.outputs.table)") %>% 
  parse() %>% eval()

model2 <- paste0("gls(", model.string, ", data = halo.outputs.table, method = 'REML')") %>% 
  parse() %>% eval()

model3 <- paste0("glm(", model.string, ", data = halo.outputs.table)") %>% 
  parse() %>% eval()

# anova analysis

dependent.variable <- cell.area.log.variable

# basic model

model.string <- paste0(dependent.variable, " ~ cell_type + Slide + SAMPLE_ID")

anova.string <- paste0("aov(", model.string, ", data = halo.outputs.table)")

anova.model1 <- anova.string %>% parse() %>% eval()

# complex model
model.string <- paste0(dependent.variable, " ~ cell_type + Slide + Combo + Position + BrNum")
anova.string <- paste0("aov(", model.string, ", data = halo.outputs.table)")
anova.model2 <- anova.string %>% parse() %>% eval()
# akt3
dependent.variable <- gene.marker.label
model.string <- paste0(dependent.variable, " ~ cell_type + Slide + Combo + Position + BrNum")
anova.model3 <- eval(parse(text = paste0("aov(", model.string, ", data = halo.outputs.table)")))
# complex model
model.string <- paste0(dependent.variable, " ~ cell_type + Slide + Combo + Position + BrNum")
anova.string <- paste0("aov(", model.string, ", data = halo.outputs.table)")
anova.model4 <- anova.string %>% parse() %>% eval()

# make residual plots
# get residuals
residuals1 <- residuals(anova.model1)
residuals2 <- residuals(anova.model2)
residuals3 <- residuals(anova.model3)
# do quantile-quantile plots
# do quantile-quantile plots
qqnorm(residuals1)
qqnorm(residuals2)
qqnorm(residuals3)
# do quantile-quantile lines
qqline(residuals1)
qqline(residuals2)
qqline(residuals3)

# explained variance analysis
explained.variance.title.string <- "Nucleus area (log10-transformed)"
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
filename <- "ggbar-perc-var_nalog10-complex_ro1-dlpfc.jpg"
path <- file.path(save.directory.path, filename)
jpeg(path, width = 3.5, height = 2.5, units = "in", res = 400); ggbar; dev.off()

ggbar <- ggplot(dfp, aes(x = variable, y = perc.var, fill = variable)) +
  geom_bar(stat = "identity") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  ylab("Percent variance\n(log10-scaled)") + 
  scale_y_log10() + ggtitle(title.str) +
  geom_text(aes(x = variable, y = perc.var, label = label))
filename <- "ggbar-perc-var-log10_nalog10-complex_ro1-dlpfc.jpg"
path <- file.path(save.directory.path, filename)
jpeg(path, width = 3.5, height = 2.5, units = "in", res = 400); ggbar; dev.off()

dfpe <- dfp[!dfp$variable=="Residuals",]
dfpe$perc.var <- 100*dfpe$ssq/sum(dfpe$ssq)
dfpe$label <- paste0(round(dfpe$perc.var, 1), "%")
ggbar <- ggplot(dfpe, aes(x = variable, y = perc.var, fill = variable)) +
  geom_bar(stat = "identity") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  ylab("Total explained variance") + ggtitle(title.str) +
  geom_text(aes(x = variable, y = perc.var, label = label))
filename <- "ggbar-perc-expl-var_nalog10-complex_ro1-dlpfc.jpg"
path <- file.path(save.directory.path, filename)
jpeg(path, width = 3.5, height = 2.5, units = "in", res = 400); ggbar; dev.off()
