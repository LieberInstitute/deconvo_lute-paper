#!/usr/bin/env R

libv <- c("ggplot2", "reshape2")
sapply(libv, library, character.only = T)

# Author: Sean Maden
#
# Plot to show direction change in expression and prediction with scale factor change.
#

plotTitleString <- "Approximate affect of cell type scale change"
yLabelString <- "New - Old\n-           +"

#---------------------
# barplots with facets
#---------------------

# get plot data as differences, new-old
dfp <- data.frame(
  cellScaleFactor = c(1, 0.5, -0.5, -1),
  referenceExpression = c(-1, -0.5, 0.5, 1),
  bulkPrediction = c(1, 0.5, -0.5, -1)
)
dfp <- melt(dfp)
levelsVector <- c("increase", "slight_increase", "slight_decrease", "decrease")
dfp$type <- rep(levelsVector, 3)
dfp$type <- factor(dfp$type, levels = levelsVector)
changeLevelsVector <- c("Increase", "Decrease")
dfp$Change <- 
  ifelse(dfp$value > 0, changeLevelsVector[1], changeLevelsVector[2])
dfp$Change <- factor(dfp$Change, levels = changeLevelsVector)

#dfp <- t(dfp) |> as.data.frame()

ggplot(dfp, aes(x = variable, y = value, fill = Change)) + 
  geom_bar(stat="identity", color = "black") + theme_bw() +
  ylab(yLabelString) + facet_wrap(~type, ncol = 1) + 
  geom_hline(yintercept = 0) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45,hjust=1)) +
  ggtitle(plotTitleString) +
  scale_fill_manual(breaks = changeLevelsVector, 
                    values=c("dodgerblue", "gold"))

ggplot(dfp, aes(x = variable, y = value, fill = Change)) + 
  geom_bar(stat="identity", color = "black") + theme_bw() +
  ylab(yLabelString) + facet_wrap(~type, nrow = 1) + 
  geom_hline(yintercept = 0) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45,hjust=1),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  ggtitle(plotTitleString) +
  scale_fill_manual(breaks = changeLevelsVector, 
                    values=c("dodgerblue", "gold"))

#---------------------------------------
# line segments with arrowheads, faceted
#---------------------------------------
# get plot data as differences, new-old
dfp <- data.frame(
  cellScaleFactor = c(1, 0.5, -0.5, -1),
  referenceExpression = c(-1, -0.5, 0.5, 1),
  bulkPrediction = c(1, 0.5, -0.5, -1)
)
dfp <- melt(dfp)
levelsVector <- c("increase", "slight_increase", "slight_decrease", "decrease")
dfp$type <- rep(levelsVector, 3)
dfp$type <- factor(dfp$type, levels = levelsVector)
changeLevelsVector <- c("Increase", "Decrease")
dfp$Change <- 
  ifelse(dfp$value > 0, changeLevelsVector[1], changeLevelsVector[2])
dfp$Change <- factor(dfp$Change, levels = changeLevelsVector)

dfp$value.start <- 0

ggplot(dfp, aes(x = variable, y = value, fill = Change)) + 
  geom_line(aes(start = value.start, )) + theme_bw() +
  ylab(yLabelString) + facet_wrap(~type, nrow = 1) + 
  geom_hline(yintercept = 0) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45,hjust=1),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  ggtitle(plotTitleString) +
  scale_fill_manual(breaks = changeLevelsVector, 
                    values=c("dodgerblue", "gold"))

