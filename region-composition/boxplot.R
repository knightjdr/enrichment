#!/usr/local/bin/Rscript

library(cowplot)
library(ggplot2)

# plot longest disordered region
listTable = read.table('longest.txt', sep="\t", header = TRUE)
columns = c('Tier1', 'Tiers12', 'Tiers123', 'Tiers1234', 'Proteome')
orderedList = listTable[, columns]

pdf('longest_boxplot_noOutliers.pdf')
boxplot(x = as.list(orderedList), data = orderedList, outline = FALSE, xlab = "Longest disordered region", ylab = "Amino acids")
dev.off()
pdf('longest_boxplot.pdf')
boxplot(x = as.list(orderedList), data = orderedList, xlab = "Longest disordered region", ylab = "Amino acids")
dev.off()

# violin plots
orderedList.m <- reshape2::melt(orderedList, id.vars = NULL)
pdf('longest_violin.pdf')
ggplot(orderedList.m, aes(x = variable, y = value)) +
  geom_violin() +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  theme_cowplot() +
  xlab('Longest disordered region') +
  ylab('Amino acids')
dev.off()

# plot fraction disordered length
listTable = read.table('fraction.txt', sep="\t", header = TRUE)
orderedList = listTable[, columns]

pdf('fraction_boxplot_noOutliers.pdf')
boxplot(x = as.list(orderedList), data = orderedList, outline = FALSE, xlab = "Fraction disordered region", ylab = "Amino acids")
dev.off()
pdf('fraction_boxplot.pdf')
boxplot(x = as.list(orderedList), data = orderedList, xlab = "Fraction disordered region", ylab = "Amino acids")
dev.off()

# violin plots
orderedList.m <- reshape2::melt(orderedList, id.vars = NULL)
pdf('fraction_violin.pdf')
ggplot(orderedList.m, aes(x = variable, y = value)) +
  geom_violin() +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  theme_cowplot() +
  xlab('Fraction disordered region') +
  ylab('Amino acids')
dev.off()
