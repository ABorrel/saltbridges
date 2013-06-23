#!/usr/bin/env Rscript



# MAIN
args <- commandArgs(TRUE)
filin = args[1]
l_plot = list()

data = read.table (filin, sep = "\t", header = FALSE)

#print (length (l_plot))
#print (l_plot)
png (paste (filin, ".png", sep = ""), 800, 800)
par(mar = c(8,4,4,4))
boxplot (data[,1]~data[,2], las = 2, cex.axis = 1.5)
dev.off()


