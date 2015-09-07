#!/usr/bin/env Rscript


source ("tool.R")

args <- commandArgs(TRUE)
file = args[1]

d = read.table(file, sep = "\t")

sumData = dim(data)[1]
#brk = seq(0,5,0.10)

color_plot = addTrans(defColor (file), 120)

print (color_plot)

svg(filename=paste(file,".svg",sep = ""), 12, 10)
par (mar = c(6,6,2,2))
hist(d[,1], breaks = 20, right=F, main="" , freq=T, col = color_plot, cex = 2, cex.lab = 2, cex.axis = 1.5, xlab = "Angle (degree)", ylab = "Number of occurencies" )
dev.off()

