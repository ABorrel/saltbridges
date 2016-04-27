#!/usr/bin/env Rscript

source("tool.R")


#######################
#      Main           #
#######################

args <- commandArgs(TRUE)
pathData = args[1]
study_type = args[2]

data = read.table(pathData , sep = "\t", header = TRUE)
color = defColor (colnames (data))

svg (filename=paste(pathData,".svg",sep = ""),bg = "transparent", 15,10)
par(xpd=T, mar=par()$mar+c(0,2,0,10)) 
barplot(t(data), ylim = c(0,1), ylab="", xlab = "", col=color, space=0.6, cex.axis=3, las=1, cex=3, axes = TRUE)
legend("right",legend=colnames(data), fill=color,inset=c(-0.2,0), cex = 2)
dev.off()

