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

png(filename=paste(pathData,".png",sep = ""),width=as.integer(600))
par(xpd=T, mar=par()$mar+c(0,0,0,10)) 
barplot(t(data), ylim = c(0,1), main=paste("Distribution neighbors ", study_type, sep = ""), ylab="Frequencies", xlab = "", col=color, space=0.6, cex.axis=1, las=1, cex=1, axes = TRUE)
legend("right",legend=colnames(data), fill=color,inset=c(-0.2,0))
dev.off()

