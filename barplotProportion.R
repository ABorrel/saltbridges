#!/usr/bin/env Rscript

source("tool.R")


#######################
#      Main           #
#######################

args <- commandArgs(TRUE)
pathData = args[1]
study_type = args[2]

d = read.table(pathData , sep = "\t", header = TRUE)

print (d)

for (li in 1:dim(d)[1]){
	tot = sum(d[li,])
	for (cl in 1:dim(d)[2]){
		d[li,cl] = d[li,cl] / tot
	}
}

print (d)

color = c("red", "orange", "blue", "grey")

png(filename=paste(pathData,".png",sep = ""),width=as.integer(600))
par(xpd=T, mar=par()$mar+c(2,0,0,6)) 
barplot(t(d), ylim = c(0,1), main = "", ylab="Frequencies", xlab = "", col=color, space=0.6, cex.axis=1, las=2, cex=1, axes = TRUE)
legend("right",legend=colnames(d), fill=color,inset=c(-0.2,0))
dev.off()

