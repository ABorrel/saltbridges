#!/usr/bin/env Rscript

source("tool.R")




#######################
#      Main           #
#######################

args <- commandArgs(TRUE)
pathData = args[1]
pathData_distance = args[2]
study_type = args[3]

data = read.table(pathData , sep = "\t", header = TRUE)
data_distance = read.table(pathData_distance , sep = "\t")
data_distance = data_distance[,-1]

color = defColor (colnames (data))

png(filename=paste(pathData,".png",sep = ""),width=as.integer(600))
par(xpd=T, mar=par()$mar+c(0,0,0,10)) 
barplot(t(data), ylim = c(0,1), main=paste("Distribution neighbors ", study_type, sep = ""), ylab="Frequencies", xlab = "", col=color, space=0.6, cex.axis=1, las=1, cex=1, axes = TRUE)
legend("right",legend=colnames(data), fill=color,inset=c(-0.2,0))
dev.off()

png(filename=paste(pathData_distance ,".png",sep = ""),width=as.integer(400))
par(mfrow = c(3,1))

plot(density(as.double(data_distance[1,])), xlim = c(2,5))
plot(density(as.double(data_distance[2,])), xlim = c(2,5))
plot(density(as.double(data_distance[3,])), xlim = c(2,5))

dev.off()

