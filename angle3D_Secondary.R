#!/usr/bin/env Rscript


source("tool.R")



classAngle = function(data){
	nbCol = dim(data)[2]
	nbLine = dim(data)[1]
	for (i in seq(1,nbLine)){
		if (data[i,2] < data[i,3]){
			key = data[i,2]
			data[i,2] = data[i,3]
			data[i,3] = key
		}
	}
	return (data)
}




#######################
#        MAIN         #
#######################

library(lattice)
library(scatterplot3d)
require(methods)

args <- commandArgs(TRUE)
file = args[1]

data = read.csv(file, sep = "\t", header = FALSE)
data = classAngle(data)
png(paste(file,"_angle3D.png", sep = ""))
cloud(data[,1] ~ data[,2] * data[,3], screen = list(z = 40, x = 650, y=20), xlab = "distance Å", ylab = "angle2", zlab = "angle3")
dev.off()


