#!/usr/bin/env Rscript


########################
#         MAIN         #
########################


library(lattice)
library(scatterplot3d)


require(methods)

args <- commandArgs(TRUE)
file = args[1]



data = read.csv(file, sep = "\t", header = FALSE)

png(paste(file,"_angle3D.png", sep = ""))
cloud(data[,2] ~ data[,3] * data[,4], screen = list(z = 40, x = 650, y=20), xlab = "angle1", ylab = "angle2", zlab = "angle3")
dev.off()
