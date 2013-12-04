#!/usr/bin/env Rscript





#######################
#      Main           #
#######################

args <- commandArgs(TRUE)
path_file = args[1]

d = read.table (path_file)
listType = c("Primary", "Secondary", "Tertiary", "Guanidium", "Imidazole", "AcidCarboxylic", "Global")


print (d)

png (paste(path_file, ".png", sep = ""), 6000, 6000)
par (mar = c(30,30,30,30))
barplot(d[,2],space=c(0.5,0.5,0.5), xlab="", ylab="",las=1, col= "#009DE0", ylim = c(0, 30), cex.axis = 10)

arrows(c(1,2.5,4,5.5,7,8.5,10,11.5,13,14.5),d[,2],c(1,2.5,4,5.5,7,8.5,10,11.5,13,14.5),d[,2]+d[,3],angle=90,length=1, lwd = 6)
arrows(c(1,2.5,4,5.5,7,8.5,10,11.5,13,14.5),d[,2],c(1,2.5,4,5.5,7,8.5,10,11.5,13,14.5),d[,2]-d[,3],angle=90,length=1, lwd = 6)

grid (NA, 6, col = "black", lwd = 6, lty = 1)

dev.off ()
