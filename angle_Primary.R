#!/usr/bin/env Rscript


source ("tool.R")






#######################
#        MAIN         #
#######################

args <- commandArgs(TRUE)
file = args[1]

data = read.csv(file, sep = "\t", header = FALSE)


nbLine = dim(data)[1]
nb_col = dim(data)[2]


l_element = unique(data[,nb_col])
color_element  = defColor(l_element)

for (element in l_element){
	
	data_plot = data[which(data[,nb_col] == element),]

	png(filename = paste(file , element, ".png" , sep = ""), width=600, height = 600)
	par(mar=c(6,6,6,6))
	plot(data[,1], data[,2], xlab = "Distance Å", ylab = "Angles en degres", type = "n", xlim = c(1.5,5), main = element)
	try(points(data_plot[,2]~data_plot[,1], col = color_element[element], pch = 19))
	dev.off()
}

brk = seq(0,180,20)
png(filename=paste(file, "_distanceDistribution.png", sep = ""))

hist(data[,2], xlab ="Distances (Å)", ylab = "Number of occurences", breaks=brk, right=F, main=paste("Distribution of angles","\n", "Primary amines", sep = ""),las=1, freq=T, col = "#D8D8D8")

dev.off()
