#!/usr/bin/env Rscript


source("tool.R")




#######################
#        MAIN         #
#######################

args <- commandArgs(TRUE)
file = args[1]
data = read.csv(file, sep = "\t", header = FALSE)


nbLine = dim(data)[1]
nb_col = dim(data)[2]

# case with several angle -> deviation with median position
print (nb_col)
if (nb_col > 3){
	data = deviationAngle(data)
}

nbLine = dim(data)[1]
nb_col = dim(data)[2]

# retrieve type of neigbors
l_element = unique(data[,nb_col])
color_element  = defColor(l_element)

for (element in l_element){
	data_plot = data[which(data[,nb_col] == element),]
	png(filename = paste(file , element, ".png" , sep = ""), width=600, height = 600)
	par(mar=c(6,6,6,6))
	plot(data[,1], data[,2], xlab = "Distance Å", ylab = "Deviation of two angles", type = "n", xlim = c(1.5,5), main = as.character(element))
	try(points(as.double(data_plot[,2])~as.double(data_plot[,1]), col = color_element[element], pch = 19))
	dev.off()
	
	png(filename = paste(file , element, "density_deviation.png" , sep = ""), , width=6000, height = 6000)
	
	par(mar=c(30,30,30,30))

	smoothScatter(data_plot, main = "", xlim = c(1.5,4), cex.lab = 10, xlab = "", ylab = "", cex.axis = 10, cex.main = 6)
	dev.off()

}

png(filename=paste(file, "_angleDistribution.png", sep = ""))
hist(as.double(data[,2]), xlab ="Angles or deviation", ylab = "Number of occurences", right=F, main=paste("Distribution of deviation or angles", sep = ""),las=1, freq=T, col = "#D8D8D8")
dev.off()
