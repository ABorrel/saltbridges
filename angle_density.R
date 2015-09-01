#!/usr/bin/env Rscript


source("tool.R")


revector<-function(n){
	out = NULL
	nb_value = length (n)
	for(i in nb_value:1){
		out = append (out,n[i])
	}
	return(out)
}


mapFig = function (main_plot, data, filin){
	nb_col = dim (data)[2]
	l_element = unique(data[,nb_col])
	color_element  = defColor(l_element)

	data_bis = data

	# color	
	Lab.palette = colorRampPalette(c("#FFFFFF", color_element), space = "Lab")

	svg(filename = paste(file, "_", main_plot, "_density.svg" , sep = ""), width = 10, height = 10)
	par (mar = c(6,6,2,2))
	smoothScatter(data_bis, colramp = Lab.palette , main = "", xlim =c(1,max(data_bis[,1])), cex = 2, cex.lab = 2, cex.axis = 1.5, xlab = "Distance (Å)", ylab = "Angle (degres)")
	for (x in seq (1, max (data_bis[,1]), 0.5)){
		segments (x, min(data_bis[,2]),x, max (data_bis[,2]), lty = 2, col = "black", lwd = 1.5 )
	}
	for (y in seq (0,max(data_bis[,2]),10)){
		segments (1, y, max (data_bis[,1]),y,lty = 2, col = "black", lwd = 1.5)
	}
	dev.off()


	svg(filename = paste(file, "_", main_plot, "_point.svg" , sep = ""), width=10, height = 10)
	par (mar = c(6,6,2,2))

	plot(data_bis[,1], data_bis[,2], xlim = c(1,max(data_bis[,1])), main = "", cex = 2, cex.lab = 2, cex.axis = 1.5, xlab = "Distance (Å)", ylab = "Angle (degres)", col = color_element[element], pch = 19)
	
	for (x in seq (1, max (data_bis[,1]), 0.5)){
		segments (x, min(data_bis[,2]),x, max (data_bis[,2]), lty = 2, col = "black", lwd = 1.5 )
	}
	for (y in seq (0,max(data_bis[,2]),10)){
		segments (1, y, max (data_bis[,1]),y,lty = 2, col = "black", lwd = 1.5)
	}
	dev.off()

}



#######################
#        MAIN         #
#######################

args <- commandArgs(TRUE)
file = args[1]

data = read.csv(file, sep = "\t", header = FALSE)


nbLine = dim(data)[1]
nb_col = dim(data)[2]


l_element = unique(data[,nb_col])

for (element in l_element){
	
	data_plot = data[which(data[,nb_col] == element),]
	try(mapFig (element, data_plot, file))

	if (dim (data)[2] == 4){
		png(filename = paste(file ,"_",  element, "_angle.png" , sep = ""), width=1000, height = 1000)
		par(mar=c(5,5,5,5))
		Lab.palette = colorRampPalette(revector(heat.colors(50)), space = "Lab")
		smoothScatter(data_plot[,c(2,3)], colramp = Lab.palette, main = element, cex.axis = 2.5, cex.main = 3, xlab = "angle1", ylab = "angle2", cex.lab = 3)
		dev.off()
	}
}

