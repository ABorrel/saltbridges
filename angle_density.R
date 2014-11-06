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
	#if (dim(data)[2] > 3 ){
	#	nb_line = dim (data)[1]
	#	nb_col = dim (data)[2]
	#	for (i in seq (1,nb_line)){
	#		a = c()
	#		for (j in 2:nb_col-1){
	#			a = append(a,data[i,j])
	#		}
	#	data_bis = rbind (data_bis, c(data[i,1], mean (a)))
	#	}
	#}
	#else {
	#	data_bis = data
	#}
	
	png(filename = paste(file, "_", main_plot, "_density.png" , sep = ""), width=1000, height = 1000)
	par(mar=c(5,5,5,5))

	smoothScatter(data_bis, col = color_element[main_plot], main = main_plot, xlim = c(1.5,5), xlab = "", ylab = "", cex.axis = 1.5, cex.main = 2)
		
	for (x in seq (1.5, 5, 0.5)){
		segments (x,0,x,max (data_bis[,2]),lwd = 2)
	}
	for (y in seq (0,max(data_bis[,2]),10)){
		segments (1.5,y,5,y,lwd = 2)
	}

	dev.off()


	png(filename = paste(file, "_", main_plot, "_point.png" , sep = ""), width=1000, height = 1000)
	par(mar=c(5,5,5,5))

	plot(data_bis[,1], data_bis[,2], xlab = "Distance Å", ylab = "Angles",xlim = c(1.5,5), main = as.character(element), cex.axis = 1.5, cex.main = 2, col = color_element[element], pch = 19)

	for (x in seq (1.5, 5, 0.5)){
		segments (x,0,x,max (data_bis[,2]),lwd = 2)
	}
	for (y in seq (0,max(data_bis[,2]),10)){
		segments (1.5,y,5,y,lwd = 2)
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

