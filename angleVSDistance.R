#!/usr/bin/env Rscript

source("tool.R")



plotCombined = function (d_in, name_filout){
	
	# color
	nb_col = dim (d_in)[2]
	l_element = unique(d_in[,nb_col])
	color_element  = defColor(l_element)

	# svg
	svg (paste(name_filout, ".svg", sep = ""), 22, 22)
	top = max(c(d_in[,1], d_in[,2]))
	nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)
	par (mar=c(8,8,0,0))
	
	plot(d_in[,1], d_in[,2], pch=16, cex = 4, col = addTrans(color_element, 120), xlab = "Distances (Å)", ylab = "Angles ()", cex.lab = 3, cex.axis = 2.8, xlim = c(1,max(d_in[1])), ylim = c(0,180))
	
	# grid
	# horizontal
	y_grid = seq(0, 180, 10)
	for (y in y_grid){
		segments (1, y, max(d_in[,1]), y, lty = 2, col = "black", lwd = 1.5)
	}


	# vertical -> need improvment with a loop
	x_grid = seq(1, max(d_in[,1]) ,0.25)
	for (x in x_grid){
		segments (x, 0, x, 180, lty = 2, col = "black", lwd = 1.5)
	}

	
	par (mar=c(1,8,1,1))
	v = seq(1, max(d_in[,1]),0.25)
	barplot (hist(d_in[which(d_in[,1] >= 1 & d_in[,1] <= v[length (v)]), 1], plot = FALSE, breaks = v)$counts, axes = FALSE, col = addTrans(color_element, 120), main = "")
	axis (2, cex.axis = 3)	
	
	par (mar=c(3,1,1,1))
	v = seq(0, 180, 10)
	barplot (hist(d_in[,2], plot = FALSE, breaks = v)$counts, axes = FALSE, col = addTrans(color_element, 120) ,horiz = TRUE, main = "")
	axis (3, cex.axis = 3)

	dev.off ()


}






############
#   MAIN   #
############


args <- commandArgs(TRUE)
file = args[1]



d = read.csv(file, sep = "\t", header = FALSE)

nbLine = dim(d)[1]
nb_col = dim(d)[2]


l_element = unique(d[,nb_col])

for (element in l_element){
	
	data_plot = d[which(d[,nb_col] == element),]
	plotCombined (data_plot, paste(file, element, "combined", sep = ""))
}


