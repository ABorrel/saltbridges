#!/usr/bin/env Rscript

source("tool.R")



plotCombined = function (d_in, name_filout){
	
	# color
	nb_col = dim (d_in)[2]
	l_element = unique(d_in[,nb_col])
	color_element  = defColor(l_element)

	# svg
	svg (paste(name_filout, ".svg", sep = ""), 13, 18)
	nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(2,1), c(1,3), TRUE)
	par (mar=c(3,8,1,1))
	
	plot(d_in[,1], d_in[,2], xaxt = "n", yaxt = "n", frame = TRUE, pch=16, cex = 4, col = addTrans(color_element, 120), xlab = "", ylab = "",xlim = c(2,max(d_in[1])), ylim = c(0,180))
	
	axis (1, cex.axis = 2.2, lty = 0)
	axis (2, cex.axis = 2.2, lty = 0)
	# append the number of subs
	text (2.4,165, label = paste("n = ", dim (d_in)[1], sep = ""), cex = 3)

	# grid
	# horizontal
	y_grid = seq(0, 180, 10)
	for (y in y_grid){
		segments (2, y, max(d_in[,1]), y, lty = 2, col = "black", lwd = 1.5)
	}


	# vertical -> need improvment with a loop
	x_grid = seq(2, max(d_in[,1]) ,0.25)
	for (x in x_grid){
		segments (x, 0, x, 180, lty = 2, col = "black", lwd = 1.5)
	}

	
	par (mar=c(1,8,1,1))
	v = seq(2, max(d_in[,1]),0.25)
	barplot (hist(d_in[which(d_in[,1] >= 2 & d_in[,1] <= v[length (v)]), 1], plot = FALSE, breaks = v)$counts, axes = FALSE, col = addTrans(color_element, 120), main = "")
	axis (2, cex.axis = 2.0)	
	
	par (mar=c(3,1,1,1))
	v = seq(0, 180, 10)
	barplot (hist(d_in[,2], plot = FALSE, breaks = v)$counts, axes = FALSE, col = addTrans(color_element, 120) ,horiz = TRUE, main = "")
	axis (3, cex.axis = 2.0)

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

# merge N basic
l_temp = c("Nim", "NaI", "Ngu")
data_plot_merged = NULL

for (element in l_element){
	if (element == "Nim" | element == "NaI" | element == "Ngu"){
	  data_plot_merged = rbind (data_plot_merged, d[which(d[,nb_col] == element),])  
	}
	data_plot = d[which(d[,nb_col] == element),]
	plotCombined (data_plot, paste(file, element, "combined", sep = "_"))
}

plotCombined (data_plot_merged , paste(file, "NimNaINgu_combined", sep = "_"))


