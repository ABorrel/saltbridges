#!/usr/bin/env Rscript
source("tool.R")


GetPercent = function (d_in){

	nb_sub = sum (d_in)
	for (j in seq (1, dim(d_in)[2])){
		d_in[sub, j] = d_in[sub, j] / nb_sub
	}

	return (d_in)
}



barplotCum = function (d_in, l_color){
	# transfom %
	d_percent = GetPercent(d_in)

	barplot (t(d_percent), col = l_color, axes = FALSE, axisnames = FALSE)
	
	#axis 
	cum = 0
	for (y in d_percent){
		axis (2, (cum + y / 2), paste (round(y*100), "%", sep = ""), las = 2, cex.axis = 2.0)
		cum = cum + y
	}
}

#######################
#      Main           #
#######################


args <- commandArgs(TRUE)
p_atleast1 = args[1]
p_noatleast1 = args[2]
p_atleast2 = args[3]
p_noatleast2 = args[4]
pr_out = args[5]

d_atleast1 = read.table(p_atleast1, header = TRUE)
d_notatleast1 = read.table (p_noatleast1, header = TRUE)

d_atleast2 = read.table(p_atleast2, header = TRUE)
d_notatleast2 = read.table (p_noatleast2, header = TRUE)

# svg

l_color = defColor(colnames (d_atleast1))
for (sub in rownames(d_atleast1)){
	svg (paste(pr_out, sub, "_combine.svg", sep = ""), 15, 8)
	nf <- layout(matrix(c(0,0,0,1,2,3),2,3,byrow=TRUE), c(1,1,4), c(0,3), TRUE)
	par (mar=c(2,6,1,3))

	barplotCum (d_atleast1[sub,], l_color)
	barplotCum (d_atleast2[sub,], l_color)
	
	d_percent1 = GetPercent (d_notatleast1[sub,])
	d_percent2 = GetPercent (d_notatleast2[sub,])
	

	d_plot = matrix(c(as.double(d_percent1), as.double(d_percent2)), 2, length (d_percent1),byrow=TRUE )
	#d_plot = t(d_plot)

	#colnames (d_plot) = 

	color_all = NULL
	for (colors in l_color){
		color_all = append (color_all, rep (colors,2)) 
	}

	barplot (d_plot, col = color_all, cex.lab = 2, ylab = "Frequencies", ylim = c(0,1), cex.names = 2.0, cex.axis = 1.8, beside=TRUE)
	
	# grid
	# horizontal
	y_grid = seq(0, 1, 0.1)
	for (y in y_grid){
		segments (0, y, length(d_percent1)*5 , y, lty = 2, col = "black", lwd = 1.5)
	}

}



