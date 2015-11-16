#!/usr/bin/env Rscript
source("tool.R")


barplotCum = function (d_in, l_color, main_plot){
	# transfom %
	d_percent = GetPercent(d_in, 0)

	barplot (t(d_percent), col = l_color, axes = FALSE, axisnames = FALSE, main = main_plot, cex.main = 3.2)
	
	#axis 
	cum = 0
	for (y in d_percent){
	  if (y > 3){
		  axis (2, (cum + y / 2), paste (round(y), "%", sep = ""), las = 2, cex.axis = 2.9)
	  }
		cum = cum + y
	}
}


testX2 = function (d1, d2){

	d_x2 = rbind (d1, d2)
	out_X2 = chisq.test(as.matrix(d_x2))

	return (out_X2$p.value)
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
	svg (paste(pr_out, sub, "_combine.svg", sep = ""), bg = "transparent", 15, 9)
	nf <- layout(matrix(c(0,0,0,1,2,3),2,3,byrow=TRUE), c(1,1,4), c(0,3), TRUE)
	par (mar=c(1,6,1,3))
	par (omar = c(0,0,0,0))

	# X2
	pval_atleast = testX2 (GetPercent(d_atleast1[sub,], 0), GetPercent(d_atleast2[sub,],0))

	barplotCum (d_atleast1[sub,], l_color, "")
	barplotCum (d_atleast2[sub,], l_color, "")
	
	nb_sub1 = sum (d_atleast1[sub,])
	nb_sub2 = sum (d_atleast2[sub,])

	d_percent1 = GetPercent (d_notatleast1[sub,], nb_sub1)
	d_percent2 = GetPercent (d_notatleast2[sub,], nb_sub2)
	
	pval_notatleast = testX2 (d_percent1, d_percent2)
	

	d_plot = matrix(c(as.double(d_percent1), as.double(d_percent2)), 2, length (d_percent1),byrow=TRUE )

	# barplot
	colnames (d_plot) = colnames (d_notatleast1)
	color_all = NULL
	for (colors in l_color){
		color_all = append (color_all, rep (colors,2)) 
	}
	barplot (d_plot, col = color_all, cex.lab = 3.2, ylab = "Frequencies", ylim = c(0,100), cex.names = 4, cex.axis = 2.4, beside=TRUE, cex.main = 3)
	text (10,96, paste("PDB1.5; n = ", nb_sub1, sep = ""), cex = 4)
	text (10,88, paste("PDB3.0; n = ", nb_sub2, sep = ""), cex = 4)

	# grid
	# horizontal
	y_grid = seq(0, 100, 10)
	for (y in y_grid){
		segments (0, y, length(d_percent1)*5 , y, lty = 2, col = "black", lwd = 1.5)
	}

}



