#!/usr/bin/env Rscript
source("tool.R")


barplotCum = function (d_in, l_color, main_plot){
	# transfom %
	d_percent = GetPercent(d_in, 0)

	barplot (t(d_percent), col = l_color, axes = FALSE, axisnames = FALSE, main = main_plot, cex.main = 2)
	
	#axis 
	cum = 0
	for (y in d_percent){
		axis (2, (cum + y / 2), paste (round(y), "%", sep = ""), las = 2, cex.axis = 2.0)
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
	svg (paste(pr_out, sub, "_combine.svg", sep = ""), 15, 9)
	nf <- layout(matrix(c(0,0,0,1,2,3),2,3,byrow=TRUE), c(1,1,4), c(0,3), TRUE)
	par (mar=c(1,6,2,3))

	# X2
	pval_atleast = testX2 (GetPercent(d_atleast1[sub,], 0), GetPercent(d_atleast2[sub,],0))

	barplotCum (d_atleast1[sub,], l_color, "PDB1.5")
	barplotCum (d_atleast2[sub,], l_color, paste("PDB3.0 ", signifPvalue (pval_atleast), sep = ""))
	
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
	barplot (d_plot, col = color_all, cex.lab = 2, ylab = "Frequencies", ylim = c(0,100), cex.names = 2.0, cex.axis = 1.8, beside=TRUE, cex.main = 3)
	text (13,96, paste("PDB1.5; n = ", nb_sub1, sep = ""), cex = 3)
	text (13,92, paste("PDB3.0; n = ", nb_sub2, sep = ""), cex = 3)

	# grid
	# horizontal
	y_grid = seq(0, 100, 10)
	for (y in y_grid){
		segments (0, y, length(d_percent1)*5 , y, lty = 2, col = "black", lwd = 1.5)
	}

}



