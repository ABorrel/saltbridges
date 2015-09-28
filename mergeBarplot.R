#!/usr/bin/env Rscript
source("tool.R")



#######################
#      Main           #
#######################


args <- commandArgs(TRUE)
p_atleast = args[1]
p_noatleast = args[2]
pr_out = args[3]

d_atleast = read.table(p_atleast, header = TRUE)
d_notatleast = read.table (p_noatleast, header = TRUE)

print (d_notatleast)
print (d_atleast)

# svg

l_color = defColor(colnames (d_atleast))
for (sub in rownames(d_atleast)){
	svg (paste(pr_out, sub, "_combine.svg", sep = ""), 13, 11)
	nf <- layout(matrix(c(0,0,1,2),2,2,byrow=TRUE), c(1,3), c(0,3), TRUE)
	par (mar=c(3,8,1,3))

	# transfom %
	nb_sub = sum (d_atleast[sub,])
	for (j in seq (1, dim(d_atleast)[2])){
		d_atleast[sub, j] = d_atleast[sub, j] / nb_sub
	}


	barplot (t(d_atleast[sub,]), col = l_color, axes = FALSE, axisnames = FALSE)
	
	#axis 
	cum = 0
	for (y in d_atleast[sub,]){
		axis (2, (cum + y / 2), paste (round(y), "%", sep = ""), las = 2, cex.axis = 2.0)
		cum = cum + y
	}


	# transfom %
	nb_sub = sum (d_notatleast[sub,])
	
	for (j in seq (1, dim(d_notatleast)[2])){
		d_notatleast[sub, j] = d_notatleast[sub, j] / nb_sub
	}

	d_temp = as.double(d_notatleast[sub,])
	names(d_temp) = colnames(d_notatleast)

	par (mar = c(5,5,2,2))

	barplot (d_temp, col = l_color, cex.lab = 2, ylab = "Frequencies", ylim = c(0,1), cex.names = 2.0, cex.axis = 1.8)
	
	# grid
	# horizontal
	y_grid = seq(0, 1, 0.1)
	for (y in y_grid){
		segments (0, y, length(d_temp) + 1.2, y, lty = 2, col = "black", lwd = 1.5)
	}

}



