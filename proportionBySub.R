#!/usr/bin/env Rscript
source("tool.R")



#######################
#      Main           #
#######################

args <- commandArgs(TRUE)
p_file = args[1]

d = read.table (p_file, sep = "\t", header = TRUE)
l_color = defColor (colnames (d))


for (sub in rownames (d)){

	# pie plot
	leg = NULL
	for (l in colnames (d)){
		leg = append (leg, paste (l, "\n", round (d[sub,l]/sum(d[sub,])*100), "%", sep = ""))
	}

	svg(filename=paste(p_file, "_", sub, "_pie.svg", sep = ""), 10, 10)
	try(pie(as.double(d[sub,]), col = l_color, label = leg, cex = 1.6))
	dev.off()
	
	# transfom %
	d_temp = GetPercent (d[sub,], 0)
	print (d_temp)

	d_temp = as.double (d_temp)
	
	names(d_temp) = colnames(d)
	
	svg(filename=paste(p_file, "_", sub, "_barplot.svg", sep = ""), 8, 10)
	par (mar = c(5,5,2,2))

	barplot (d_temp, col = l_color, cex.lab = 2, ylab = "Frequencies", ylim = c(0,100), cex.axis = 1.6)
	
	# grid
	# horizontal
	y_grid = seq(0, 100, 10)
	for (y in y_grid){
		segments (0, y, length(d_temp) + 1.2, y, lty = 2, col = "black", lwd = 1.5)
	}


	# vertical -> need improvment with a loop
	#x_grid = seq(1, length(d_temp) + 1 ,5)
	#for (x in x_grid){
	#	segments (x, 0, x, 1, lty = 2, col = "black", lwd = 1.5)
	#}

	dev.off()
}



