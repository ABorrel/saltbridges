#!/usr/bin/env Rscript

# By BORREL Alexandre
# 04-2012

library (MASS)
require(plotrix)

source("tool.R")


# distribution
histData = function (d, p_png){

	l_type_neighbor = unique(d[,2])
	l_color = defColor (unique(d[,2]))

	l = list ()
	i = 1
	for (type_neighbor in l_type_neighbor){
		print (d[which (d[,2] == type_neighbor),1])
		l[[i]] = d[which (d[,2] == type_neighbor),1]
		i = i + 1
	}

	print (l)
	print (l_color)
	png (p_png, 1000, 800)
	multhist(l, col = l_color, cex.names = 1.5, freq = TRUE, cex.axis = 1.5)
	legend("topright", col=l_color, legend = names (l_color), pch=rep (16,length (l_color)),  cex = 1.5)
	dev.off()
}





# main #
########
args <- commandArgs(TRUE)
p_file = args[1]

d = read.table (p_file, sep = "\t", header = FALSE)
histData(d,paste(p_file, ".png", sep = ""))

