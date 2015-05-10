#!/usr/bin/env Rscript
source("tool.R")


mapFig = function (d, name, p_file){

	l_combi = unique (d[,3])

	l_col = defColorGrep (l_combi)

	print (p_file)

	png(filename = paste(p_file, name , "density.png", sep = "_"), width=1000, height = 1000)
	par(mar=c(5,5,5,5))

	smoothScatter(d, colramp = colorRampPalette (c("white", l_col[1])), main = "", xlim = c(1.5,5),  xlab = "", ylab = "", cex.axis = 1.5, cex.main = 2)
		
	for (x in seq (1.5, 5, 0.5)){
		segments (x,0,x,max (d[,2]),lwd = 2)
	}
	for (y in seq (0,max(d[,2]),10)){
		segments (1.5,y,5,y,lwd = 2)
	}

	dev.off()


	png(filename = paste(p_file,name, "point.png" , sep = "_"), width=1000, height = 1000)
	par(mar=c(5,5,5,5))

	plot(d[,1], d[,2], xlab = "Distance Å", ylab = "Angles",xlim = c(1.5,5), main = "", cex.axis = 1.5, cex.main = 2, col = l_col, pch = 19)

	for (x in seq (1.5, 5, 0.5)){
		segments (x,0,x,max (d[,2]),lwd = 2)
	}
	for (y in seq (0,max(d[,2]),10)){
		segments (1.5,y,5,y,lwd = 2)
	}


	dev.off()

}



#######################
#        MAIN         #
#######################

args <- commandArgs(TRUE)
p_file = args[1]


d = read.table (p_file, sep = "\t")


nb_line = dim(d)[1]
nb_col = dim(d)[2]
l_combi = unique(d[,nb_col])

print (l_combi)

for (combi in l_combi){
	d_plot = d[which(d[,nb_col] == combi),]
	mapFig (d_plot, combi, p_file)
}


