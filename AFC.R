#!/usr/bin/env Rscript


# By BORREL Alexandre
# 11-2013


library (FactoMineR)
source("tool.R")


factorAFC = function (xplot, yplot, xplotdata, yplotdata ){
	#return (1)
	factor = 1
	while ( abs(max(xplot)) < max(xplotdata) && abs(max(yplot)) < max(yplotdata) ){
		factor = factor + 0.5
		xplot = xplot * factor
		yplot = yplot* factor

	}
	return (factor)	
}



AFC = function (d, path_file){

	r = CA (d, graph = FALSE)

	svg (file = paste (path_file, "_AFC.svg", sep = ""), 15, 15)
	par(mar=c(8,8,8,8))

	# descriptors
	plot (r$row$coord[,1], r$row$coord[,2], type = "n", xlab = paste("DIM 1 : ", round(r$eig[1,2],1), "%", sep = ""), ylab = paste("DIM 2 : ", round(r$eig[2,2],1), "%", sep = ""), cex.lab = 2.4, ylim = c(min(r$row$coord[,2]) - 0.1, max(r$row$coord[,2])))
	col_des = defColor(names(r$col$coord[,1]))

	factor = factorAFC (r$col$coord[,1], r$col$coord[,2], r$row$coord[,1], r$row$coord[,2] )

	arrows (0,0,r$col$coord[,1]*factor, r$col$coord[,2]*factor, col = as.character(col_des), lwd = 3 )
	#text (r$col$coord[,1]*factor, r$col$coord[,2]*factor, labels = names(r$col$coord[,1]), col = as.character(col_des), cex = 1.4)

	# data
	text (r$row$coord[,1], r$row$coord[,2], col = "#009DE0", label = names (r$row$coord[,1]), cex = 2)
	abline(h=0,v=0)

	dev.off()
}




