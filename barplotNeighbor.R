#!/usr/bin/env Rscript
source("tool.R")




#######################
#      Main           #
#######################

args <- commandArgs(TRUE)
pathData = args[1]
nb_type = as.integer(args[2])


png(filename=paste(pathData ,".png",sep = ""), width=as.integer(800), heigh = 1000)
par(mfrow = c(nb_type/2+1,2))


for (i in 1:nb_type){
	d = read.table(paste(pathData, "_",i, sep = ""), sep = "\t", header = TRUE)
	
	indiv =  (sum(d))
	color = defColor (colnames (d))
	barplot(t(d), ylab= "Number of occurences", col = color, xlab = "distance (Å)", names.arg = rownames (d),space = 0.5, main =paste(i, indiv, sep = "   "))

	if (i == 2){
		legend("left", legend=colnames(d), fill=color)

	}

}

dev.off()
