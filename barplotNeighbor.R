#!/usr/bin/env Rscript
source("tool.R")




#######################
#      Main           #
#######################

args <- commandArgs(TRUE)
pathData = args[1]
study_type = args[2]


png(filename=paste(pathData ,".png",sep = ""), width=as.integer(800), heigh = 1000)
par(mfrow = c(3,1))


for (i in 1:3){
	d = read.table(paste(pathData, "_",i, sep = ""), sep = "\t", header = TRUE)
	print (d)
	d = d[,-1]	

	color = defColor (colnames (d))
	barplot(t(d), ylab= "Number of occurences", col = color, xlab = "distance (Å)", names.arg = rownames (d),space = 0.5, main =i)

	if (i == 2){
		legend("left", legend=colnames(d), fill=color )
	}

}

dev.off()
