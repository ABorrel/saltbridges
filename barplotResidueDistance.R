#!/usr/bin/env Rscript

source("tool.R")


calculmatrix = function(inMatrix){

	nbCol = dim(inMatrix)[2]
	nbLine = dim(inMatrix)[1]
	for (i in seq(1,nbLine)){
		for (j in seq(2,nbCol)){
			temp = j - 1
			while(temp >= 2){
				inMatrix[i,j] =  inMatrix[i,j] - inMatrix[i,temp]
				temp = temp -1 
			}
		}
	}
	return (inMatrix)
}


##################
#     MAIN       #
##################


args <- commandArgs(TRUE)
file = args[1]

d = read.table (file, sep = "\t", header = TRUE)


print (d)
l_c = defColor (rownames(d))
print (l_c)


svg (filename=paste(file,".svg", sep = ""), 20, 10)
barplot(as.matrix(d), ylab= "Number of amino acids", beside = TRUE, col = l_c, xlab = "Distance from query group (Å)", names.arg = colnames (d))
legend( "topleft", legend = rownames (d), bty = "n", fill = l_c)
dev.off()

svg (filename=paste(file,"_cum.svg", sep = ""), 20, 10)
barplot(as.matrix(d), ylab= "Number of amino acids", beside = FALSE, col = l_c, xlab = "Distance from query group (Å)")
legend( "topleft", legend = rownames (d), bty = "n", fill = l_c)
dev.off()