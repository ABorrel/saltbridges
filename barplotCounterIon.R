#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
file = args[1]
distance = args[2]

data = read.table(file, header = FALSE,sep = "\t")
nameX = c("Primary", "Secondary", "Tertiary", "Imidazole")

nbline = dim(data)[1]
nbcol = dim(data)[2]


for (i in 1:nbline){
	sumLine = sum(data[i,])
	for(j in 1:nbcol){
		data[i,j]=data[i,j]/sumLine
	}
}


sumRetrieve = sum(data[,1]) + sum(data[,2])

png(filename=paste(file,".png",sep = ""))

barplot(t(data),main="Proportion Nref with a least one counter-ions ", ylim = c(0,1), ylab="Frequency", xlab = "Groups", col=heat.colors(2), space=0.5, cex.axis=0.8, las=1, names.arg=nameX, cex=0.8, axes = TRUE)

legend(3.7,0.30,legend=c("Carbonyl of E or D " , "Others"), cex=0.8, fill=heat.colors(2))

dev.off()

warnings()
