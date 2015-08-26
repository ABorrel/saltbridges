#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
file = args[1]

data = read.table(file, sep = "\t")

sumData = dim(data)[1]
#brk = seq(0,max(data[,1]),1)


png(filename=paste(file,".png",sep = ""), 800, 800)
par (mar = c(5,5,5,5))
hist(data[,1], ylab ="Quantity", xlab = "Number of neighbors", xlim=c(0,max(data[,1])),  right=F, main=paste (round(mean (data[,1]),2), " +/- ", round(sd (data[,1]),2), sep = "") , freq=T, col = "#009DE0", cex = 1, cex.lab = 1, cex.axis = 1)#, yaxt = "n", xaxt = "n" )
max_count = max (hist(data[,1], plot = FALSE)$count)

#horizontal
for (i in seq (0,max_count,10)){
	segments (0,i,max(data[,1]),i,lwd = 2, lty = 2)
}

# vertical
for (i in seq (0,max(data[,1]),5)){
	segments (i,0,i,max_count,lwd = 2, lty = 2)
}

dev.off()


