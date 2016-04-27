#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
file = args[1]

data = read.table(file, sep = "\t")

sumData = dim(data)[1]
#brk = seq(0,max(data[,1]),1)


svg(filename=paste(file,".svg",sep = ""), 6, 10, bg = "transparent")
par (mar = c(5,5,5,5))
hist(data[,1], ylab ="Number of occurrences", xlab = "Number of neighbors", xlim=c(0,max(data[,1])), main=paste ("Number of Neighbors ",round(mean (data[,1]),2), " +/- ", round(sd (data[,1]),2), sep = "") , freq=T, col = "grey", cex = 2, cex.lab = 1.5, cex.axis = 1.5, breaks = 7)
max_count = max (hist(data[,1], plot = FALSE)$count)

#horizontal
#for (i in seq (0,max_count,10)){
#	segments (0,i,max(data[,1]),i,lwd = 2, lty = 2)
#}

# vertical
#for (i in seq (0,max(data[,1]),5)){
#	segments (i,0,i,max_count,lwd = 2, lty = 2)
#}
#
dev.off()


