#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
file = args[1]

data = read.table(file, sep = "\t")

sumData = dim(data)[1]
brk = seq(0,max(data[,1]),1)


png(filename=paste(file,".png",sep = ""), 6000, 6000)
par (mar = c(30,30,30,30))
hist(data[,1], ylab ="", xlab = "", breaks=brk, xlim=c(0,max(data[,1])),  right=F, main="" , freq=T, col = "#009DE0", cex = 20, cex.lab = 20, cex.axis = 10)#, yaxt = "n", xaxt = "n" )
max_count = max (hist(data[,1], plot = FALSE)$count)

#horizontal
for (i in seq (0,max_count,10)){
	segments (0,i,max(data[,1]),i,lwd = 7)
}

# vertical
for (i in seq (0,max(data[,1]),5)){
	segments (i,0,i,max(data[,1]),lwd = 7)
}

dev.off()


