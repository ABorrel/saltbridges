#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
file = args[1]

data = read.table(file, sep = "\t")

sumData = dim(data)[1]
nameX = data[,2]
brk = seq(0,5,0.10)


png(filename=paste(file,".png",sep = ""), 6000, 6000)
par (mar = c(30,30,30,30))

#hist(data[,1], xlab ="Oxygen (COO-) of nitrogen (Nref) distances (Å)", ylab = "Number of occurences", breaks=brk, xlim=c(0,5), right=F, main="Distribution of distance amine counter-Ion",las=1, freq=T)

hist(data[,1], ylab ="", xlab = "", breaks=brk, xlim=c(0,5),  right=F, main="" , freq=T, col = "#009DE0", cex = 20, cex.lab = 20, cex.axis = 10)#, yaxt = "n", xaxt = "n" )
grid (5, 5, col = "black", lwd = 6, lty = 1)
dev.off()

#png(filename=paste(file,"density.png",sep = ""), 6000, 6000)
#par (mar = c(30,30,30,30), oma = c (5,5,5,5))

#plot (density(data[,1]))

#dev.off()


