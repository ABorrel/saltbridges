#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
file = args[1]

d = read.table(file, sep = "\t")

sumData = dim(data)[1]
#brk = seq(0,5,0.10)


png(filename=paste(file,".png",sep = ""), 600, 600)
par (mar = c(6,6,2,2))
hist(d[,1], breaks = 20, xlim=c(0,max(d[,1])),  right=F, main="" , freq=T, col = "#009DE0", cex = 2, cex.lab = 2, cex.axis = 1, xlab = "Distance (Ang.)", ylab = "Number of occurencies" )
#max_count = max (hist(d[,1], plot = FALSE)$count)
#vertical
#segments (0,0,0,max_count,lwd = 6)
#segments (1,0,1,max_count,lwd = 6)
#segments (2,0,2,max_count,lwd = 6)
#segments (3,0,3,max_count,lwd = 6)
#segments (4,0,4,max_count,lwd = 6)
#segments (5,0,5,max_count,lwd = 6)

#horizontal
#print (max_count)
#for (i in seq (0,max_count,100)){
#        print (i)
#	segments (0,i,5,i,lwd = 6)
#}



dev.off()

#png(filename=paste(file,"density.png",sep = ""), 6000, 6000)
#par (mar = c(30,30,30,30), oma = c (5,5,5,5))

#plot (density(data[,1]))

#dev.off()


