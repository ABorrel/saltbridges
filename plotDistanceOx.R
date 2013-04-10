#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
file = args[1]

data = read.table(file, sep = "\t")

sumData = dim(data)[1]
nameX = data[,2]
brk = seq(0,5,0.05)


png(filename=paste(file,".png",sep = ""))

hist(data[,1], xlab ="Oxygen (COO-) of nitrogen (Nref) distances (Å)", ylab = "Number of occurences", breaks=brk, xlim=c(0,5), right=F, main="Distribution of distance amine counter-Ion",las=1, freq=T)

dev.off()
