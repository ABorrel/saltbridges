#!/usr/bin/env Rscript




#######################
#        MAIN         #
#######################

args <- commandArgs(TRUE)
file = args[1]
data = read.csv(file, sep = "\t", header = FALSE)
print (data)


png(filename=paste(file, ".png", sep = ""))
par (mfrow = c(2,2))
hist(data[,1], xlab ="Angles between neighbor 1-2", ylab = "Number of occurences", right=F, main=paste("Angles", sep = ""),las=1, freq=T, col = "#D8D8D8")
hist(data[,2], xlab ="Angles between neighbor 1-3", ylab = "Number of occurences", right=F, main=paste("Angles", sep = ""),las=1, freq=T, col = "#D8D8D8")
hist(data[,3], xlab ="Angles between neighbor 2-3", ylab = "Number of occurences", right=F, main=paste("Angles", sep = ""),las=1, freq=T, col = "#D8D8D8")

hist(data[,3] + data[,1] + data[,2], xlab ="Sum angles", ylab = "Number of occurences", right=F, main=paste("Angles", sep = ""),las=1, freq=T, col = "#D8D8D8")


dev.off()
