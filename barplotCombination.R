#!/usr/bin/env Rscript

# By BORREL Alexandre
# 05-2013


############
#   MAIN   #
############

args <- commandArgs(TRUE)
path_file = args[1]


d = read.table(path_file, sep = "\t", header = FALSE)
print (d[seq (1,10),])
d = d[order(d[,2],decreasing = TRUE),]
print (d[seq (1,10),])

png (paste(path_file, ".png", sep = ""), 1820, 800)
par(mar = c(30,4,4,4))
barplot(d[seq (1,100),2], names.arg = d[seq (1,100),1], las = 2, cex.names= 1.5,cex.axis= 1.5, col = grey, main = "")
dev.off()





