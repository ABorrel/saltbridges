#!/usr/bin/env Rscript




#######################
#        MAIN         #
#######################

args <- commandArgs(TRUE)
file = args[1]

d = read.csv(file, sep = "\t", header = FALSE)

d = na.omit (d)

png(filename = paste(file, "_hist.png" , sep = ""), width=800, height = 2800)
par (mfrow = c (4,1))
hist (d[,1], main = "angle1VS2", col = "blue")
hist (d[,2], main = "angle1VS3", col = "blue")
hist (d[,3], main = "angle2VS3", col = "blue")
hist (d[,1] + d[,2] + d[,3], main = "Sum", col = "red")

dev.off()


png(filename = paste(file, "_cor.png" , sep = ""), width=800, height = 2800)
par (mfrow = c (4,1))
plot (d[,1], d[,2], main = "", pch = 19, col = "blue", xlim = c(0,180), ylim = c(0,180))
par (new = TRUE)
plot (d[,2], d[,3], main = "", pch = 19, col = "red", xlim = c(0,180), ylim = c(0,180))
par (new = TRUE)
plot (d[,1], d[,3], main = "", pch = 19, col = "green", xlim = c(0,180), ylim = c(0,180))


plot (d[,1], d[,2], main = "angle1VS2", pch = 19)
plot (d[,1], d[,3], main = "angle1VS3", pch = 19)
plot (d[,2], d[,3], main = "angle2VS3", pch = 19)

dev.off()
