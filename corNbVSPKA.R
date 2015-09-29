#!/usr/bin/env Rscript

source("tool.R")


############
#   MAIN   #
############

args <- commandArgs(TRUE)
p_file = args[1]

d = read.table(p_file, header = TRUE)

print (d)

v_middle = NULL
for (i in seq (1,dim(d)[1])){
  v_middle = append (v_middle, MiddlePoint (d[i,3], d[i,4]))
}
print (v_middle)

svg (paste(p_file, ".svg", sep = ""),  bg = "transparent", 12, 10)
par (mar=c(4,5,1,2))
plot (v_middle, d[,1], pch = 8, ylim = c(min(d[,1] - d[,2]), max (d[,1] + d[,2])), xlim = c(min (d[,3]), max (d[,4])), cex = 3, cex.lab = 2, xlab = "pKa", ylab = "Number of neighbours", cex.axis = 1.5)
text (v_middle + 0.4,d [,1] + 1, labels = rownames(d), cex = 1.6)

#error bar
arrows(v_middle, d[,1] - d[,2],  v_middle, d[,1] + d[,2],  angle=90, lwd = 2)
arrows(v_middle, d[,1] + d[,2],v_middle, d[,1] - d[,2],  angle=90, lwd = 2)

arrows(d[,4], d[,1], d[,3], d[,1], angle=90, lwd = 2)
arrows(d[,3], d[,1], d[,4], d[,1], angle=90, lwd = 2)

dev.off()

