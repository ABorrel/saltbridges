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
  v_middle = append (v_middle, MiddlePoint (d[i,5], d[i,4]))
}
print (v_middle)


svg (paste(p_file, "_CI.svg", sep = ""), 16, 10)
par (mar=c(4,5,1,2))
plot (d[,1], v_middle, pch = 8, ylim = c(min(d[,4]), max (d[,5])), cex = 3, cex.lab = 2, ylab = "pKa", xlab = "% of counter ion", cex.axis = 1.5)
text (d[,1] - 2, v_middle, labels = rownames(d), cex = 1.6)

print (rownames(d))

#error bar
arrows(d[,1], d[,4], d[,1], d[,5], angle=90, lwd = 2)
arrows(d[,1], d[,5], d[,1], d[,4], angle=90, lwd = 2)
dev.off()


svg (paste(p_file, "_water.svg", sep = ""), 16, 10)
par (mar=c(4,5,1,2))
plot (d[,2], v_middle, pch = 8, ylim = c(min(d[,4]), max (d[,5])), cex = 3, cex.lab = 2, ylab = "pKa", xlab = "% of water", cex.axis = 1.5)
text (d[,2] - 2, v_middle, labels = rownames(d), cex = 1.6)

#error bar
arrows(d[,2], d[,4], d[,2], d[,5], angle=90, lwd = 2)
arrows(d[,2], d[,5], d[,2], d[,4], angle=90, lwd = 2)
dev.off()


svg (paste(p_file, "_CIwater.svg", sep = ""), 16, 10)
par (mar=c(4,5,1,2))
plot (d[,3], v_middle, pch = 8, ylim = c(min(d[,4]), max (d[,5])), cex = 3, cex.lab = 2, ylab = "pKa", xlab = "% of counter ion and water", cex.axis = 1.5)
text (d[,3] - 2, v_middle, labels = rownames(d), cex = 1.6)

#error bar
arrows(d[,3], d[,4], d[,3], d[,5], angle=90, lwd = 2)
arrows(d[,3], d[,5], d[,3], d[,4], angle=90, lwd = 2)
dev.off()