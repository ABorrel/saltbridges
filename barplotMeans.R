#!/usr/bin/env Rscript
source("tool.R")




#######################
#      Main           #
#######################

args <- commandArgs(TRUE)
p_file = args[1]

d = read.table (p_file, header = TRUE)
d_sub_type = d[,seq(3,dim(d)[2], 2)]
d_M = d[,1]
d_SD = d[,2]
names (d_M) = rownames(d)
print (d_M)

svg (paste(p_file, ".svg", sep = ""), bg = "transparent", 16, 10)
par (mar=c(4,5,1,2))
barplot(d_M, xlab="", ylab="", col= "grey", ylim = c(0, max(d_M + d_SD)), cex.axis = 2, cex = 1.6, space = 0.5)

#error bar
arrows(seq(1, by = 1.5, length.out = length(d_SD)),d_M - d_SD, seq(1, by = 1.5, length.out = length(d_SD)) ,d_M + d_SD, angle=90, length = 0.5, lwd = 2)
arrows(seq(1, by = 1.5, length.out = length(d_SD)),d_M + d_SD, seq(1, by = 1.5, length.out = length(d_SD)) ,d_M - d_SD, angle=90, length = 0.5, lwd = 2)

# grid
# horizontal
y_grid = seq(0, max(d_M + d_SD), 5)
if (length (y_grid <= 5)){
  y_grid = seq(0, max(d_M + d_SD), 2)
}
for (y in y_grid){
  segments (0, y, length(d_M) + 0.5*length(d_M), y, lty = 2, col = "black", lwd = 1.5)
}

dev.off ()

# to do with atom type
print (d_sub_type)
l_color = defColor(colnames (d_sub_type))

svg (paste(p_file, "_type.svg", sep = ""), 16, 10)
par (mar=c(4,5,1,2))

barplot(t(d_sub_type), xlab="", ylab="", col= l_color, ylim = c(0, max(d_M + d_SD)), cex.axis = 2, cex = 1.6, space = 0.5)

#error bar
arrows(seq(1, by = 1.5, length.out = length(d_SD)),d_M - d_SD, seq(1, by = 1.5, length.out = length(d_SD)) ,d_M + d_SD, angle=90, length = 0.5, lwd = 2)
arrows(seq(1, by = 1.5, length.out = length(d_SD)),d_M + d_SD, seq(1, by = 1.5, length.out = length(d_SD)) ,d_M - d_SD, angle=90, length = 0.5, lwd = 2)

# grid
# horizontal
y_grid = seq(0, max(d_M + d_SD), 5)
if (length (y_grid <= 5)){
  y_grid = seq(0, max(d_M + d_SD), 2)
}
for (y in y_grid){
  segments (0, y, length(d_M) + 0.5*length(d_M), y, lty = 2, col = "black", lwd = 1.5)
}

dev.off ()
