#!/usr/bin/env Rscript
source("tool.R")

#######################
#      Main           #
#######################

args <- commandArgs(TRUE)
p_filin = args[1]

d_count = read.table (p_filin, header = TRUE)
i_null = which(apply (d_count, 2, sum) == 0)

# remove colum with only 0
if (is.integer0 (i_null) == FALSE){
	d_count = d_count[,-i_null]
}


l_sub = rownames (d_count)


AFC (d_count, p_filin)

# barplot by proportion res
l_color = defColor (colnames (d_count))

for (sub in l_sub){
  d_count[sub,] = GetPercent (d_count[sub,], 0)
}

svg(paste(p_filin, ".svg", sep = ""), 16, 10)
par (mar=c(4,5,1,2))
barplot (as.matrix(t(d_count)), col = l_color, cex.lab = 2, ylab = "Frequencies", ylim = c(0,1), cex.names = 1.7, cex.axis = 1.8, beside=FALSE, cex.main = 2.5)
  
# grid
# horizontal
y_grid = seq(0, 1, 0.1)
for (y in y_grid){
  segments (0, y, length(d_count)*5 , y, lty = 2, col = "black", lwd = 1.5)
}

legend ("topright",legend = names(l_color), fill = l_color )
dev.off()
  