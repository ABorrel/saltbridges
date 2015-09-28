#!/usr/bin/env Rscript

source("tool.R")


############
#   MAIN   #
############

args <- commandArgs(TRUE)
p_interact = args[1]
p_nb_neighbor = args[2]
pr_result = args[3]

d_interact = read.table(p_interact, header = TRUE)
d_neighbor = read.table(p_nb_neighbor, header = TRUE)


d_percentage = GetPercent (d_interact, 0)
rownames(d_percentage) = rownames(d_interact)
colnames(d_percentage) = colnames(d_interact)
print (d_percentage)
print (d_neighbor[,2])

yplot = NULL
for (sub in rownames(d_percentage)){
  if (sub == "COO"){
    yplot = append (yplot, d_percentage[sub,"N"])
  }else{
    yplot = append (yplot, d_percentage[sub,"COO"])
  }
}

names(yplot) = rownames(d_percentage)


svg (paste(pr_result, "NbneighbourVSInteract.svg", sep = ""), 16, 10)
par (mar=c(4,5,1,2))


plot (d_neighbor[,1], yplot, pch = 8, xlim = c(min(d_neighbor[,1] - d_neighbor[,2]), max (d_neighbor[,1] + d_neighbor[,2])), ylim = c(0,100), cex = 3, cex.lab = 2, ylab = "% of counter ion", xlab = "Number of neighbours", cex.axis = 1.5)
text (d_neighbor[,1] + 1, yplot + 1, labels = names(yplot), cex = 1.6)

#error bar
arrows(d_neighbor[,1] - d_neighbor[,2], yplot, d_neighbor[,1] + d_neighbor[,2], yplot, angle=90, lwd = 2)
arrows(d_neighbor[,1] + d_neighbor[,2], yplot, d_neighbor[,1] - d_neighbor[,2], yplot, angle=90, lwd = 2)

#arrows(d[,3], d[,1], d[,4], d[,1], angle=90, lwd = 2)
#arrows(d[,4], d[,1], d[,3], d[,1], angle=90, lwd = 2)

dev.off()

