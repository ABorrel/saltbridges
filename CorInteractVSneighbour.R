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

print (d_interact)
print (d_neighbor)


d_percentage = GetPercent (d_interact, 0)
rownames(d_percentage) = rownames(d_interact)
colnames(d_percentage) = colnames(d_interact)

yplot = NULL
for (sub in rownames(d_neighbor)){
  if (sub == "COO"){
    yplot = append (yplot, d_percentage[sub,"N"])
  }else{
    yplot = append (yplot, d_percentage[sub,"COO"])
  }
}

names(yplot) = rownames(d_percentage)


svg (paste(pr_result, "NbVSCI.svg", sep = ""), bg = "transparent", 12, 10)
par (mar=c(4,5,1,2))
plot ( yplot, d_neighbor[,1], pch = 8, ylim = c(min(d_neighbor[,1] - d_neighbor[,2]), max (d_neighbor[,1] + d_neighbor[,2])), xlim = c(0,100), cex = 3, cex.lab = 2, xlab = "% of charged groups", ylab = "Number of neighbours", cex.axis = 1.5)
text (yplot + 1,d_neighbor[,1] + 1, labels = rownames(d_neighbor), cex = 1.6)

#error bar
arrows(yplot, d_neighbor[,1] - d_neighbor[,2], yplot, d_neighbor[,1] + d_neighbor[,2], angle=90, lwd = 2)
arrows(yplot, d_neighbor[,1] + d_neighbor[,2], yplot, d_neighbor[,1] - d_neighbor[,2], angle=90, lwd = 2)

#arrows(d[,3], d[,1], d[,4], d[,1], angle=90, lwd = 2)
#arrows(d[,4], d[,1], d[,3], d[,1], angle=90, lwd = 2)

dev.off()


yplot = NULL
for (sub in rownames(d_neighbor)){
  if (sub == "COO"){
    yplot = append (yplot, d_percentage[sub,"N"] + d_percentage[sub,"HOH"] + d_percentage[sub,"NH"])
  }else{
    yplot = append (yplot, d_percentage[sub,"COO"] + d_percentage[sub,"HOH"] + d_percentage[sub,"OH"])
  }
}

names(yplot) = rownames(d_percentage)

svg (paste(pr_result, "NbVSCIALL.svg", sep = ""), bg = "transparent", 12, 10)
par (mar=c(4,5,1,2))
plot ( yplot, d_neighbor[,1], pch = 8, ylim = c(min(d_neighbor[,1] - d_neighbor[,2]), max (d_neighbor[,1] + d_neighbor[,2])), xlim = c(0,100), cex = 3, cex.lab = 2, xlab = "% of all ionizable groups", ylab = "Number of neighbours", cex.axis = 1.5)
text (yplot + 1,d_neighbor[,1] + 1, labels = rownames(d_neighbor), cex = 1.6)

#error bar
arrows(yplot, d_neighbor[,1] - d_neighbor[,2], yplot, d_neighbor[,1] + d_neighbor[,2], angle=90, lwd = 2)
arrows(yplot, d_neighbor[,1] + d_neighbor[,2], yplot, d_neighbor[,1] - d_neighbor[,2], angle=90, lwd = 2)

#arrows(d[,3], d[,1], d[,4], d[,1], angle=90, lwd = 2)
#arrows(d[,4], d[,1], d[,3], d[,1], angle=90, lwd = 2)

dev.off()

