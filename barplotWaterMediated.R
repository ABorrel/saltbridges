#!/usr/bin/env Rscript
source("tool.R")



barplotCum = function (d_in, l_color, main_plot){
  # transfom %
  d_percent = GetPercent(d_in, 0)
  
  barplot (t(d_percent), col = l_color, axes = FALSE, axisnames = FALSE, main = main_plot, cex.main = 2)
  
  #axis 
  cum = 0
  for (y in d_percent){
    if (y > 3){
      axis (2, (cum + y / 2), paste (round(y), "%", sep = ""), las = 2, cex.axis = 2.9)
    }
    cum = cum + y
  }
}



#######################
#      Main           #
#######################


args <- commandArgs(TRUE)
p_filin = args[1]

d = read.table(p_filin, header = TRUE)
l_color = defColor(colnames (d))

svg (paste(p_filin, "_combine.svg", sep = ""), bg = "transparent" , 15, 9)
nf <- layout(matrix(c(0,0,0,0,0,0,1,2,3,4,5,6),2,6,byrow=TRUE), c(1,1,1,1,1,1), c(0,3), TRUE)
par (mar=c(1,6,2,3))


for (sub in rownames(d)){
  barplotCum (d[sub,], l_color, sub)

}



dev.off()