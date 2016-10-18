#!/usr/bin/env Rscript
source("tool.R")



barplotCum = function (d_in, l_color, main_plot){
  # transfom %
  d_percent = GetPercent(d_in, 0)
  
  barplot (t(d_percent), col = l_color, axes = FALSE, axisnames = FALSE, main = main_plot, cex.main = 3.2)
  
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

d_count = read.table (p_filin, header = TRUE)
i_null = which(apply (d_count, 2, sum) == 0)

# remove colum with only 0
if (is.integer0 (i_null) == FALSE){
	d_count = d_count[,-i_null]
}


l_sub = rownames (d_count)


#AFC (d_count, p_filin)

# barplot by proportion res
l_color = defColor (colnames (d_count))

for (sub in l_sub){
  d_count[sub,] = GetPercent (d_count[sub,], 0)
}

# wirte matrix percentage
write.csv(x = d_count, file = paste(p_filin,"percent.csv", sep = ""))


svg (paste(p_filin, ".svg", sep = ""), bg = "transparent" , 15, 9)
nf <- layout(matrix(c(0,0,0,0,0,0,0,1,2,3,4,5,6,7),2,7,byrow=TRUE), c(1,1,1,1,1,1,1), c(0,3), TRUE)
par (mar=c(1,6,2,3))


for (sub in rownames(d_count)){
  if (sub == "global"){
    barplotCum (d_count[sub,], l_color, "Any atoms")
  }else{
    barplotCum (d_count[sub,], l_color, sub)
  }
}

dev.off()