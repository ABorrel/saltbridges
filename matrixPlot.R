#!/usr/bin/env Rscript


require(graphics)

generateLegend = function (length_tot, cut){
  
  list_out = seq(0,length_tot, cut)
  if (list_out[length(list_out)] != length_tot){
    list_out = append(list_out, length_tot)
  }
  return (list_out)
}

generatePosition = function (list_legend, value_ecart){
  
  list_out = c()
  for (element in list_legend){
    list_out = append(list_out,element * value_ecart)
  }
  return (list_out)
}



addTrans <- function(color,trans)
{
  # This function adds transparancy to a color.
  # Define transparancy with an integer between 0 and 255
  # 0 being fully transparant and 255 being fully visable
  # Works with either color and trans a vector of equal length,
  # or one of the two of length 1.
  
  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
  
  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}





cardMatrix = function(matrixIN, name_file){
  
  nb_col = dim(matrixIN)[2]
  nb_line = dim(matrixIN)[1]
  
  dim_x = nb_line
  dim_y = nb_col
  
  if (nb_col == 1){
    return ()
  }
  if (nb_line == 1){
    return ()
  }
  
  if (nb_col < 30){
    dim_x = 30
  }
  
  if (nb_line < 30){
    dim_y = 30
  }
  
  bk = c(0,20,40,60,80,100) 
  
  svg (file = paste (name_file, ".svg", sep = ""), 15 + (dim_x * 0.4), 15 + (dim_y * 0.4))
  par( mar=c(15,15,0.5,0.5))
  image(as.matrix(matrixIN), yaxt = "n", xaxt = "n", breaks = bk, col = c("white", "cyan","deepskyblue", "red", "black"))
  grid(nx = nb_line, ny = nb_col, col = "black", lwd = 1, lty = 1)
  box()
  # place les petites barres 
  axis(1,seq(0,1,(1/(nb_line-1))), labels = FALSE)
  axis(2,seq(0,1,(1/(nb_col-1))), labels = FALSE)
  
  # place les positions en fonction du cut
  ecart1 = 1/(nb_line-1)	
  ecart2 = 1/(nb_col-1)
  list_L1 = generateLegend (nb_line,1)
  list_L2 = generateLegend (nb_col,1)
  
  # place les legendes
  posX = generatePosition(list_L1, ecart1)
  posY = generatePosition(list_L2, ecart2)
  axis(1,seq(0,1,(1/(nb_line-1))),rownames (matrixIN), cex.axis = 2.25, las = 2)
  axis(2,seq(0,1,(1/(nb_col-1))),rownames (matrixIN), cex.axis = 2.25, las = 2)
  
  #legend ("right", fill = c("black", "darkred", "red", "darkmagenta", "darkorchid","deepskyblue", "cyan", "white"), legend = c("0-2", "2-3", "3-4", "4-5", "5-6","6-7", "7-8", "> 10"), bg = addTrans ("#FFFFFF", 120) )
  dev.off()
}


###########
#  MAIN   #
###########

args = commandArgs(TRUE)
path_matrix = args[1]

d = read.table (path_matrix, header = T, sep = "\t")

cardMatrix (d, path_matrix)
