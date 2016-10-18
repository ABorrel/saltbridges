#!/usr/bin/env Rscript


is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}


############
# PIE plot #
############

pieType = function (d, path_out){
	
	#print (d)
	colors = defColor(names(d))


	par (lwd = 1000)
	png(filename=paste(path_out,".png",sep = ""),400, 400)
	try(pie(as.double(d), col = colors, label = names(d), lwd = 10))
	dev.off()

	svg(filename=paste(path_out,".svg",sep = ""))
	try(pie(as.double(d), col = colors, label = names(d)))
	dev.off()
}


############
# AFC-plot #
############

library (FactoMineR)

factorAFC = function (xplot, yplot, xplotdata, yplotdata ){
	# max window
	print (yplotdata)

	max_window_x = max (abs (xplotdata))
	max_window_y = max (abs (yplotdata))

	# max coordinate
	max_x = max (abs (xplot))
	max_y = max (abs (yplot))


	factor = 1
	while ( max_x < max_window_x  && max_y < max_window_y ){
		factor = factor + 0.1
		max_x = max_x * factor
		max_y = max_y * factor

	}


	if (factor == 1 ){
		while ( max_x > max_window_x || max_y > max_window_y ){

			#print (max_x)
			#print (max_y)
			#print (max_window_x)
			#print (max_window_y)
			#print ("********")

			factor = factor - 0.2
			max_x = max_x * factor
			max_y = max_y * factor
			
		}
	}
	return (factor)	
}



AFC = function (d, path_file){

  # genere afc coords
	r = CA (d, graph = FALSE)
  
	# plot
	xplot = c(r$row$coord[,1], r$col$coord[,1])
	yplot = c(r$row$coord[,2], r$col$coord[,2])
	l_pch = c(rep (20,length (r$row$coord[,1])), rep (18,length (r$col$coord[,2])))
	l_col_all = c(rep ("#000000", length (r$row$coord[,1])), defColor(names(r$col$coord[,1])))
	l_col = defColor(names(r$col$coord[,1]))
	
	lim_x = max (abs (xplot))
	lim_y = max (abs (yplot))
	
	
	svg (file = paste (path_file, "_AFC_text.svg", sep = ""), bg = "transparent" ,10, 10)
	par(mar=c(5,5,2,2))
	
	plot (xplot, yplot, xlab = paste("DIM 1 : ", round(r$eig[1,2],1), "%", sep = ""), ylab = paste("DIM 2 : ", round(r$eig[2,2],1), "%", sep = ""), cex.lab = 2.2, cex.axis = 1.8, xlim = c(-lim_x, lim_x), ylim = c(-lim_y, lim_y), pch = l_pch, col = l_col_all, cex = 2.5)
	text (r$col$coord[,1], r$col$coord[,2], label = names(r$col$coord[,1]), col = l_col, cex = 2, pos = 4)
	text (r$row$coord[,1], r$row$coord[,2], col = "black", label = names (r$row$coord[,1]), cex = 2.5, pos = 4)
	abline(h=0,v=0, lwd = 2)

	dev.off()
	
	svg (file = paste (path_file, "_AFC_point.svg", sep = ""), bg = "transparent", 10, 10)
	par(mar=c(5,5,2,2))
	plot (xplot, yplot, type = "n", xlab = paste("DIM 1 : ", round(r$eig[1,2],1), "%", sep = ""), ylab = paste("DIM 2 : ", round(r$eig[2,2],1), "%", sep = ""), cex.lab = 2.2, cex.axis = 1.8, xlim = c(-lim_x, lim_x), ylim = c(-lim_y, lim_y))
	points (r$col$coord[,1], r$col$coord[,2], col = l_col, cex = 4, pch = 16)
	text (r$row$coord[,1], r$row$coord[,2], col = "black", label = names (r$row$coord[,1]), cex = 2.5, pos = 4)
	abline(h=0,v=0, lwd = 2)
	
	dev.off()
	
	
	
}





#########
# color #
#########

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



defColor = function (l_name){

	out = c()
	colorCharge = "#FF0000"
	colorPolar = "#FF6600"
	colorAr = "#9900CC"
	colorOther = "#808080"
	colorApolar =  "#66FF33"
	colorCharge2 = "#000099"
	colorWater = "#33FFFF"
	colorTyr = "#FFFF00"
	
	for (element in l_name){
	  if (attr(regexpr("CI_1",element),"match.length") >= 4){
	    out = append (out, "#FF0000")
	  }
	  else if (attr(regexpr("CI_2",element),"match.length") >= 4){
	    out = append (out, "#FF0000")
	  }	  
	  else if (attr(regexpr("CI_3",element),"match.length") >= 4){
	    out = append (out, "#FF0000")
	  }	  
	  else if (attr(regexpr("CI_4",element),"match.length") >= 4){
	    out = append (out, "#FF0000")
	  }
	  else if (attr(regexpr("More_CI",element),"match.length") >= 4){
	    out = append (out, "#FF0000")
	  }
	  else if (attr(regexpr("Oox_long",element),"match.length") >= 4){
	    out = append (out, "#990000")
	  }
	  else if (attr(regexpr("Water_Mediated",element),"match.length") >= 6){
	    out = append (out, "#A67E2E")
	  }
	  else if (attr(regexpr("Other", element),"match.length") >= 4 ){
	    out = append (out, "#D3D3D3")
	  }
	  else if (attr(regexpr("Oph", element),"match.length") >= 3 ){
	    out = append (out, "#F9E79F")
	  }
	  else if (attr(regexpr("Su", element),"match.length") >= 1 ){
	    out = append (out, "#626567")
	  }
	  else if (attr(regexpr("Cox",element),"match.length") >= 1){
			out = append (out, "#5F021F")
	  }
		else if (attr(regexpr("Oox",element),"match.length") >= 1){
			out = append (out, "#FF0000")
		}
		else if (attr(regexpr("Oh",element, ignore.case = FALSE),"match.length") == 2){
		  print (element)
			out = append (out, "#FF6600")
		}
		else if (attr(regexpr("Sulfur",element),"match.length") >= 1){
			out = append (out, "#FFD700")
		}
		else if (attr(regexpr("Nim",element),"match.length") >= 2){
			out = append (out, "#7FB6E2")
		}
		else if (attr(regexpr("Ngu",element),"match.length") >= 2){
			out = append (out, "#006BFF")
		}
	  else if (attr(regexpr("NaI",element),"match.length") >= 1){
	    out = append (out, "#0000FF")
		}
		else if (attr(regexpr("Car",element),"match.length") >= 1){
			out = append (out, "#9900CC")
		}
		else if (attr(regexpr("Oc",element),"match.length") >= 1){
			out = append (out, "#000000")
		}
		else if (attr(regexpr("Nam",element),"match.length") >= 1){
			out = append (out, "#00b200")
		}
		else if (attr(regexpr("Oc",element),"match.length") >= 1){
			out = append (out, "#CDB79E")
		}
		else if (attr(regexpr("Ow",element),"match.length") >= 1){
			out = append (out, "#33FFFF")
		}
		else if (attr(regexpr("Xot",element),"match.length") >= 1){
			out = append (out, "#D3D3D3")
		}
		else if (attr(regexpr("Cgu",element),"match.length") >= 1){
			out = append (out, "#003366")
		}
	  else if (is.integer0 (grep("COO", element))== FALSE){
	    out = append (out, "#FF0000")
	  }
	  else if (attr(regexpr("H2O",element),"match.length") == 3){
	    out = append (out, "#33FFFF")
	  }
	  else if (attr(regexpr("HOH",element, ignore.case = FALSE),"match.length") == 3){
	    out = append (out, "#33FFFF")	    
	  }
	  else if (attr(regexpr("OC",element),"match.length") >= 1){
	    out = append (out, "#000000")
	  }

	  else if (attr(regexpr("OH",element, ignore.case = FALSE),"match.length") == 2){
	    out = append (out, "#FF6600")
	  }
	  else if (is.integer0 (grep("NH", element))== FALSE){
	    out = append (out, "#00b200")
	  }
	  else if (is.integer0 (grep("N", element))== FALSE){
	    out = append (out, "#0000FF")
	  }
	  else if (is.integer0 (grep("AR", element))== FALSE){
	    out = append (out, "#d869ff")	  
	  }
		else if (is.integer0 (grep("GLU", element))== FALSE){
			out = append (out, colorCharge)
		}
		else if (is.integer0 (grep("ASP", element))== FALSE){
			out = append (out, colorCharge)
		}	
		else if (is.integer0 (grep("ARG", element))== FALSE){
			out = append (out, colorCharge2)
		}	
		else if (is.integer0 (grep("THR", element))== FALSE){
			out = append (out, colorPolar)
		}	
		else if (is.integer0 (grep("SER", element))== FALSE){
			out = append (out, colorPolar)
		}	
		else if (is.integer0 (grep("CYS", element))== FALSE){
			out = append (out, colorApolar)
		}	
		else if (is.integer0 (grep("ASN", element))== FALSE){
			out = append (out, colorPolar)
		}
		else if (is.integer0 (grep("LYS", element))== FALSE){
			out = append (out, colorCharge2)
		}
		else if (is.integer0 (grep("GLN", element))== FALSE){
			out = append (out, colorPolar)
		}	
		else if (is.integer0 (grep("PHE", element))== FALSE){
			out = append (out, colorAr)
		}
		else if (is.integer0 (grep("HIS", element))== FALSE){
			out = append (out, colorTyr)
		}
		else if (is.integer0 (grep("TRP", element))== FALSE){
			out = append (out, colorAr)
		}
		else if (is.integer0 (grep("TYR", element))== FALSE){
			out = append (out, colorTyr)
		}
		else if (is.integer0 (grep("ALA", element))== FALSE){
			out = append (out, colorApolar)
		}
		else if (is.integer0 (grep("ILE", element))== FALSE){
			out = append (out, colorApolar)
		}
		else if (is.integer0 (grep("LEU", element))== FALSE){
			out = append (out, colorApolar)
		}
		else if (is.integer0 (grep("MET", element))== FALSE){
			out = append (out, colorApolar)
		}
		else if (is.integer0 (grep("PRO", element))== FALSE){
			out = append (out, colorOther)
		}
		else if (is.integer0 (grep("VAL", element))== FALSE){
			out = append (out, colorApolar)
		}
	  else if (is.integer0 (grep("GLY", element))== FALSE){
	    out = append (out, colorOther)
	  }	  
	  else if (is.integer0 (grep("E", element))== FALSE){
	    out = append (out, colorCharge)
	  }
	  else if (is.integer0 (grep("D", element))== FALSE){
	    out = append (out, colorCharge)
	  }	
	  else if (is.integer0 (grep("R", element))== FALSE){
	    out = append (out, colorCharge2)
	  }	
	  else if (is.integer0 (grep("T", element))== FALSE){
	    out = append (out, colorPolar)
	  }	
	  else if (is.integer0 (grep("S", element))== FALSE){
	    out = append (out, colorPolar)
	  }	
	  else if (is.integer0 (grep("C", element))== FALSE){
	    out = append (out, colorApolar)
	  }	
	  else if (is.integer0 (grep("N", element))== FALSE){
	    out = append (out, colorPolar)
	  }
	  else if (is.integer0 (grep("K", element))== FALSE){
	    out = append (out, colorCharge2)
	  }
	  else if (is.integer0 (grep("Q", element))== FALSE){
	    out = append (out, colorPolar)
	  }	
	  else if (is.integer0 (grep("F", element))== FALSE){
	    out = append (out, colorAr)
	  }
	  else if (is.integer0 (grep("H", element))== FALSE){
	    out = append (out, colorTyr)
	  }
	  else if (is.integer0 (grep("W", element))== FALSE){
	    out = append (out, colorAr)
	  }
	  else if (is.integer0 (grep("Y", element))== FALSE){
	    out = append (out, colorTyr)
	  }
	  else if (is.integer0 (grep("A", element))== FALSE){
	    out = append (out, colorApolar)
	  }
	  else if (is.integer0 (grep("I", element))== FALSE){
	    out = append (out, colorApolar)
	  }
	  else if (is.integer0 (grep("L", element))== FALSE){
	    out = append (out, colorApolar)
	  }
	  else if (is.integer0 (grep("M", element))== FALSE){
	    out = append (out, colorApolar)
	  }
	  else if (is.integer0 (grep("P", element))== FALSE){
	    out = append (out, colorOther)
	  }
	  else if (is.integer0 (grep("V", element))== FALSE){
	    out = append (out, colorApolar)
	  }
	  else if (is.integer0 (grep("G", element))== FALSE){
	    out = append (out, colorOther)
	  }
		else {
			out = append (out, "#000000")
		}
	}
	names(out) = l_name
	return (out)
}



defColorSubstruct = function (l_name){

	out = c()
	for (element in l_name){
		#print (element)
		if (element == "I"){
			out = append (out, "red")
		}
		else if (element == "II"){
			out = append (out, "orange")
		}
		else if (element == "III"){
			out = append (out, "yellow")
		}
		else if (element =="Diamine"){
			out = append (out, "cyan")
		}
		else if (element == "GAI"){
			out = append (out, "blue")
		}
		else if (element == "IMD"){
			out = append (out, "green")
		}
		else if (element == "Pyridine"){
			out = append (out, "purple")
		}
		else if (element == "COO"){
			out = append (out, "black")
		}
		else {
			out = append (out, "grey")
		}
	}
	names(out) = l_name
	return (out)
}


deviationAngle = function(matrixAngle){
	nbLine = dim(matrixAngle)[1]
	nb_col = dim(matrixAngle)[2]
	matrix_out = NULL
	for(i in seq(1,nbLine)){
		dev_angle = 0
		for (j in seq(1, nb_col-1)){
			#print (j)
			dev_angle = dev_angle + abs(120-matrixAngle[i,j])
			#print (dev_angle)
		}
		#print (matrixAngle[i,nb_col])
		line_table = c(matrixAngle[i,1],dev_angle,as.character(matrixAngle[i,nb_col]))
		matrix_out = rbind (matrix_out, line_table)
		}
	return (as.matrix(matrix_out))

}

###############################
# significiance of the pvalue #
###############################


signifPvalue = function (a){
	if (a < 0.001){
		return ("***")
	}
	else if (a < 0.01){
		return ("**")
	}
	else if (a < 0.05){
		return ("*")
	}
	else {
		return ("-")
	}	
}


#####################
# useful function   #
#####################


GetPercent = function (d_in, sum_in){
  
  d_out = data.frame()
  for (i in seq (1, dim(d_in)[1])){
    if (sum_in == 0){
      sum_subs = sum (d_in [i,])
      print (sum_subs)
    }else{
      sum_subs = sum_in
    }
    for (j in seq (1, dim(d_in)[2])){
      d_out[i, j] = (d_in[i, j] / sum_subs)*100
    }
  }
  return (d_out)
}


MiddlePoint = function (v1, v2){
  
  difference = v2-v1
  
  return (v1 + (difference/2))
  
}


