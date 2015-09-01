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

	print (dim(d))

	r = CA (d, graph = FALSE)

	svg (file = paste (path_file, "_AFC.svg", sep = ""), 15, 15)
	par(mar=c(8,8,8,8))

	# descriptors
	plot (c(r$row$coord[,1], r$col$coord[,1]), c(r$row$coord[,2], r$col$coord[,2]), type = "n", xlab = paste("DIM 1 : ", round(r$eig[1,2],1), "%", sep = ""), ylab = paste("DIM 2 : ", round(r$eig[2,2],1), "%", sep = ""), cex.lab = 2.4)
	col_des = defColor(names(r$col$coord[,1]))
	
	print (col_des)	
	
	text (r$col$coord[,1], r$col$coord[,2], label = names(r$col$coord[,1]), col = col_des, cex = 2)
	# data
	text (r$row$coord[,1], r$row$coord[,2], col = "black", label = names (r$row$coord[,1]), cex = 2.5)
	abline(h=0,v=0)

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
		if (attr(regexpr("OxAcid",element),"match.length") >= 1){
			out = append (out, "#FF0000")
		}
		else if (attr(regexpr("ODonAcc",element),"match.length") >= 1){
			out = append (out, "#FF6600")
		}
		else if (attr(regexpr("Sulfur",element),"match.length") >= 1){
			out = append (out, "#FFD700")
		}
		else if (attr(regexpr("Nhis",element),"match.length") >= 1){
			out = append (out, "#000080")
		}
		else if (attr(regexpr("Nbasic",element),"match.length") >= 1){
			out = append (out, "#0000FF")
		}
		else if (attr(regexpr("Carom",element),"match.length") >= 1){
			out = append (out, "#9900CC")
		}
		else if (attr(regexpr("OxPep",element),"match.length") >= 1){
			out = append (out, "#FF7F50")
		}
		else if (attr(regexpr("Ndonnor",element),"match.length") >= 1){
			out = append (out, "#00b200")
		}
		else if (attr(regexpr("OxAccept",element),"match.length") >= 1){
			out = append (out, "#CDB79E")
		}
		else if (attr(regexpr("H2O",element),"match.length") >= 1){
			out = append (out, "#33FFFF")
		}
		else if (attr(regexpr("CPep",element),"match.length") >= 1){
			out = append (out, "#663300")
		}
		else if (attr(regexpr("others",element),"match.length") >= 1){
			out = append (out, "#D3D3D3")
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
		else {
			out = append (out, "#000000")
		}
	}
	names(out) = l_name
	return (out)
}



defColorGrep = function (l_name){
	print (l_name)
	out = c()

	for (element in l_name){
		#print (element)
		if (is.integer0 (grep("OxAcid", element)) == FALSE){
			out = append (out, "red")
		}
		else if (is.integer0 ( grep("ODonAcc", element))== FALSE){
			out = append (out, "orange")
		}
		else if (is.integer0 (grep("Sulfur", element))== FALSE){
			out = append (out, "gold")
		}
		else if (is.integer0 (grep("Nhis", element))== FALSE){
			out = append (out, "darkblue")
		}
		else if (is.integer0 (grep("Nbasic", element))== FALSE){
			out = append (out, "blue")
		}
		else if (is.integer0 ( grep("Carom", element))== FALSE){
			out = append (out, "purple")
		}
		else if (is.integer0 (grep("OxPep", element))== FALSE){
			out = append (out, "coral")
		}
		else if (is.integer0 (grep("Ndonnor", element))== FALSE){
			out = append (out, "green")
		}
		else if (is.integer0 (grep("OxAccept", element))== FALSE){
			out = append (out, "bisque3")
		}
		else if (is.integer0 (grep("H2O", element))== FALSE){
			out = append (out, "cyan")
		}
		else if (is.integer0 (grep("CPep", element))== FALSE){
			out = append (out, "brown")
		}
		else if (is.integer0 (grep("others", element))== FALSE){
			out = append (out, "grey")
		}
		else {
			out = append (out, "black")
		}
	}
	names(out) = l_name
	return (out)
}


defColorSubstruct = function (l_name){

	out = c()
	for (element in l_name){
		#print (element)
		if (element == "Primary"){
			out = append (out, "red")
		}
		else if (element == "Secondary"){
			out = append (out, "orange")
		}
		else if (element == "Tertiary"){
			out = append (out, "yellow")
		}
		else if (element =="Diamine"){
			out = append (out, "cyan")
		}
		else if (element == "Guanidium"){
			out = append (out, "blue")
		}
		else if (element == "Imidazole"){
			out = append (out, "green")
		}
		else if (element == "Pyridine"){
			out = append (out, "purple")
		}
		else if (element == "AcidCarboxylic"){
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

