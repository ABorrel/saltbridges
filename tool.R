#!/usr/bin/env Rscript

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

	#print (d)
	r = CA (d, graph = FALSE)

	print (r)

	svg (file = paste (path_file, "_AFC.svg", sep = ""), 15, 15)
	par(mar=c(8,8,8,8))

	# descriptors
	plot (r$row$coord[,1], r$row$coord[,2], type = "n", xlab = paste("DIM 1 : ", round(r$eig[1,2],1), "%", sep = ""), ylab = paste("DIM 2 : ", round(r$eig[2,2],1), "%", sep = ""), cex.lab = 2.4, ylim = c(-max(abs(r$row$coord[,2])), max(abs(r$row$coord[,2]))), xlim = c(-max(abs(r$row$coord[,1])), max(abs(r$row$coord[,1]))))
	col_des = defColor(names(r$col$coord[,1]))
	print (col_des)
	factor = factorAFC (r$col$coord[,1], r$col$coord[,2], r$row$coord[,1], r$row$coord[,2] )

	arrows (0,0,r$col$coord[,1]*factor, r$col$coord[,2]*factor, col = as.character(col_des), lwd = 3 )

	# data
	text (r$row$coord[,1], r$row$coord[,2], col = "#009DE0", label = names (r$row$coord[,1]), cex = 2)
	abline(h=0,v=0)

	nameGroup = colnames(data)
	name_angle = rownames(data)
	color = defColor (nameGroup)
	legend( "bottomleft" , legend=names(r$col$coord[,1]), bty="n", fill=col_des)

	dev.off()
}





#########
# color #
#########


defColor = function (l_name){
	out = c()
	for (element in l_name){
		#print (element)
		if (element == "OxAcid"){
			out = append (out, "red")
		}
		else if (element == "ODonAcc"){
			out = append (out, "orange")
		}
		else if (element == "Sulfur"){
			out = append (out, "gold")
		}
		else if (element == "Nhis"){
			out = append (out, "darkblue")
		}
		else if (element == "Nbasic"){
			out = append (out, "blue")
		}
		else if (element == "Carom"){
			out = append (out, "purple")
		}
		else if (element == "OxPep"){
			out = append (out, "coral")
		}
		else if (element == "Ndonnor"){
			out = append (out, "green")
		}
		else if (element == "OxAccept"){
			out = append (out, "bisque3")
		}
		else if (element == "H2O"){
			out = append (out, "cyan")
		}
		else if (element == "CPep"){
			out = append (out, "brown")
		}
		else if (element == "others"){
			out = append (out, "grey")
		}
		else {
			out = append (out, "black")
		}
	}
	names(out) = l_name
	return (out)
}



defColorGrep = function (l_name){
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


is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}

