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
	#return (1)
	factor = 1
	while ( abs(max(xplot)) < max(abs(xplotdata)) && abs(max(yplot)) < max(abs(yplotdata)) ){
		factor = factor + 0.2
		xplot = xplot * factor
		yplot = yplot* factor

	}
	return (factor)	
}



AFC = function (d, path_file){

	#print (d)
	r = CA (d, graph = FALSE)

	svg (file = paste (path_file, "_AFC.svg", sep = ""), 15, 15)
	par(mar=c(8,8,8,8))

	# descriptors
	plot (r$row$coord[,1], r$row$coord[,2], type = "n", xlab = paste("DIM 1 : ", round(r$eig[1,2],1), "%", sep = ""), ylab = paste("DIM 2 : ", round(r$eig[2,2],1), "%", sep = ""), cex.lab = 2.4, ylim = c(-max(r$row$coord[,2]), max(r$row$coord[,2])), xlim = c(-max(r$row$coord[,1]), max(r$row$coord[,1])))
	col_des = defColor(names(r$col$coord[,1]))

	factor = factorAFC (r$col$coord[,1], r$col$coord[,2], r$row$coord[,1], r$row$coord[,2] )

	arrows (0,0,r$col$coord[,1]*factor, r$col$coord[,2]*factor, col = as.character(col_des), lwd = 3 )
	#text (r$col$coord[,1]*factor, r$col$coord[,2]*factor, labels = names(r$col$coord[,1]), col = as.character(col_des), cex = 1.4)

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
			out = append (out, "lemonchiffon")
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


