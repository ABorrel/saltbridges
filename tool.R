#!/usr/bin/env Rscript

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
