#!/usr/bin/env Rscript


revector<-function(n){
	out = NULL
	nb_value = length (n)
	for(i in nb_value:1){
		out = append (out,n[i])
	}
	return(out)
}


mapFig = function (main_plot, data, filin){
	data_bis = NULL
	if (dim(data)[2] > 3 ){
		nb_line = dim (data)[1]
		nb_col = dim (data)[2]
		for (i in seq (1,nb_line)){
			a = c()
			for (j in 2:nb_col-1){
				a = append(a,data[i,j])
			}
		data_bis = rbind (data_bis, c(data[i,1], mean (a)))
		}
	}
	else {
		data_bis = data
	}

	png(filename = paste(file, main_plot, "_density.png" , sep = ""), width=2000, height = 2000)
	par(mar=c(6,6,6,6))
	Lab.palette = colorRampPalette(revector(heat.colors(50)), space = "Lab")
	smoothScatter(data_bis, colramp = Lab.palette, main = main_plot, xlim = c(1.5,5), cex.axis = 2.5, cex.main = 3, xlab = "distance", ylab = "angle", cex.lab = 3)
	dev.off()
}



#######################
#        MAIN         #
#######################

args <- commandArgs(TRUE)
file = args[1]

data = read.csv(file, sep = "\t", header = FALSE)


nbLine = dim(data)[1]
nb_col = dim(data)[2]


l_element = unique(data[,nb_col])

for (element in l_element){
	
	data_plot = data[which(data[,nb_col] == element),]
	try(mapFig (element, data_plot, file))

	if (dim (data)[2] == 4){
		png(filename = paste(file , element, "angle_density.png" , sep = ""), width=4000, height = 2000)
		par(mar=c(6,6,6,6))
		Lab.palette = colorRampPalette(revector(heat.colors(50)), space = "Lab")
		smoothScatter(data_plot[,c(2,3)], colramp = Lab.palette, main = element, cex.axis = 2.5, cex.main = 3, xlab = "angle1", ylab = "angle2", cex.lab = 3)
		dev.off()
	}
}

