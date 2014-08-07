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

	png(filename = paste(file, main_plot, "_density.png" , sep = ""), width=6000, height = 6000)
	par(mar=c(30,30,30,30))
	Lab.palette = colorRampPalette(revector(heat.colors(50)), space = "Lab")

	#if (main_plot == "Imidazole"){

		smoothScatter(data_bis, col = blues9, main = main_plot, xlim = c(1.5,4), ylim = c(0, 180), cex.lab = 10, xlab = "", ylab = "", cex.axis = 10, cex.main = 6)
		# vertical	
		segments (2,0,2,180,lwd = 6)
		segments (2.5,0,2.5,180,lwd = 6)
		segments (3,0,3,180,lwd = 6)
		segments (3.5,0,3.5,180,lwd = 6)
	
		# horizontal
		segments (1.5,0,4,0,lwd = 6)
		segments (1.5,20,4,20,lwd = 6)
		segments (1.5,40,4,40,lwd = 6)
		segments (1.5,60,4,60,lwd = 6)
		segments (1.5,80,4,80,lwd = 6)
		segments (1.5,100,4,100,lwd = 6)
		segments (1.5,120,4,120,lwd = 6)
		segments (1.5,140,4,140,lwd = 6)
		segments (1.5,160,4,160,lwd = 6)



	#}else {

		#smoothScatter(data_bis, col = blues9, main = main_plot, xlim = c(1.5,4), ylim = c(40, 180), cex.lab = 10, xlab = "", ylab = "", cex.axis = 10, cex.main = 6)
		# vertical	
		#segments (2,40,2,180,lwd = 6)
		#segments (2.5,40,2.5,180,lwd = 6)
		#segments (3,40,3,180,lwd = 6)
		#segments (3.5,40,3.5,180,lwd = 6)
	
		# horizontal
		#segments (1.5,40,4,40,lwd = 6)
		#segments (1.5,60,4,60,lwd = 6)
		#segments (1.5,80,4,80,lwd = 6)
		#segments (1.5,100,4,100,lwd = 6)
		#segments (1.5,120,4,120,lwd = 6)
		#segments (1.5,140,4,140,lwd = 6)
		#segments (1.5,160,4,160,lwd = 6)
	#}


	#grid (5, 7, col = "black", lwd = 6, lty = 1)
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
		par(mar=c(10,10,10,10))
		Lab.palette = colorRampPalette(revector(heat.colors(50)), space = "Lab")
		smoothScatter(data_plot[,c(2,3)], colramp = Lab.palette, main = element, cex.axis = 2.5, cex.main = 3, xlab = "angle1", ylab = "angle2", cex.lab = 3)
		dev.off()
	}
}

