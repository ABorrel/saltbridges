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



listOxAcid = c()
listOxAmphi = c()
listNbasic = c()
listNdonnor = c()
listCar = c()
listOUT = c()
listOxAccept = c()
listH2O = c()

nbLine = dim(data)[1]
nb_col = dim(data)[2]

for (i in seq(1,nbLine)){
	if (data[i,nb_col] == "OxAcid"){
		listOxAcid = rbind(listOxAcid,data[i,])
	}
	if (data[i,nb_col] == "amphiprotic"){
		listOxAmphi = rbind(listOxAmphi,data[i,])
	}
	if (data[i,nb_col] == "Nbasic"){
		listNbasic = rbind(listNbasic,data[i,])
	}
	if (data[i,nb_col] == "Ndonnor"){
		listNdonnor = rbind(listNdonnor,data[i,])
	}
	if (data[i,nb_col] == "Carom"){
		listCar = rbind(listCar,data[i,])
	}
	if (data[i,nb_col] == "others"){
		listOUT = rbind(listOUT,data[i,])
	}
	if (data[i,nb_col] == "OxAccept"){
		listOxAccept = rbind(listOxAccept,data[i,])
	}
	if (data[i,nb_col] == "H2O"){
		listH2O = rbind(listH2O,data[i,])
	}
}

# angle vs angle
if (dim (data)[2] == 4){
	png(filename = paste(file , "angle_density.png" , sep = ""), width=4000, height = 2000)
	par(mar=c(6,6,6,6))
	par(mfrow=c(2,4))
	Lab.palette = colorRampPalette(revector(heat.colors(50)), space = "Lab")
	smoothScatter(listOxAcid[,c(2,3)], colramp = Lab.palette, main = "OxAcid", cex.axis = 2.5, cex.main = 3, xlab = "angle1", ylab = "angle2", cex.lab = 3)
	smoothScatter(listOxAmphi[,c(2,3)], colramp = Lab.palette, main = "amphiprotic", cex.axis = 2.5, cex.main = 3, xlab = "angle1", ylab = "angle2", cex.lab = 3)
	smoothScatter(listNbasic[,c(2,3)], colramp = Lab.palette, main = "Nbasic", cex.axis = 2.5, cex.main = 3, xlab = "angle1", ylab = "angle2", cex.lab = 3)
	smoothScatter(listNdonnor[,c(2,3)], colramp = Lab.palette, main = "Ndonnor", cex.axis = 2.5, cex.main = 3, xlab = "angle1", ylab = "angle2", cex.lab = 3)
	smoothScatter(listCar[,c(2,3)], colramp = Lab.palette, main = "Carom", cex.axis = 2.5, cex.main = 3, xlab = "angle1", ylab = "angle2", cex.lab = 3)
	smoothScatter(listOUT[,c(2,3)], colramp = Lab.palette, main = "others", cex.axis = 2.5, cex.main = 3, xlab = "angle1", ylab = "angle2", cex.lab = 3)
	smoothScatter(listOxAccept[,c(2,3)], colramp = Lab.palette, main = "OxAccept", cex.axis = 2.5, cex.main = 3, xlab = "angle1", ylab = "angle2", cex.lab = 3)
	smoothScatter(listH2O[,c(2,3)], colramp = Lab.palette, main = "H2O", cex.axis = 2.5, cex.main = 3, xlab = "angle1", ylab = "angle2", cex.lab = 3)
	dev.off()
}




mapFig ("OxAcid", listOxAcid, file)
mapFig ("amphiprotic", listOxAmphi, file)
mapFig ("Nbasic", listNbasic, file)
mapFig ("Ndonnor", listNdonnor, file)
mapFig ("Carom", listCar, file)
mapFig ("others", listOUT, file)
mapFig ("OxAccept", listOxAccept, file)
mapFig ("H2O", listH2O, file)



png(filename = paste(file , "density.png" , sep = ""), width=4000, height = 2000)
par(mfrow=c(2,4))
Lab.palette = colorRampPalette(revector(heat.colors(50)), space = "Lab")

smoothScatter(listOxAcid, colramp = Lab.palette, main = "OxAcid", xlim = c(1.5,5), cex.axis = 2.5, cex.main = 3, xlab = "distance", ylab = "angle2", cex.lab = 3)
smoothScatter(listOxAmphi, colramp = Lab.palette, main = "amphiprotic", xlim = c(1.5,5), cex.axis = 2.5, cex.main = 3, xlab = "distance", ylab = "angle2", cex.lab = 3)
smoothScatter(listNbasic, colramp = Lab.palette, main = "Nbasic", xlim = c(1.5,5), cex.axis = 2.5, cex.main = 3, xlab = "distance", ylab = "angle2", cex.lab = 3)
smoothScatter(listNdonnor, colramp = Lab.palette, main = "Ndonnor", xlim = c(1.5,5), cex.axis = 2.5, cex.main = 3, xlab = "distance", ylab = "angle2", cex.lab = 3)
smoothScatter(listCar, colramp = Lab.palette, main = "Carom", xlim = c(1.5,5), cex.axis = 2.5, cex.main = 3, xlab = "distance", ylab = "angle2", cex.lab = 3)
smoothScatter(listOUT, colramp = Lab.palette, main = "others", xlim = c(1.5,5), cex.axis = 2.5, cex.main = 3, xlab = "distance", ylab = "angle2", cex.lab = 3)
smoothScatter(listOxAccept, colramp = Lab.palette, main = "OxAccept", xlim = c(1.5,5), cex.axis = 2.5, cex.main = 3, xlab = "distance", ylab = "angle2", cex.lab = 3)
smoothScatter(listH2O, colramp = Lab.palette, main = "H2O", xlim = c(1.5,5), cex.axis = 2.5, cex.main = 3, xlab = "distance", ylab = "angle2", cex.lab = 3)

dev.off()

