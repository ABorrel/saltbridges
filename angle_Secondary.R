#!/usr/bin/env Rscript

classAngle = function(data){
	nbCol = dim(data)[2]
	nbLine = dim(data)[1]

	for (i in seq(1,nbLine)){

		if (data[i,2] < data[i,3]){
			key = data[i,2]
			data[i,2] = data[i,3]
			data[i,3] = key
		}
	}
	return (data)
}

deviation = function(matrixAngle){
	nbLine = dim(matrixAngle)[1]
	diff = NULL
	for(i in seq(1,nbLine)){
		diff = cbind(diff,(matrixAngle[i,2]-matrixAngle[i,3])^2)
	}
	return (diff)

}


#######################
#        MAIN         #
#######################

library(lattice,lib.loc="/home/student10/R_libs/")
library(scatterplot3d,lib.loc="/home/student10/R_libs/")

require(methods)

args <- commandArgs(TRUE)
file = args[1]

data = read.csv(file, sep = "\t", header = FALSE)

#data = classAngle(data)


png(paste(file,"_angle3D.png", sep = ""))
cloud(data[,1] ~ data[,2] * data[,3], screen = list(z = 40, x = 650, y=20), xlab = "distance Å", ylab = "angle2", zlab = "angle3")
dev.off()




listOxAcid = c()
listOxAmphi = c()
listNbasic = c()
listNdonnor = c()
listCar = c()
listOUT = c()
listOxAccept = c()
listH2O = c()

nbLine = dim(data)[1]

for (i in seq(1,nbLine)){
	if (data[i,4] == "OxAcid"){
		listOxAcid = rbind(listOxAcid,data[i,])
	}
	if (data[i,4] == "amphiprotic"){
		listOxAmphi = rbind(listOxAmphi,data[i,])
	}
	if (data[i,4] == "Nbasic"){
		listNbasic = rbind(listNbasic,data[i,])
	}
	if (data[i,4] == "Ndonnor"){
		listNdonnor = rbind(listNdonnor,data[i,])
	}
	if (data[i,4] == "Carom"){
		listCar = rbind(listCar,data[i,])
	}
	if (data[i,4] == "others"){
		listOUT = rbind(listOUT,data[i,])
	}
	if (data[i,4] == "OxAccept"){
		listOxAccept = rbind(listOxAccept,data[i,])
	}
	if (data[i,4] == "H2O"){
		listH2O = rbind(listH2O,data[i,])
	}
}


calibrate = data[,1]
calibrate = c(calibrate,data[,1])

temp = data[,2]
temp = c(temp, data[,3])

calibrate = cbind(calibrate,temp) 

png(filename = paste(file , "_typeAtom.png" , sep = ""), width=800, height = 1200)
par(mfrow=c(2,4))


plot(data[,1], data[,2], xlab = "Distance Å", ylab = "Angles en degres", type = "n", xlim = c(1.5,5))
try(points(listOxAmphi[,2]~listOxAmphi[,1], col = "orange", type = "p"))
try(points(listOxAmphi[,3]~listOxAmphi[,1], col = "orange", type = "p"))

plot(data[,1], data[,2], xlab = "Distance Å", ylab = "Angles en degres", type = "n", xlim = c(1.5,5))
try(points(listOxAcid[,3]~listOxAcid[,1], col = "red", type = "p"))
try(points(listOxAcid[,2]~listOxAcid[,1], col = "red", type = "p"))

plot(data[,1], data[,2], xlab = "Distance Å", ylab = "Angles en degres", type = "n", xlim = c(1.5,5))
try(points(listOxAccept[,3]~listOxAccept[,1], col = "cyan", type = "p"))
try(points(listOxAccept[,2]~listOxAccept[,1], col = "cyan", type = "p"))

plot(data[,1], data[,2], xlab = "Distance Å", ylab = "Angles en degres", type = "n", xlim = c(1.5,5))
try(points(listH2O[,3]~listH2O[,1], col = "yellow", type = "p"))
try(points(listH2O[,2]~listH2O[,1], col = "yellow", type = "p"))


plot(data[,1], data[,2], xlab = "Distance Å", ylab = "Angles en degres", type = "n", xlim = c(1.5,5))
try(points(listNdonnor[,3]~listNdonnor[,1], col = "green", type = "p"))
try(points(listNdonnor[,2]~listNdonnor[,1], col = "green", type = "p"))

plot(data[,1], data[,2], xlab = "Distance Å", ylab = "Angles en degres", type = "n", xlim = c(1.5,5))
try(points(listNbasic[,3]~listNbasic[,1], col = "blue", type = "p"))
try(points(listNbasic[,2]~listNbasic[,1], col = "blue", type = "p"))


plot(data[,1], data[,2], xlab = "Distance Å", ylab = "Angles en degres", type = "n", xlim = c(1.5,5))
try(points(listCar[,3]~listCar[,1], col = "purple", type = "p"))
try(points(listCar[,2]~listCar[,1], col = "purple", type = "p"))

plot(data[,1], data[,2], xlab = "Distance Å", ylab = "Angles en degres", type = "n", xlim = c(1.5,5))
try(points(listOUT[,3]~listOUT[,1], col = "grey", type = "p"))
try(points(listOUT[,2]~listOUT[,1], col = "grey", type = "p"))




color = c("red","orange","yellow","cyan","blue","green","purple","grey")
nameGroup = c( "O (COOH)", "O (Tyr, SER, THR), S (CYS)","O (H2O)", "O (main chain) Side chain ASN, GLN", "N (HIS, LYS, ARG) and Nxt", "N (main chain) ASN, GLN", "Others", "C (side chain TYR, PHE, TRP)")
legend( 1.5 , 85, legend=nameGroup, bty="n", fill=color)

dev.off()


