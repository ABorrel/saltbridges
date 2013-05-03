#!/usr/bin/env Rscript


classAngleDistance = function(data){

	nbLine = dim(data)[1]
	vect02 = NULL
	vect225 = NULL
	vect253 = NULL
	vect335 = NULL
	vect354 = NULL
	vect445 = NULL

	for (i in seq(0,nbLine)){
		if (data[i,]<=2){
		vect02 = cbind(vect02,data[i,])
		}
		else if (data[i,]<=2.5){
		vect225 = cbind(vect225,data[i,])
		}
		else if (data[i,]<=3){
		vect253 = cbind(vect253,data[i,])
		}
		else if (data[i,]<=3.5){
		vect335 = cbind(vect335,data[i,])
		}
		else if (data[i,]<=4){
		vect354 = cbind(vect354,data[i,])
		}
		else {
		vect445 = cbind(vect445,data[i,])
		}
	}
	return (vect02,vect225,vect253,vect335,vect354,vect445)
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

for (i in seq(1,nbLine)){
	if (data[i,3] == "OxAcid"){
		listOxAcid = rbind(listOxAcid,data[i,])
	}
	if (data[i,3] == "amphiprotic"){
		listOxAmphi = rbind(listOxAmphi,data[i,])
	}
	if (data[i,3] == "Nbasic"){
		listNbasic = rbind(listNbasic,data[i,])
	}
	if (data[i,3] == "Ndonnor"){
		listNdonnor = rbind(listNdonnor,data[i,])
	}
	if (data[i,3] == "Carom"){
		listCar = rbind(listCar,data[i,])
	}
	if (data[i,3] == "others"){
		listOUT = rbind(listOUT,data[i,])
	}
	if (data[i,3] == "OxAccept"){
		listOxAccept = rbind(listOxAccept,data[i,])
	}
	if (data[i,3] == "H2O"){
		listH2O = rbind(listH2O,data[i,])
	}
}


png(filename = paste(file , ".png" , sep = ""), width=400, height = 200)
par(mfrow=c(2,4))

plot(data[,1], data[,2], xlab = "Distance Å", ylab = "Angles en degres", type = "n", xlim = c(1.5,5))
try(points(listOxAmphi[,2]~listOxAmphi[,1], col = "orange", type = "p"))

plot(data[,1], data[,2], xlab = "Distance Å", ylab = "Angles en degres", type = "n", xlim = c(1.5,5))
try(points(listOxAcid[,2]~listOxAcid[,1], col = "red", type = "p"))

plot(data[,1], data[,2], xlab = "Distance Å", ylab = "Angles en degres", type = "n", xlim = c(1.5,5))
try(points(listH2O[,2]~listH2O[,1], col = "yellow", type = "p"))

plot(data[,1], data[,2], xlab = "Distance Å", ylab = "Angles en degres", type = "n", xlim = c(1.5,5))
try(points(listOxAccept[,2]~listOxAccept[,1], col = "cyan", type = "p"))

plot(data[,1], data[,2], xlab = "Distance Å", ylab = "Angles en degres", type = "n", xlim = c(1.5,5))
try(points(listNbasic[,2]~listNbasic[,1], col = "blue", type = "p"))

plot(data[,1], data[,2], xlab = "Distance Å", ylab = "Angles en degres", type = "n", xlim = c(1.5,5))
try(points(listNdonnor[,2]~listNdonnor[,1], col = "green", type = "p"))

plot(data[,1], data[,2], xlab = "Distance Å", ylab = "Angles en degres", type = "n", xlim = c(1.5,5))
try(points(listOUT[,2]~listOUT[,1], col = "grey", type = "p"))

plot(data[,1], data[,2], xlab = "Distance Å", ylab = "Angles en degres", type = "n", xlim = c(1.5,5))
try(points(listCar[,2]~listCar[,1], col = "purple", type = "p"))



color = c("red","orange","yellow","cyan","blue","green","grey","purple")
nameGroup = c( "O (COOH)", "O (Tyr, SER, THR), S (CYS)","O (H2O)", "O (main chain) Side chain ASN, GLN", "N (HIS, LYS, ARG) and Nxt", "N (main chain) ASN, GLN", "Others", "C (side chain TYR, PHE, TRP)")
legend( 1.5 , 85, legend=nameGroup, bty="n", fill=color)

dev.off()



brk = seq(0,180,20)

png(filename=paste(file, "_distanceDistribution.png", sep = ""))

hist(data[,2], xlab ="Distances (Å)", ylab = "Number of occurences", breaks=brk, right=F, main=paste("Distribution of angles","\n", "Primary amines", sep = ""),las=1, freq=T, col = "#D8D8D8")

dev.off()
