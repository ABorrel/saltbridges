#!/usr/bin/env Rscript



args <- commandArgs(TRUE)
fileGlobal = args[1]
distanceMax = as.double(args[2])


listDistance = NULL

for (i in seq(2.00,distanceMax,0.5)){
	listDistance = c(listDistance,i)

}

nbDistance = length(listDistance)


nbGrapheLine = as.integer((nbDistance+1)/2)

png (paste(fileGlobal,"_distance.png", sep = ""),width=1000, height = nbGrapheLine * 300)
par(mfrow=c(nbGrapheLine,2))

color = c("red","orange","yellow","cyan","blue","green","purple","grey")
nameGroup = c( "O (COOH)", "O (Tyr, SER, THR), S (CYS)","O (H2O)", "O (main chain) Side chain ASN, GLN", "N (HIS, LYS, ARG) and Nxt", "N (main chain) ASN, GLN", "C (side chain TYR, PHE, TRP)", "Others")


i=0
distanceTemp = 0
for (distance in listDistance){
	i = i +1 
	distance = as.character(distance)
	lengthChar = nchar(distance)
	if (lengthChar == 1){
		distance = paste(distance,".0",sep = "")
	}
	file = paste(fileGlobal, "_",distance, sep = "")

	data = read.table(file, sep = "\t",header = F)
	data = data[order(data[,1],decreasing = F),]
	nameX = data[,1]
	data = data[,-1]
	barplot(t(data), ylab= "Number of occurences", col = color, xlab = "Angles (degres)", names.arg = nameX,space = 0.5, main = paste(distanceTemp," < d < ",distance,sep = ""))
	if (i == 3){
		yLegend = 0
		for(i in seq(1,dim(data)[1])){
			sumLine = sum(data[i,])
			if (sumLine > yLegend){
				yLegend = sumLine
			}

		}
		legend( 2 , yLegend , legend=nameGroup, bty="n", fill=color)

	}
	distanceTemp = distance

}

dev.off()

warnings()
