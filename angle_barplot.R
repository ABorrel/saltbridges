#!/usr/bin/env Rscript

source ("tool.R")

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
	
	# Read table
	data = read.table(file, sep = "\t",header = T)
	nameGroup = colnames(data)
	name_angle = rownames(data)
	color = defColor (nameGroup)

	data = data[order (as.double(name_angle)),] # order with angle values
	barplot(t(data), ylab= "Number of occurences", col = color, xlab = "Angles (degres)", names.arg = rownames (data),space = 0.5, main = paste(distanceTemp," < d < ",distance,sep = ""))
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
