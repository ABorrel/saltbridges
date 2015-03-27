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

# case area
# limit 3.5

limit_zone = 3.5
d_area1 = 0
d_area2 = 0

for (distance in listDistance){
	distance = as.character(distance)
	lengthChar = nchar(distance)
	if (lengthChar == 1){
		distance = paste(distance,".0",sep = "")
	}

	file = paste(fileGlobal, "_",distance, sep = "")
	
	# Read table
	d = read.table(file, sep = "\t",header = T)
	
	if (distance <= limit_zone){
		
		print (distance)
		if (d_area1 == 0){
			d_area1 = d
		}
		else {
			d_area1 = d_area1 + d
		}
	}else{
		print (distance)
		if (d_area2 == 0){
			d_area2 = d
		}
		else {
			d_area2 = d_area2 + d
		}
	}
}


nameGroup = colnames(d_area1)
name_angle = rownames(d_area1)
color = defColor (nameGroup)

d_area1 = d_area1[order (as.double(name_angle)),] # order with angle values
d_area2 = d_area2[order (as.double(name_angle)),] # order with angle values	


svg (paste(fileGlobal,"_area1_", limit_zone, ".svg", sep = ""), width = 10, height = 12)
barplot(t(d_area1), ylab= "Number of occurences", col = color, xlab = "Angles (degres)", names.arg = rownames (d_area1),space = 0.5, main = "")
legend( "right", legend = nameGroup, bty="n", fill=color)
dev.off ()

svg (paste(fileGlobal, "_area2_", limit_zone, ".svg", sep = ""), width = 10, height = 12)
barplot(t(d_area2), ylab= "Number of occurences", col = color, xlab = "Angles (degres)", names.arg = rownames (d_area2),space = 0.5, main = "")
legend( "right", legend = nameGroup, bty="n", fill=color)
dev.off ()




