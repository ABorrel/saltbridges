#!/usr/bin/env Rscript



modifMatrixProportion = function(file){

	data = read.table(file, header = FALSE, sep = "\t")
	
	nbLine = dim(data)[1]
	nbCol = dim(data)[2]
	data = data[,-nbCol]
	nbCol = nbCol -1	
	
	
	for (col in seq(1,nbCol,2)){
		for (line in seq(1,nbLine)){
		data[line,col] = data[line,col] /(data[line,col]+data[line,col+1]) 
		}
	
	}

	for (i in seq(2,nbCol/2 +1)){
		data = data[,-i]
	}

	return (data)
}



###########################################
#               MAIN                      #
###########################################


args <- commandArgs(TRUE)
file = args[1]
fileGlobal = args[2]
distanceMax = as.numeric(args[3])

print (distanceMax)

listDistance = NULL

for (i in seq(2.00,distanceMax,0.5)){
	listDistance = c(listDistance,i)
}

listName = c()
for (distance in listDistance){
	print (distance)
	distance = paste("<",distance,"Å", sep = "")
	listName = c(listName,distance)
}

#distance = c("<2Å", "<2.5Å", "<3Å", "<3.5Å", "<4Å", "<4.5Å")
nameGroup = c("Primary", "Secondary", "Tertiary", "Imidazole", "All atoms")
color = c("red","blue", "yellow", "green", "gray")

data = modifMatrixProportion(file)
global = modifMatrixProportion(fileGlobal)
data = rbind(data,global)
data = cbind(data,nameGroup)

#data = data[order(data[,4],decreasing = T),]

nbCol = dim(data)[2]
#nameGroupe = data[,nbCol]
data = data[,-nbCol]

large = 600

png(filename=paste(file,".png",sep = ""),width=as.integer(large))

barplot(as.matrix(data), main=paste("Percentage of aminergic and imidazole","\n"," Nitrogen atoms (Nref)at least one carboxylic oxygen","\n","at selected distances", sep = ""), ylab= "% of nitrogen atoms", beside=TRUE, col=color, xlab = "Oxygen (COO-) to nitrogen (Nref) distances (Å)", names.arg = listName, ylim = c(0,1))

legend( 2 , 0.9 , legend=nameGroup, bty="n", fill=color)

dev.off()
