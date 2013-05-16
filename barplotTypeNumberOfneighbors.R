#!/usr/bin/env Rscript



source ("tool.R")

args <- commandArgs(TRUE)
file = args[1]
type = args[2]
distance = args[3]
flag = 0

if(type == "Global"){
	nameGraphe = paste ("Surronding of all atoms (",distance, " Å)", sep = "")
flag = 1
}
if(type == "GlobalAmine"){
	nameGraphe = paste ("Surronding of all Nref for all groups (",distance, " Å)", sep = "")
flag = 1
}

if(type == "Imidazole"){
	nameGraphe = paste ("Surronding of all Nref","\n", "Imidazoles (",distance, " Å)", sep = "")
	flag = 1
}

if(type == "AcidCarboxylic"){
	nameGraphe = paste ("Surronding of all Oref","\n", "AcidCarboxylic (",distance, " Å)", sep = "")
	flag = 1
}


if(flag == 0){
	nameGraphe = paste ("Surronding of Nref", "\n",type," amines (",distance, " Å)", sep = "")
}



data = read.table(file, header = TRUE)
print (data)

color = defColor (colnames(data))

data = data[order (as.double(rownames(data))),] 

legendHisto = rownames(data)

lenData = dim(data)[1]
large = lenData*37

if (large < 300){
	large = 300
}

maxY = 0


for (i in 1:lenData){
	sumLine = sum(data[i,])
	if (sumLine > maxY){
		maxY = sumLine
	}
}



png(filename=paste(file,".png",sep = ""),width=as.integer(large))
barplot(t(data),main=nameGraphe, ylab="Quantity", xlab = "Number of neighboring atoms", col=color, space=0.6, cex.axis=0.8, las=1, names.arg=legendHisto, cex=0.6, axes = TRUE)

legend(1,0.75 * maxY,legend=colnames(data), cex=0.8, fill=color)

dev.off()
warnings()

