#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
file = args[1]
type = args[2]
distance = args[3]
flag = 0

if(type == "Global"){
	nameGraphe = paste ("Surronding of all atoms", sep = "")
flag = 1
}
if(type == "GlobalAmine"){
	nameGraphe = paste ("Surronding of all Nref for all groups", sep = "")
flag = 1
}

if(type == "Imidazole"){
	nameGraphe = paste ("Surronding of all Nref","\n", "Imidazoles", sep = "")
flag = 1
}

if(flag == 0){
	nameGraphe = paste ("Surronding of Nref", "\n",type," amines", sep = "")
}

if(type == "AcidCarboxylic"){
	nameGraphe = paste ("Surronding of all Oref","\n", "AcidCarboxylic (",distance, " Å)", sep = "")
	flag = 1
}



color = c("red","orange","yellow","cyan","blue","green","purple", "grey")



data = read.csv(file, sep = "\t", header = FALSE)
print (data)



nbDistance = dim(data)[1] -1

distance = c()

for(i in 0:nbDistance){
	deb = 2 + 0.5*i
	distance = c(distance,paste(deb,"<D", sep = ""))

}
print (distance)



lenData = dim(data)[1]
large = lenData*37

if (large < 300){
	large = 300
}

maxY = 0


for (i in 1:lenData){
	sumLine = data[i,1] + data[i,2] + data[i,3] + data[i,4] + data[i,5]
	if (sumLine > maxY){
		maxY = sumLine
	}
}



png(filename=paste(file,".png",sep = ""), width=as.integer(large))
barplot(t(data), main=nameGraphe, ylab="Quantity", xlab = "Distance", col=color, space=0.6, cex.axis=0.9, las=1, cex=0.6, axes = TRUE, names.arg=distance)

legend(2,0.75 * maxY,legend=c("OxAcid", "Amphiprotic","H2O", "OxAccept", "Nbasic", "Ndonnor", "Carbon aromatic","Others"), cex=0.8, fill=color)

dev.off()

warnings()
