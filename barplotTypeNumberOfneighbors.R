#!/usr/bin/env Rscript

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

if(flag == 0){
	nameGraphe = paste ("Surronding of Nref", "\n",type," amines (",distance, " Å)", sep = "")
}

color = c("red","orange","yellow","cyan","blue","green","purple","grey")




data = read.table(file)
data = data[order(data[,1],decreasing = F),]

legendHisto = data[,1]
data = data[,-1]

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

legend(1,0.75 * maxY,legend=c("O (COOH)", "O (Tyr, SER, THR), S (CYS)","O (H2O)", "O (main chain) Side chain ASN, GLN", "N (HIS, LYS, ARG) and Nxt", "N (main chain) ASN, GLN", "C (side chain TYR, PHE, TRP)", "Others"), cex=0.8, fill=color)

dev.off()

warnings()
