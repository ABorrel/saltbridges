#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
file = args[2]
type = args[1]
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

if(type == "IMD"){
	nameGraphe = paste ("Surronding of all Nref","\n", "IMD (",distance, " Å)", sep = "")
flag = 1
}

if(flag == 0){
	nameGraphe = paste ("Surronding of Nref", "\n",type," amines (",distance, " Å)", sep = "")
}

if(type == "COO"){
	nameGraphe = paste ("Surronding of all Oref","\n", "COO (",distance, " Å)", sep = "")
	flag = 1
}


color = c("#B8B8B8", "#FF3300", "#66CCFF", "#FFFF66", "#CC99CC","#CCFF99")




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
	sumLine = data[i,1] + data[i,2] + data[i,3] + data[i,4] + data[i,5]
	if (sumLine > maxY){
		maxY = sumLine
	}
}



png(filename=paste(file,".png",sep = ""),width=as.integer(large))
barplot(t(data), ylim = c(0,maxY),main=nameGraphe, ylab="Quantity", xlab = "Number of neighboring atoms", col=color, space=0.6, cex.axis=0.8, las=1, names.arg=legendHisto, cex=0.6, axes = TRUE)

legend(lenData-1,0.75 * maxY,legend=c("C", "O", "N", "S", "others"), cex=0.8, fill=color)

dev.off()

warnings()
