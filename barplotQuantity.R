#!/usr/bin/env Rscript


args <- commandArgs(TRUE)
file = args[1]
distance = args[2]
type = args[3]
struct = args[4]

print (type)
if (type == "Atoms"){
	nameX = "Atoms"
}else if (type == "Ligands"){
	nameX = "Ligands"

}else if (type == "Residues"){
	nameX = "Residues"
}



if (type == "Residues"){
	
	color = c("#A8A8A8","#FF3300", "#0066FF") 
	data = read.table(file)
	

	lenData = dim(data)[1]
	large = dim(data)[1]*37

	line0 = rep(0,21)
	data = cbind(data,line0)
	
	for(i in 1:lenData){
		if (data[i,1]=="HOH"){
		data[i,4] = data[i,2]
		data=data[order(data[,2]+data[,3],decreasing = T),]
		
		}
	} 

	for(i in 1:lenData){
		if (data[i,1]=="HOH"){
		data[i,2] = 0
		}
	} 
	

	if (large < 300){
		large = 300
	}
	
	
	legendHisto = data[,1]
	data = data[,-1]

	
	png(filename=paste(file,".png",sep = ""),width=as.integer(large))
	barplot(t(data),main=paste("Environment in amino acid around Nref of ",struct," amine","\n",distance,"Å"), ylab="Quantity", xlab = nameX, col=color, space=0.6, cex.axis=0.8, las=1, names.arg=legendHisto, cex=0.6, axes = TRUE)

	legend(20,0.5 * max(data),legend=c("Main chain", "Side chain", "Water"), cex=0.8, fill=color)
	dev.off()

}else{
	print (file)
	data = read.table(file)
	data=data[order(data[,2],decreasing = T),]

	sumRetrieve = sum(data[,2])


	large = dim(data)[1]*40

	if (large < 300){
		large = 300
	}

	png(filename=paste(file,".png",sep = ""), width=as.integer(large))
	barplot(data[,2], names.arg=data[,1], main=paste("Distribution of residue around Nref","\n",distance,"Å"), xlab=nameX, ylab="Quantity", axes=TRUE, cex.axis=1, cex.lab=1, cex.main=1.5, cex.names = 0.6, col = "green", las=1,space = 0.7, cex = 0.6)

	dev.off()

}

warnings()
