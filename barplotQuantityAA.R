#!/usr/bin/env Rscript


args <- commandArgs(TRUE)
file = args[1]
distance = args[2]
aa = args[3]
distanceMax = args[4]

nameX = paste("Atom of ",aa , sep = "")

data = read.table(file)
data = data[order(data[,2]+data[,3],decreasing = T),]

nameX = data[,1]
data = data[,-1]


sumRetrieve = sum(data[,1]) + sum(data[,2])


lenData = dim(data)[1]
for (i in 1:lenData){
	data[i,2] = data[i,2]
	data[i,1] = data[i,1]
}


large = dim(data)[1]*37

if (large < 300){
	large = 300
}

png(filename=paste(file,".png",sep = ""),width=as.integer(large))


barplot(t(data),main=paste("Distribution of atoms","\n","Amino Acid:",aa,sep = " "), ylab="Number of occurences",xlab = "Atoms (reference name in PDB)", col=heat.colors(2), space=0.6, cex.axis=0.8, las=1, names.arg=nameX, cex=0.6, axes = TRUE)

legend(lenData-3,0.5*max(data),legend=c("1.7 Å < distance < 3.5 Å", paste("3.5 Å < distance < ",distanceMax,"Å", sep = "")), cex=0.7, fill=heat.colors(2))

dev.off()

warnings()
