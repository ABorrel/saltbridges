#!/usr/bin/env Rscript


openFile = function(path, type, line){
	
	file = paste(path,"proportionType", type, sep = "")
	#print (file)
	data = read.csv(file , sep = "\t", header = FALSE)
	nbCol = dim(data)[2]
	data = data[,-nbCol]
	#print (data[line,])
	
	return (data[line,])


}


frequencyMatrix = function (data){
	nbCol = dim(data)[2]
	nbLine = dim(data)[1]

	for (i in 1:nbLine){
		sumLine = sum(data[i,])
		if(sumLine != 0){
			for(j in 1:nbCol){
				data[i,j] = data[i,j]/sumLine
			}
		}
	}

	return (data)


}

plotGlobal = function (data, distance, listType, path){

	#print (data)
	color = c("red","orange","yellow","cyan","blue","green","purple","grey")

	legendHisto = listType



	png(filename=paste(path,"porportionAllType",distance,".png",sep = ""),width=as.integer(600))
	barplot(t(data), ylim = c(0,1),main=paste("Space arround type structure\n",distance), ylab="Frequencies", xlab = "Type Structure", col=color, space=0.6, cex.axis=0.8, las=1, names.arg=listType, cex=0.6, axes = TRUE)

	legend(1,0.9,legend=c("O (COOH)", "O (Tyr, SER, THR), S (CYS)","O (H2O)", "O (main chain) Side chain ASN, GLN", "N (HIS, LYS, ARG) and Nxt", "N (main chain) ASN, GLN", "Others", "C (side chain TYR, PHE, TRP)"), cex=0.8, fill=color)

	dev.off()

}


#######################
#      Main           #
#######################

args <- commandArgs(TRUE)
nbDistance = args[1]
pathData = args[2]



####Retrieve list distance####



listDistance = c()
for(i in 0:nbDistance){
	deb = 2 + 0.5*i
	listDistance = c(listDistance,paste("D<",deb, sep = ""))

}
################################

listType = c("Primary", "Secondary", "Tertiary", "Diamine", "Guanidium","Imidazole","Pyridine", "Global")


countDistance = 0


for (distance in listDistance){
	data = NULL
	countDistance =  countDistance + 1
	for (type in listType){
		#print (paste(type, countDistance))
		data = rbind(data,openFile(pathData, type,countDistance))
	}
	data = frequencyMatrix (data)
	#print (data)
	plotGlobal(data,distance, listType, pathData)
}




warnings()
