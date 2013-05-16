#!/usr/bin/env Rscript

source("tool.R")

openFile = function(path, type, line){
	
	file = paste(path,"proportionType", type, sep = "")
	#print (file)
	data = read.csv(file , sep = "\t", header = TRUE)
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
	color = defColor(colnames(data))

	legendHisto = listType



	png(filename=paste(path,"porportionAllType",distance,".png",sep = ""),width=as.integer(600))
	#par(mar = c(6,6,6,10))	
	par(xpd=T, mar=par()$mar+c(0,0,0,10)) 
	barplot(t(data), ylim = c(0,1),main=paste("Space arround type structure\n",distance), ylab="Frequencies", xlab = "", col=color, space=0.6, cex.axis=1, las=2, names.arg=listType, cex=1, axes = TRUE)

	legend("right",legend=colnames(data), fill=color,inset=c(-0.2,0))

	dev.off()

}


#######################
#      Main           #
#######################

args <- commandArgs(TRUE)
pathData = args[1]


####Retrieve list distance####

listType = c("Primary", "Secondary", "Tertiary", "Diamine", "Guanidium","Imidazole","Pyridine", "AcidCarboxylic", "Global")



list_distance = rownames(read.table(paste(pathData, "proportionType", listType[1], sep = "")))
i_line_distance = 0
print (list_distance)
for (distance in list_distance){
	print (distance)
	data = NULL
	i_line_distance =  i_line_distance + 1
	for (type in listType){
		#print (paste(type, i_line_distance))
		data = rbind(data,openFile(pathData, type,i_line_distance))
		#print (data)
	}
	data = frequencyMatrix (data)
	#print (data)
	plotGlobal(data,distance, listType, pathData)
}




warnings()
