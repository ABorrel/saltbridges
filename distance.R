#!/usr/bin/env Rscript




grapheCN = function(file, base){
	data = read.table(file)

	#print (data) 	

	brk = seq(0,max(data[,1]),0.01)

	png(filename=paste(file,".png",sep = ""))
	hist(data[,1], xlab ="Length of bond C,N (Å)", ylab = "Number of occurences", xlim=c(1,2), breaks = brk, main=paste("Distribution of lengths for C,N bonds","\n", base," dataset",sep = ""),las=1, freq=T)
	dev.off()

}



graphePolar = function(file, base){


	data = read.table(file)
	data = data[!is.na(data)]
	brk = seq(0,10,0.05)
	print (brk)

	png(filename=paste(file, ".png", sep = ""))

	hist(data, xlab ="orthogonal distances (Å)", ylab = "Number of occurences", breaks=brk, xlim=c(0,3), right=F, main=paste("Distribution of distances to orthogonal plane","\n", "tertiary amines, ",base," dataset",sep = ""),las=1, freq=T)

	dev.off()


}


###############main#################


args <- commandArgs(TRUE)
file = args[1]
type = args[2]
base = args[3]


print (file)
print (type)

if (type == "CN"){
	grapheCN(file, base)
}
if(type == "coplar"){
	graphePolar(file,base)

}
