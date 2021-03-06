#!/usr/bin/env Rscript




grapheBond = function(file, type_bond){
	d = read.table(file)
	brk = 20
	svg(filename=paste(file,".svg",sep = ""), 10, 8)
	par (mar = c(5,5,2,2))
	hist(d[,1], xlab = paste (type_bond, "bonds (Å)", sep = ""), ylab = "Number of occurences", breaks = 100, main="", las=1, cex = 0.6, cex.lab = 2, xlim = c(0,max (d[,1])), col = "grey")
	dev.off()

}



graphePolar = function(file){


	d = read.table(file)

	svg(filename=paste(file, ".svg", sep = ""), 10, 8)
	par (mar = c(5,5,2,2))
	hist(d[,1], xlab ="Distances to plane (Å)", ylab = "Number of occurences", breaks = 100, main="", las=1, cex = 0.6, cex.lab = 2, xlim = c(0,max (d[,1])), col = "grey")
	dev.off()
}



grapheOther = function(file, type){


	data = read.table(file)
	data = data[!is.na(data)]
	brk = seq(0,max(data)+0.01,0.010)

	png(filename=paste(file, ".png", sep = ""))

	hist(data, breaks=brk, right=F, main=type ,las=1, freq=T)

	dev.off()
}



###############main#################


args <- commandArgs(TRUE)
file = args[1]
type = args[2]


print (file)
print (type)


if(type == "coplar"){
	graphePolar(file)
}
if (type == "CN" | type == "CO" | type == "CC"){
  grapheBond(file, type)
}else{
	grapheOther (file, type)
}


