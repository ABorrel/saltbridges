#!/usr/bin/env Rscript
source("tool.R")



pieType = function (d, path_out){
	
	#print (data)
	colors = defColor(names(d))
	par (lwd = 1000)
	png(filename=paste(path_out,".png",sep = ""),4000, 4000)
	try(pie(as.double(d), col = colors, label = "", lwd = 10))
	dev.off()

	svg(filename=paste(path_out,".svg",sep = ""))
	try(pie(as.double(d), col = colors, label = names(d)))
	dev.off()
}





#######################
#      Main           #
#######################

args <- commandArgs(TRUE)
pathData = args[1]
nb_type = as.integer(args[2])



png(filename=paste(pathData ,".png",sep = ""), width=as.integer(800), heigh = 1000)
par(mfrow = c(nb_type/2+1,2))


for (i in 1:nb_type){
	d = read.table(paste(pathData, "_",i, sep = ""), sep = "\t", header = TRUE)
	
	indiv =  (sum(d))
	color = defColor (colnames (d))
	barplot(t(d), ylab= "Number of occurences", col = color, xlab = "distance (Å)", names.arg = rownames (d),space = 0.5, main =paste(i, indiv, sep = "   "))

	

	if (i == 2){
		legend("left", legend=colnames(d), fill=color)

	}

}

dev.off()




for (i in 1:nb_type){
	d = read.table(paste(pathData, "_",i, sep = ""), sep = "\t", header = TRUE)
	d_global = read.table(paste("/home/borrel/saltBridgesProject/result/PDB50/3.00_angle_morecomplexe/neigbhor/barplot_global_",i, sep = ""), sep = "\t", header = TRUE)

	
	data_pie = apply(d,2,"sum")
	data_global_pie = apply(d_global,2,"sum")
	d_x2 = rbind (data_pie, data_global_pie)
	print (chisq.test(d_x2)$p.value)
	pieType  (data_pie, paste(pathData, "_pie_", i, sep = ""))
	
}




