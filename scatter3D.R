#!/usr/bin/env Rscript

library("vrmlgen")
require(lattice)
library(scatterplot3d)
source ("tool.R")


#############
#  MAIN     #
#############


args <- commandArgs(TRUE)
p_file = args[1]

d = read.table (p_file, sep = "\t")



l_type = unique(d[,4])
#print (l_type)


for (t in l_type){
	d_temp1 = d[which(d[,4] == t),]
	d_temp = rbind (d_temp1,d[which(d[,4] == "REF"),] )

	# VRML
	if (dim (d_temp1)[1] > 500){	
		print ("1111")
		s_500 = c(sample (dim (d_temp1)[1])[1:450],which(d[,4] == "REF"))
		d_vrml = d_temp[s_500,]
	}else {
		d_vrml = d_temp
	}
	
	print (dim(d_vrml))	
	
	#cloud3d(d_vrml[,1:3], vrml_showdensity = TRUE, labels = paste ("type", as.numeric(d_vrml[,4])), filename = paste(p_file, t, ".wrl", sep = ""))

	# GIF
	col_group = defColor (d_temp[,4])
	png(file=paste(p_file, "3DPlot%03d.png", sep = ""), width=600, height=600)

	for (i in seq(0, 350 ,10)){
		print(cloud(d_temp[,1]~d_temp[,2]*d_temp[,3], col = col_group, pch = 16, screen = list(z = i, x = -80)))
	}
	dev.off()

	# convert pngs to one gif using ImageMagick
	system(paste("convert -delay 50 ", p_file,"*.png ", p_file, t, ".gif", sep = ""))

}



