#!/usr/bin/env Rscript


formatForMean = function (vect_in){
	vect_in = na.omit(vect_in)
	return (vect_in)
}



meanWater = function (data){
	output = NULL
	for (resolution in seq(0.5,6,0.2)){
		i_value = which (data[,2]<resolution & data[,2]>(resolution - 0.2))
		m_n_water = mean(formatForMean(data[i_value,3])) 
		
		line_add = c (resolution, m_n_water)
		output = rbind (output, line_add)
	}
	return (output)
}



##################
#     MAIN       #
##################


args <- commandArgs(TRUE)
path_file = args[1]
data = read.table (path_file, sep = "\t")



mean_data = meanWater (data)
data = na.omit(data)

png(filename=paste(path_file,"mean.png",sep = ""), 600, 600)
plot (mean_data[,1], mean_data[,2], main = "Number of H2O \n function of resolution", xlab = "Resolution Å", ylab = "Number of H2O", col = "black", type = "l")
dev.off()
