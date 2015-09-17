#!/usr/bin/env Rscript



formatForMean = function (vect_in){
	vect_in = na.omit(vect_in)
	return (vect_in)
}


meanWater = function (data){


	output = NULL

	for (resolution in seq(0.5,6,0.2)){

		
		i_value = which (data[,2]<resolution & data[,2]>(resolution - 0.2))
		m_nb_water_exposed = mean(formatForMean(data[i_value,6])) 
		#print (data[i_value,6])
		m_nb_water = mean(formatForMean(data[i_value,7]))
		
		line_add = c (resolution, m_nb_water_exposed, m_nb_water)
		output = rbind (output, line_add)
	}
	return (output)
}




##################
#     MAIN       #
##################


args <- commandArgs(TRUE)
path_file = args[1]
data = read.table (path_file, sep = "\t", header = TRUE)

data = cbind (data, data[,5]/data[,4])
data = data[-(which(data[,6] == "Inf")),]
#data = cbind (data, data[,5]/data[,3]) -> exposed residues not used now
#data = data[-(which(data[,6] == "Inf")),]
#data = data[-(which(data[,7] == "Inf")),]
#data = data[-(which(data[,2] >= 6)),]

# global
svg (paste(path_file, "RXVSH2O.svg", sep = ""), 10, 10)
par (mar = c(5,5,2,2))
plot(data[,2], data[,6], xlim = c(0.5,4), ylim = c(0,0.4), pch = 19, xlab = "Resolution (Å)", ylab = "Mean of water molecules by amino acid", cex = 0.6, cex.lab = 2)

m = c(mean (data[which (data[,2] < 0.5),6]), mean (data[which (data[,2] < 1 & data[,2] > 0.5),6]), mean (data[which (data[,2] < 1.5 & data[,2] > 1),6]), mean (data[which (data[,2] < 2 & data[,2] > 1.5),6]), mean (data[which (data[,2] < 2.5 & data[,2] > 2),6]), mean (data[which (data[,2] < 3 & data[,2] > 2.5),6]), mean (data[which (data[,2] < 3.5 & data[,2] > 3.0),6]), mean (data[which (data[,2] < 4 & data[,2] > 3.5),6]))
print (m)
points (c(0.5,1,1.5,2,2.5,3,3.5,4), m, lwd = 4, type = "l", col = "red")
dev.off()

# global
png (paste(path_file, "RXVSH2O.png", sep = ""), 800, 800)
par (mar = c(5,5,2,2))
plot(data[,2], data[,6], xlim = c(0.5,4), ylim = c(0,0.4), pch = 19, xlab = "Resolution (Å)", ylab = "Mean of water molecules by amino acid", cex = 0.6, cex.lab = 2)

m = c(mean (data[which (data[,2] < 0.5),6]), mean (data[which (data[,2] < 1 & data[,2] > 0.5),6]), mean (data[which (data[,2] < 1.5 & data[,2] > 1),6]), mean (data[which (data[,2] < 2 & data[,2] > 1.5),6]), mean (data[which (data[,2] < 2.5 & data[,2] > 2),6]), mean (data[which (data[,2] < 3 & data[,2] > 2.5),6]), mean (data[which (data[,2] < 3.5 & data[,2] > 3.0),6]), mean (data[which (data[,2] < 4 & data[,2] > 3.5),6]))
print (m)
points (c(0.25,0.75,1.25,1.75,2.25,2.75,3.25,3.75), m, lwd = 4, type = "l", col = "red")
dev.off()
