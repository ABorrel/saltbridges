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
data = cbind (data, data[,5]/data[,3])
#data = data[-(which(data[,6] == "Inf")),]
#data = data[-(which(data[,7] == "Inf")),]
#data = data[-(which(data[,2] >= 6)),]


# global
svg (paste(path_file, "RXfunctionNbH2O.svg", sep = ""), 10, 10)
par (mar = c(5,5,5,5))
plot(data[,2], data[,6], xlim = c(0.5,4), ylim = c(0,0.4), pch = 19, xlab = "Resolution (Å)", ylab = "Mean of water molecules by residue", cex = 0.6, cex.lab = 2)

m = c(mean (data[which (data[,2] < 0.5),6]), mean (data[which (data[,2] < 1 & data[,2] > 0.5),6]), mean (data[which (data[,2] < 1.5 & data[,2] > 1),6]), mean (data[which (data[,2] < 2 & data[,2] > 1.5),6]), mean (data[which (data[,2] < 2.5 & data[,2] > 2),6]), mean (data[which (data[,2] < 3 & data[,2] > 2.5),6]), mean (data[which (data[,2] < 3.5 & data[,2] > 3.0),6]), mean (data[which (data[,2] < 4 & data[,2] > 3.5),6]))
print (m)
points (c(0.5,1,1.5,2,2.5,3,3.5,4), m, lwd = 4, type = "l", col = "red")
dev.off()

mean_data = meanWater (data)
png(filename=paste(path_file,"mean.png",sep = ""), 600, 1200)
par(mfrow = c(2,1))
plot (mean_data[,1], mean_data[,2], xlim = c(0,6) , ylim = c(0,5), main = "Ratio number of proton by exposed residues \n function of resolution", xlab = "Resolution Å", ylab = "Number of proton by exposed residues", col = "black", type = "l")
plot (mean_data[,1], mean_data[,3], xlim = c(0,6) , ylim = c(0,5), main = "Ratio number of proton by residues \n function of resolution", xlab = "Resolution Å", ylab = "Number of proton by residues", col = "black", type = "l")
dev.off()

if ((sum (data[,3])) != 0){
	png(filename=paste(path_file,"exposed.png",sep = ""), 600, 1200)
	par(mfrow = c(2,1))
	hist(data[,2], main = "Distribution function of resolution", xlab = "Resolution Å", ylab = "Number of PDB", col = "gray", xlim =c(0,6))
	plot (data[,2], data[,6], xlim = c(0,6) , ylim = c(0,10), main = "Ratio number of proton by exposed residues \n function of resolution", xlab = "Resolution Å", ylab = "Number of proton by exposed residues", col = "gray")
	dev.off()
}

png(filename=paste(path_file,".png",sep = ""), 600, 1200)
par(mfrow = c(2,1))
hist(as.vector(data[,2]), main = "Distribution function of resolution", xlab = "Resolution Å", ylab = "Number of PDB", col = "gray")
plot(data[,2], data[,7], xlim = c(0,6),ylim = c(0,10), main = "Ratio number of proton by residues\n function of resolution", xlab = "Resolution Å", ylab = "Number of proton by exposed residues", col = "gray")
dev.off()



