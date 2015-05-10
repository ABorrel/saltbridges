#!/usr/bin/env Rscript


source("tool.R")




#######################
#        MAIN         #
#######################

args <- commandArgs(TRUE)
file = args[1]
data = read.csv(file, sep = "\t", header = FALSE)


nbLine = dim(data)[1]
nb_col = dim(data)[2]

# case with several angle -> deviation with median position
# print (nb_col)
if (nb_col > 3){
	data = deviationAngle(data)
}

nbLine = dim(data)[1]
nb_col = dim(data)[2]

# retrieve type of neigbors
l_element = unique(data[,nb_col])
color_element  = defColor(l_element)

png(filename=paste(file, "_angleDistribution.png", sep = ""))
hist(as.double(data[,2]), xlab ="Angles or deviation", ylab = "Number of occurences", right=F, main=paste("Distribution of deviation or angles", sep = ""),las=1, freq=T, col = "#D8D8D8")
dev.off()
