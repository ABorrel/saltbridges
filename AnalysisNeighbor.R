#!/usr/bin/env Rscript

source("tool.R")




#######################
#      Main           #
#######################

args <- commandArgs(TRUE)
pathData = args[1]
pathData_distance = args[2]
study_type = args[3]




data = read.table(pathData , sep = "\t", header = TRUE)


data_distance = read.csv(pathData_distance , sep = "\t", header = FALSE)

data_distance = data_distance[,-1]

d1 = data_distance [1,-(which(data_distance[1,] == ""))]
d1_classe = data_distance [2,-(which(data_distance[1,] == ""))]
if (dim(d1)[2] == 0){
	d1 = data_distance [1,]
	d1_classe = data_distance [2,]
}

d2 = data_distance [3,-(which(data_distance[3,] == ""))]
d2_classe = data_distance [4,-(which(data_distance[3,] == ""))]
if (dim(d2)[2] == 0){
	d2 = data_distance [3,]
	d2_classe = data_distance [4,]
}


d3 = data_distance [5, -(which(data_distance[5,] == ""))]
d3_classe = data_distance [6, -(which(data_distance[5,] == ""))]
if (dim(d3)[2] == 0){
	d3 = data_distance [5,]
	d3_classe = data_distance [6,]
}

color = defColor (colnames (data))

png(filename=paste(pathData,".png",sep = ""),width=as.integer(600))
par(xpd=T, mar=par()$mar+c(0,0,0,10)) 
barplot(t(data), ylim = c(0,1), main=paste("Distribution neighbors ", study_type, sep = ""), ylab="Frequencies", xlab = "", col=color, space=0.6, cex.axis=1, las=1, cex=1, axes = TRUE)
legend("right",legend=colnames(data), fill=color,inset=c(-0.2,0))
dev.off()

# density by type
#png(filename=paste(pathData_distance ,"_type.png",sep = ""),width=as.integer(800), heigh = 1000)
#par(mfrow = c(3,1))

#plot(density(as.double(t(d1))), ylim = c(1,6))
#for (class in names(color)){
#	try (lines (density(as.double(t(d1[which(d1_classe == class)]))), col = color[class]))
#}

#plot(density(as.double(t(d2))), ylim = c(1,6))
#for (class in names(color)){
#	try (lines (density(as.double(t(d2[which(d2_classe == class)]))), col = color[class]))
#}

#plot(density(as.double(t(d3))), ylim = c(1,6))
#for (class in names(color)){
#	try (lines (density(as.double(t(d3[which(d3_classe == class)]))), col = color[class]))
#}


#dev.off()

# density distance
#png(filename=paste(pathData_distance ,".png",sep = ""),width=as.integer(800), heigh = 1000)
#par(mfrow = c(3,1))

#plot(density(as.double(t(d1))), xlim = c(1,5))
#plot(density(as.double(t(d2))), xlim = c(1,5))
#plot(density(as.double(t(d3))), xlim = c(1,5))


#dev.off()
