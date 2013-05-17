#!/usr/bin/env Rscript



source ("tool.R")



###########################################
#               MAIN                      #
###########################################


args <- commandArgs(TRUE)
file = args[1]
typeStudy = args[2]

data = read.table (file, sep = "\t", header = TRUE)
data = data[,-dim(data)[2]]
print (data)

nameGroup = rownames (data)
colorGrp = defColorSubstruct(nameGroup)

print (nameGroup)
print (colorGrp)

png(filename=paste(file,".png",sep = ""),width=800)
barplot(as.matrix(data), beside=TRUE, ylim = c(0,1), col = colorGrp, main=paste("Percentage of reference with at least one ",typeStudy,"\nat selected distances", sep = ""))#, ylab= "% of references atoms", beside=TRUE, col=colorGrp, xlab = paste(typeStudy," to reference distances (Å)", sep = ""), names.arg = nameGroup, ylim = c(0,1))
legend( 2 , 1 , legend=nameGroup, bty="n", fill=colorGrp)

dev.off()
