#!/usr/bin/env Rscript


#######################
#      Main           #
#######################

args <- commandArgs(TRUE)
p_file = args[1]

d = read.table (p_file, sep = "\t", header = TRUE)

for (sub in rownames (d)){
	png (paste(p_file, "_", sub, ".png", sep = ""), 800, 800)

	leg = NULL
	for (l in colnames (d)){
		#print (d[sub,l])
		#print (sum(d[sub,]))
		leg = append (leg, paste (l, "\n", round (d[sub,l]/sum(d[sub,])*100), "%", sep = ""))
	}

	colors = c("red", "orange", "blue", "grey")
	pie(as.double(d[sub,]), col = colors, label = leg, lwd = 10, cex = 2)
	dev.off ()
	
	svg(filename=paste(p_file, "_", sub, ".svg", sep = ""), 10, 10)
	try(pie(as.double(d[sub,]), col = colors, label = leg, cex = 1.6))
	dev.off()
}



