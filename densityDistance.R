#!/usr/bin/env Rscript
source("tool.R")

#######################
#      Main           #
#######################

args <- commandArgs(TRUE)
p_density = args[1]

d = read.csv(p_density , sep = "\t", header = FALSE)

l_class = unique (d[,2])
l_color = defColor (l_class)
print (l_class)
print (l_color)

png(filename=paste(p_density ,".png",sep = ""),width=as.integer(1000), heigh = 600)

dist_temp = d[which (d[,2] == l_class[1]),1]
plot (density (dist_temp), lwd = 3, xlim = c(2,6), ylim = c(0,1), col = l_color[l_class[1]], main ="")
l_class = l_class[-1]

for (x in seq (2,6,0.25)){
	segments (x,0,x,1,lwd = 2)
}
for (y in seq (0,1,0.25)){
	segments (2,y,6,y,lwd = 2)
}

for (class in l_class){

	dist_temp = d[which (d[,2] == class),1]
	lines (density(dist_temp), col = l_color[class], lwd = 3)

}
legend("topleft", col=l_color, legend = names (l_color), pch=rep (16,length (l_color)),  cex = 1.5)
dev.off()
