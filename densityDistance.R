#!/usr/bin/env Rscript
source("tool.R")

#######################
#      Main           #
#######################

args <- commandArgs(TRUE)
p_density = args[1]

d = read.csv(p_density , sep = "\t", header = FALSE)

# define classes and color
l_class = unique (d[,2])
l_color = defColor (l_class)

# remove row with less than 2 element
for (class in l_class){
	i_temp = which (d[,2] == class)
	if (length (i_temp) < 2){
		d = d[-i_temp,]
	}
}

# grade the plot
all_y = NULL
all_x = NULL
for (class in l_class){
	all_y = append (all_y, density(d[which (d[,2] == class),1])$y)
	all_x = append (all_x, density(d[which (d[,2] == class),1])$x)
}

max_x = max (all_x)
max_y = max (all_y)

svg (filename=paste(p_density ,".svg",sep = ""), bg = "transparent",16 , 12)

plot (density (d[which (d[,2] == class),1]), main ="", xlim = c(1, 6), ylim = c(0, max_y), type = "n", xlab = "", cex.axis = 2, ylab = "")


for (class in l_class){
	dist_temp = d[which (d[,2] == class),1]
	lines (density(dist_temp), col = l_color[class], lwd = 3.5)
}


# grid
# horizontal
y_grid = seq(0, max_y, 0.2)
for (y in y_grid){
	segments (1, y, max_x, y, lty = 2, col = "black", lwd = 1.5)
}

# vertical -> need improvment with a loop
x_grid = seq(1, max_x ,0.5)
for (x in x_grid){
	segments (x, 0, x, max_y, lty = 2, col = "black", lwd = 1.5)
}


legend("topleft", col=l_color, legend = names (l_color), pch=rep (16,length (l_color)),  cex = 2)
dev.off()
