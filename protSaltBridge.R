#!/usr/bin/env Rscript



MultiplotDens = function (v_in, number_run, color_lines){
  # graphic window open
  nb_elem = length  (v_in)
  
  if (nb_elem > 100000){
    div_thresold = 100000
  }
  if (nb_elem > 50000){
    div_thresold = 50000
  }else {
    div_thresold = round(nb_elem/5, 0)
  }
  
  for (i in 1:number_run){
    lines (density (v_in[sample (length (v_in))[1:div_thresold]]), lwd = 2, col = color_lines)
  }
}




##### MAIN #####
################



args <- commandArgs(TRUE)
p_file = args[1]


d = read.table(p_file, sep = "\t", header = TRUE)
print ("#### quantity of AA considered ####")
print (dim (d))
d_filter = d[-which(d[,3] == 0),]
print ("#### quantity of AA in salt bridge ####")
print (dim(d_filter))

div_thresold = 80000


###### summary ######

print ("#### % of salt bridges ####")
percent_salt_bridges = sum(d[,3])/dim(d)[1]
print (percent_salt_bridges)
print ("")
print ("#### % of salt bridges by aa####")
unique_AA = unique(d_filter[,1])
for (aa in unique_AA){
  print (paste("=> ", aa, sep = ""))
  print (dim (d_filter[which (d_filter[,1] == aa),])[1]/dim(d_filter)[1])
}

# density distance
v_dist = as.numeric(as.character(d_filter[,2]))

# hist
svg(filename=paste(p_file,"_histdist.svg",sep = ""),  bg = "transparent", 16, 12)
par(mar = c(5,5,1,1), mfrow = c(2,2))
hist (v_dist[sample (length (v_dist))[1:div_thresold]] , main = "", xlim = c(1,4.5), cex.axis = 2)
hist (v_dist[sample (length (v_dist))[1:div_thresold]] , main = "", xlim = c(1,4.5), cex.axis = 2)
hist (v_dist[sample (length (v_dist))[1:div_thresold]] , main = "", xlim = c(1,4.5), cex.axis = 2)
hist (v_dist[sample (length (v_dist))[1:div_thresold]] , main = "", xlim = c(1,4.5), cex.axis = 2)
dev.off()

#Density

max_y = max(density(v_dist)$y)
svg(filename=paste(p_file,"_densitydist.svg",sep = ""),  bg = "transparent", 16, 12)
par(mar = c(5,5,1,1), mfrow = c(1,1))
plot (density (v_dist), main = "", xlim = c(1,4.5), cex.axis = 2, type = NULL)

MultiplotDens (v_dist, 200, 1)

# grid
# horizontal
y_grid = seq(0, max_y, 0.2)
for (y in y_grid){
  segments (1, y, 4.5, y, lty = 2, col = "black", lwd = 1.5)
}

# vertical -> need improvment with a loop
x_grid = seq(1, 4.5 ,0.5)
for (x in x_grid){
  segments (x, 0, x, max_y, lty = 2, col = "black", lwd = 1.5)
}

dev.off()



# density by AA

l_dist = list()
max_y_sp = 0
l_aa_plot = NULL
for (aa in unique_AA){
  v_dis = as.numeric(as.character(d_filter[which(d_filter[,1] == aa),2]))
  vd_dist = density(v_dis)
  if (max(vd_dist$y) > max_y_sp){
    max_y_sp = max(vd_dist$y)
    d_max = v_dis
    aa_max = aa
  }
}

l_aa_plot = as.character(unique_AA[-which(aa_max == unique_AA)])
print (l_aa_plot)
print (aa_max)

svg(filename=paste(p_file,"_densitybyaa.svg",sep = ""),  bg = "transparent", 16, 12)
par(mar = c(5,5,1,1), mfrow = c(1,1))
plot (density(d_max), main = "", xlim = c(1,4.5), cex.axis = 2, type = NULL)

MultiplotDens (d_max, 200, 1)

i_col = 1
for (aa in l_aa_plot){
  i_col = i_col + 1
  v_plot = as.numeric(as.character(d_filter[which(d_filter[,1] == aa),2]))
  MultiplotDens (v_plot, 200, i_col)
  
}

# legend
legend("topleft", col = seq (i_col), legend = c(aa_max, l_aa_plot), pch=16,  cex = 2)


# grid
# horizontal
y_grid = seq(0, max_y, 0.2)
for (y in y_grid){
  segments (1, y, 4.5, y, lty = 2, col = "black", lwd = 1.5)
}

# vertical -> need improvment with a loop
x_grid = seq(1, 4.5 ,0.5)
for (x in x_grid){
  segments (x, 0, x, max_y, lty = 2, col = "black", lwd = 1.5)
}

dev.off()



