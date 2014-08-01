#!/usr/bin/env Rscript

source("tool.R")
source("AFC.R")





#######################
#      Main           #
#######################

args <- commandArgs(TRUE)
pathDir = args[1]
nb_neighbor = as.integer(args[2])


listType = c("Primary", "Secondary", "Tertiary", "Guanidium", "Imidazole", "AcidCarboxylic", "global")

for (n in seq (7)){
	d_count = NULL
	for (struct_type in listType){
		p_data = paste (pathDir, "neighbor_count_", struct_type, sep = "")
		print (p_data)
		d = read.table (p_data, header = TRUE)
		d_count = rbind (d_count, d[n,])
	}
	rownames (d_count) = listType
	AFC (d_count, paste (pathDir,"neighbor_count_", n, sep = "" ))
}





