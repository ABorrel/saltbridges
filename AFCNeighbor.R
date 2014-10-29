#!/usr/bin/env Rscript

source("tool.R")

#######################
#      Main           #
#######################

args <- commandArgs(TRUE)
p_filin = args[1]

d_count = read.table (p_filin, header = TRUE)
i_null = which(apply (d_count, 2, sum) == 0)

# remove colum with only 0
if (is.integer0 (i_null) == FALSE){
	d_count = d_count[,-i_null]
}


l_sub = rownames (d_count)


AFC (d_count, p_filin)

for (sub in l_sub){
	d_pie = d_count[sub,]
	pieType  (d_pie, paste(p_filin, "_", sub , sep = ""))
}

