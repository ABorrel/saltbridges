#!/usr/bin/env Rscript

source("tool.R")

#######################
#      Main           #
#######################

args <- commandArgs(TRUE)
p_filin = args[1]

d_count = read.table (p_filin, header = TRUE)
l_sub = rownames (d_count)


AFC (d_count, p_filin)

for (sub in l_sub){
	d_pie = d_count[sub,]
	pieType  (d_pie, paste(p_filin, "_", sub , sep = ""))
}

