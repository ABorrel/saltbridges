

for (i in 1:nb_type){
	d = read.table(paste(pathData, "_",i, sep = ""), sep = "\t", header = TRUE)
	d_global = read.table(paste("/home/borrel/saltBridgesProject/result/PDB50/3.00_angle_morecomplexe/neigbhor/barplot_global_",i, sep = ""), sep = "\t", header = TRUE)

	
	data_pie = apply(d,2,"sum")
	data_global_pie = apply(d_global,2,"sum")
	d_x2 = rbind (data_pie, data_global_pie)
	print (chisq.test(d_x2)$p.value)
	pieType  (data_pie, paste(pathData, "_pie_", i, sep = ""))
	
}



