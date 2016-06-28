#!/usr/bin/env Rscript

library(DescTools)

source ("toolStat.R")

##################
#     MAIN       #
##################


args <- commandArgs(TRUE)
file_table = args[1]

d = read.table (file_table, sep = "\t", header = TRUE)
#print (d)
t_out = data.frame()

ref = d["global",]
d = d[-c(which(rownames(d) == "global")),]

l_query = rownames(d)
l_attype = colnames(d)
df_corrected = correctionDF (dim (d)[2])

for (query in l_query){
  sum_q = sum (d[query,])
  #print (sum_q)
  
  # build table 2*2
  for (attype in l_attype){
    d_in = c(d[query, attype], sum (d[query, c(l_attype[-which(l_attype == attype)])]))
    ref_in = c(ref[, attype], sum (ref[, c(l_attype[-which(l_attype == attype)])]))
    #print ("Table 2*2")
    #print (query)
    #print (attype)
    #print (d_in)
    #print ("*****")
    #print (ref_in)
    d_in = rbind (d_in, ref_in)
    
    
    # test choice
    #X2 test
    if (sum_q > 1000){
      
      #print (chisq.test(d_in))
      pval = chisq.test(d_in)$p.val
      pval = signif(pval, digits = 2)
      #print (pval)
      t_out[query, attype] = paste(pval, signifPvalueCorrected(pval,df_corrected), sep = " ")
      
    }
    #Gtest
    else {
      # small pop
      #print (GTest(d_in))
      #print (chisq.test(d_in))
      pval = GTest(d_in)$p.val
      pval = signif(pval, digits = 2)
      t_out[query, attype] = paste(pval, signifPvalueCorrected(pval,df_corrected), sep = " ")
      
    }
  }
}

write.csv(t_out, file = paste(file_table, "_signif.csv", sep = ""))
print (t_out)





