#!/usr/bin/env Rscript

library(DescTools)

source ("toolStat.R")


##########
#  MAIN  #
##########


args <- commandArgs(TRUE)
file_table1 = args[1]
file_table2 = args[2]


d1 = read.table (file_table1, sep = "\t", header = TRUE)
d2 = read.table (file_table2, sep = "\t", header = TRUE)


#print (d)
t_out = data.frame()


###################
#  check quality  #
###################

if (dim (d1) != dim (d2)){
  print ("=> ERROR Size of the data")
}


###check qulity
l_query = rownames(d1)
l_attype = colnames(d1)
df_corrected = correctionDF (dim (d1)[2])





for (query in l_query){
  sum_q1 = sum (d1[query,])
  sum_q2 = sum (d2[query,])
  
  # build table 2*2
  for (attype in l_attype){
    d1_in = c(d1[query, attype], sum (d1[query, c(l_attype[-which(l_attype == attype)])]))
    d2_in = c(d2[query, attype], sum (d2[query, c(l_attype[-which(l_attype == attype)])]))
    print ("Table 2*2")
    print (query)
    print (attype)
    print (d2_in)
    print (d1_in)
    print ("*****")
    d_in = rbind (d1_in, d2_in)
    
    
    # test choice
    #X2 test
    if (sum_q1 > 1000 && sum_q2 > 1000){
      
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

write.csv(t_out, file = paste(file_table1, "_2tablesignif.csv", sep = ""))
print (t_out)


