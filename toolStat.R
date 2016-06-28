#!/usr/bin/env Rscript




##########################
#  Bufferoni correction  #  
##########################


correctionDF = function (df){
  i = 0
  df = df - 1
  while (df > 0) {
    i = i + df
    df = df - 1
  }
  return (i)  
  
}


###############
#   Signif    #
###############



signifPvalueCorrected = function (a, n){
  if (a < (0.001/n)){
    return ("***")
  }
  else if (a < (0.01/n)){
    return ("**")
  }
  else if (a < (0.05/n)){
    return ("*")
  }
  else {
    return ("-")
  }	
}


