#!/usr/bin/env Rscript


# with update 12-08-2014 ot used
calculmatrix = function(inMatrix){

	nbCol = dim(inMatrix)[2]
	nbLine = dim(inMatrix)[1]
	for (i in seq(1,nbLine)){
		for (j in seq(2,nbCol)){
			temp = j - 1
			while(temp >= 2){
				inMatrix[i,j] =  inMatrix[i,j] - inMatrix[i,temp]
				temp = temp -1 
			}
		}
	}
	return (inMatrix)
}


##################
#     MAIN       #
##################


args <- commandArgs(TRUE)
file = args[1]
structure = args[2]

if (structure == "all"){
	structure = "all atoms"
}

colorCharge = "red"
colorPolar = "#FF6600"
colorAr = "#9900CC"
colorOther = "#808080"
colorApolar =  "#66FF33"
colorCharge2 = "#000099"
colorWater = "#33FFFF"
colorTyr = "#FFFF00"
AA = c("HOH","ASP","GLU","THR","SER","ASN","GLN","TYR","HIS","LYS","ARG","PHE","TRP","ALA","ILE","LEU","MET","VAL","CYS","GLY","PRO")

data = read.table(file, header = FALSE,, sep = "\t")

#data = calculmatrix (data)

colorV = NULL
colorL = NULL
for(i in 1:21){
	if(data[i,1] == "GLU"){
		colorV = append(colorV,colorCharge)
	}
	if(data[i,1] == "ASP"){
		colorV = append(colorV,colorCharge)
	}
	if(data[i,1] == "ARG"){
		colorV = append(colorV, colorCharge2)
	}
	if(data[i,1] == "THR"){
		colorV = append(colorV,colorPolar)
	}
	if(data[i,1] == "SER"){
		colorV = append(colorV, colorPolar)
	}
	if(data[i,1] == "CYS"){
		colorV = append(colorV,colorApolar)
	}
	if(data[i,1] == "ASN"){
		colorV = append(colorV,colorPolar)
	}
	if(data[i,1] == "LYS"){
		colorV = append(colorV,colorCharge2)
	}
	if(data[i,1] == "GLN"){
		colorV = append(colorV,colorPolar)
	}
	if(data[i,1] == "PHE"){
		colorV = append(colorV,colorAr)
	}
	if(data[i,1] == "HIS"){
		colorV = append(colorV,colorTyr)
	}
	if(data[i,1] == "TRP"){
		colorV = append(colorV,colorAr)
	}
	if(data[i,1] == "TYR"){
		colorV = append(colorV,colorTyr)
	}
	if(data[i,1] == "ALA"){
		colorV = append(colorV,colorApolar)
	}
	if(data[i,1] == "GLY"){
		colorV = append(colorV,colorOther)
	}
	if(data[i,1] == "ILE"){
		colorV = append(colorV,colorApolar)
	}
	if(data[i,1] == "LEU"){
		colorV = append(colorV,colorApolar)
	}
	if(data[i,1] == "MET"){
		colorV = append(colorV,colorApolar)
	}
	if(data[i,1] == "PRO"){
		colorV = append(colorV,colorOther)
	}
	if(data[i,1] == "VAL"){
		colorV = append(colorV,colorApolar)
	}
	if(data[i,1] == "HOH"){
		colorV = append(colorV,colorWater)
	}
}


for(i in 1:21){
	if(AA[i] == "GLU"){
		colorL = append(colorL,colorCharge)
	}
	if(AA[i] == "ASP"){
		colorL = append(colorL,colorCharge)
	}
	if(AA[i] == "ARG"){
		colorL = append(colorL, colorCharge2)
	}
	if(AA[i] == "THR"){
		colorL = append(colorL,colorPolar)
	}
	if(AA[i] == "SER"){
		colorL = append(colorL, colorPolar)
	}
	if(AA[i] == "CYS"){
		colorL = append(colorL,colorApolar)
	}
	if(AA[i] == "ASN"){
		colorL = append(colorL,colorPolar)
	}
	if(AA[i] == "LYS"){
		colorL = append(colorL,colorCharge2)
	}
	if(AA[i] == "GLN"){
		colorL = append(colorL,colorPolar)
	}
	if(AA[i] == "PHE"){
		colorL = append(colorL,colorAr)
	}
	if(AA[i] == "HIS"){
		colorL = append(colorL,colorTyr)
	}
	if(AA[i] == "TRP"){
		colorL = append(colorL,colorAr)
	}
	if(AA[i] == "TYR"){
		colorL = append(colorL,colorTyr)
	}
	if(AA[i] == "ALA"){
		colorL = append(colorL,colorApolar)
	}
	if(AA[i] == "GLY"){
		colorL = append(colorL,colorOther)
	}
	if(AA[i] == "ILE"){
		colorL = append(colorL,colorApolar)
	}
	if(AA[i] == "LEU"){
		colorL = append(colorL,colorApolar)
	}
	if(AA[i] == "MET"){
		colorL = append(colorL,colorApolar)
	}
	if(AA[i] == "PRO"){
		colorL = append(colorL,colorOther)
	}
	if(AA[i] == "VAL"){
		colorL = append(colorL,colorApolar)
	}
	if(AA[i] == "HOH"){
		colorL = append(colorL,colorWater)
	}
}

data = cbind(data, colorV)

nbDistance = dim(data)[2] - 4

distance = c()

for(i in 0:nbDistance){
	deb = 2 + 0.5*i
	distance = c(distance,paste(deb,"<D<",deb+0.5, sep = ""))

}

print (distance)
#distance = c("2Å<D<2.5Å", "2.5Å<D<3Å", "3Å<D<3.5Å", "3.5Å<D<4Å", "4Å<D<4.5Å")

#data = data[order(data[,5],decreasing = T),]

nameRes = data[,1]
data = data[,-1]
nbCol = dim(data)[2] 
colorFig = as.vector(data[,nbCol])
data = data[,-nbCol]
#length (colorFig)

data = data[,-1]

title = paste("Amino acid types having at least one side chain atom within selected distances (D) of Nref from ",structure," amine" ,sep = "")

large = dim(data)[1]*50
png(filename=paste(file,".png",sep = ""),width=as.integer(large))
barplot(as.matrix(data), main = title, ylab= "Number of amino acids", beside = TRUE, col=colorFig, xlab = "Distance from Nref (Å)", names.arg = distance)
legend( 2 , max(data) , legend=AA, bty="n", fill=colorL)

dev.off()
