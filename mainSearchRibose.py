import searchPDB
import loadFile
import repertory
import tool

searchPDB.ligands()
listLigandInPDB = loadFile.resultLigandPDB(repertory.result() + "resultLigandInPDB")

listLigand = listLigandInPDB.keys()
nbLigand = len(listLigand)
listLigandWithRibose = []

i = 0
while i < nbLigand :
    print listLigand[i], i
    ligandID = listLigand[i]
    listAtomLigand = loadFile.ligandInPDB(listLigandInPDB[ligandID][0], ligandID) 
    if searchPDB.riboseFromADPRibose(listAtomLigand) == 1 : 
        listLigandWithRibose.append(listLigand[i])
    
    i = i +1 
        
fileResultRibose = open(repertory.result()+"riboseList","w")

print len(listLigandWithRibose)

for ligand in listLigandWithRibose : 
    fileResultRibose.write(ligand + "\n")

fileResultRibose.close()   





