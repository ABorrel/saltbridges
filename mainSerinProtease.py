from os import listdir
import managePDB
import repertory
import loadFile
import searchPDB
import tool
from string import lower
import retrieveAtom
import writePDBfile

repData = repertory.datasetSerineProtease()
listRep = listdir(repData)

#for rep in listRep :
#    rep = repData + rep + "/"
#    managePDB.appendExtention(rep)

###to do fonction that retrieve ligands in pdb files#######


ligandWithPDB = loadFile.resultLigandPDB(repertory.resultSerineProtease() + "resultLigandInPDB")
listLigand = ligandWithPDB.keys()
print listLigand

nbLigand = len(listLigand)
print nbLigand

j = 0
while j < nbLigand :
    ligand = listLigand[j]
    print ligand , j
    pdb = ligandWithPDB[ligand][0].split("_")[2]
    pdb = lower(pdb)
    atomLigand = loadFile.ligandInPDBConnectMatrixLigand(pdb, ligand)
    writePDBfile.globalStruct(ligand + "_" + pdb + ".pdb", atomLigand)
    nbAtom = len(atomLigand)
    #print atomLigand

    searchPDB.cycleGlobal(atomLigand)
#    print atomLigand
    listResult = []
    for atom in atomLigand : 
        if atom["element"] == "N" or atom["element"] == "O" or atom["element"] == "C":# impact point
            if atom["cycle"] == 1 : 
                cycleRetrieve = retrieveAtom.cycle(atom, atomLigand)
                listResult.append(cycleRetrieve)
                continue
            for serialAtom in atom["connect"] : 
                atomRetrieve = retrieveAtom.serial(serialAtom, atomLigand)
                if atomRetrieve == 0 :
                    continue
                if atomRetrieve["cycle"] == 1 : 
                    cycleRetrieve = retrieveAtom.cycle(atom, atomLigand)
                    listResult.append(cycleRetrieve)
                    continue
                
    tool.checkListResult(listResult)
    nbStruct = 1
    for element in listResult : 
        writePDBfile.globalStruct(ligand + "_" + pdb + "_" + str(nbStruct) + ".pdb", element)
        nbStruct = nbStruct + 1
    
    j = j + 1

