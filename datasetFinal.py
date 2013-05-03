import checkPDBfile
import loadFile
import log
import searchPDB
import structure
import writeFile
import repertory
import retrieveAtom
import calcul
import toolSubstructure
import runScriptR
import parsing


def construction(name_folder_dataset):
    """
    Dataset construction
    in : - open file result of filter ligand PDB
    out : - log file
          - dataset file -> ligand with associated PDB
    """

    
    rep_dataset = repertory.result(name_folder_dataset)
    # check dataSet exist !!!!!!
    list_file_dataset = repertory.retriveDataSetFile (rep_dataset)
    if len(list_file_dataset) != 0 : 
        return list_file_dataset


    start, logFile = log.initAction("Dataset construction")
    
    ligandInPDB = loadFile.resultLigandPDB(rep_dataset + "resultLigandInPDB")
    resultFilterPDB = structure.resolutionFilter()
 
    nbLigand = len(ligandInPDB.keys())
    listligand = ligandInPDB.keys()
    listDistanceCN = []
    listDistanceCoplanar = []
     
    i = 0
    while (i < nbLigand):
        nameLigand = listligand[i]
        PDBFile = ligandInPDB[nameLigand][0]
        print nameLigand, PDBFile, i
 
        listAtomLigand = loadFile.ligandInPDBConnectMatrixLigand(PDBFile, nameLigand)
        controlLenCNBond(listAtomLigand, listDistanceCN)  # only one PDB by ligands
        controlCoplanarTertiaryAmine(listAtomLigand, listDistanceCoplanar)  # only one PDB by ligands
        listStruct = searchPDB.interestStructure(listAtomLigand) # search interest structure
         
        if listStruct == []:
            i = i + 1
            continue
        else:
            checkPDBfile.checkPDB(ligandInPDB[nameLigand], nameLigand)  # seq check + ligand hooked
 
        if ligandInPDB[nameLigand] == []:
            i = i + 1
            continue
 
        for file_pdb in ligandInPDB[nameLigand]:  # append PDB file_pdb
            appendStruct(file_pdb, nameLigand, resultFilterPDB)
 
        i = i + 1

    list_files_dataset = writeFile.resultFilterLigandPDB(resultFilterPDB, rep_dataset )
    
    dir_distance = repertory.resultDistance(rep_dataset)
    writeFile.resultLengthCNBond(listDistanceCN, "lengthCNallLigand", dir_distance )
    runScriptR.histDistance("lengthCNallLigand", "CN", "PDB", dir_distance )
    writeFile.resultCoplanar(listDistanceCoplanar, "distanceCoplanar", dir_distance )
    runScriptR.histDistance("distanceCoplanar", "coplar", "PDB", dir_distance )
    
    log.endAction("Dataset construction", start, logFile)
    
    return  list_files_dataset


def controlLenCNBond (listAtomLigand, listDistance):
    """For each CN bond in list atom ligand, retrieve distance between C and N
    in: list atom ligand, list distance retrieve
    out: append in list ligand the distance"""


    for atom in listAtomLigand:
        if atom["element"] == "N":
            matrixConnect = atom["connect"]

            for serial in matrixConnect:
                atomConect = retrieveAtom.serial(serial, listAtomLigand)

                if atomConect != 0  and atomConect["element"] == "C":
                        distance = calcul.distanceTwoatoms(atom, atomConect)
                        listDistance.append(distance)



def controlCoplanarTertiaryAmine(listAtomLigand, listDistanceCoplanar) : 
    """For each matrix connect N, C, C, C the coplanar distance
    in: list atom in ligand, list distance coplar retrieve
    out: append distance in list distance coplanar"""
    
    listSerialNitrogen = searchPDB.listAtomType(listAtomLigand, "N")
    for serialNitrogen in listSerialNitrogen:
        listAtomConnectNitrogen, conect = retrieveAtom.atomConnect(listAtomLigand, serialNitrogen)
        
        connectMatrixElement = toolSubstructure.matrixElement(listAtomConnectNitrogen)
        if connectMatrixElement == ["N", "C", "C", "C"]:
            if toolSubstructure.checkConectOnlyC(listAtomConnectNitrogen[1], listAtomLigand) == 1 and toolSubstructure.checkConectOnlyC(listAtomConnectNitrogen[2], listAtomLigand) == 1 and toolSubstructure.checkConectOnlyC(listAtomConnectNitrogen[3], listAtomLigand) == 1:
                try :
                    distance = calcul.coplanar(listAtomConnectNitrogen[0], listAtomLigand)
                    if distance != None :
                        listDistanceCoplanar.append(distance)
                except :
                    pass
        
        
def appendStruct(filePDB, nameLigand, structDataset):
    """Append name filePDB in structure dataset
    in: filePDB, name ligand, structure save dataset
    out: append new ligand in structure dataset"""

    type = parsing.methods(filePDB)
    resolution = parsing.resolution(filePDB)

    if type == "SOLUTION":
        try:
            structDataset["NMR"][nameLigand].append(filePDB)
        except:
            structDataset["NMR"][nameLigand] = []
            structDataset["NMR"][nameLigand].append(filePDB)

    else:
        if resolution > 3.00:
            try:
                structDataset["OUT"][nameLigand].append(filePDB)
            except:
                structDataset["OUT"][nameLigand] = []
                structDataset["OUT"][nameLigand].append(filePDB)

        if resolution <= 3.00:
            try:
                structDataset["3.00"][nameLigand].append(filePDB)
            except:
                structDataset["3.00"][nameLigand] = []
                structDataset["3.00"][nameLigand].append(filePDB)

        if resolution <= 2.50:
            try:
                structDataset["2.50"][nameLigand].append(filePDB)
            except:
                structDataset["2.50"][nameLigand] = []
                structDataset["2.50"][nameLigand].append(filePDB)


        if resolution <= 2.00:
            try:
                structDataset["2.00"][nameLigand].append(filePDB)
            except:
                structDataset["2.00"][nameLigand] = []
                structDataset["2.00"][nameLigand].append(filePDB)
        
        if resolution <= 1.50:
            try:
                structDataset["1.50"][nameLigand].append(filePDB)
            except:
                structDataset["1.50"][nameLigand] = []
                structDataset["1.50"][nameLigand].append(filePDB)
                
