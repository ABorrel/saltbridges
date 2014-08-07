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


def construction(name_folder, RX = 3.00, RFree = 0.25):
    """
    Dataset construction
    in : - open file result of filter ligand PDB
    out : - log file
          - dataset file -> ligand with associated PDB
    """

    
    pr_result = repertory.result(name_folder)
    # check dataSet exist !!!!!!
    # short cut
    list_file_dataset = repertory.retriveDataSetFile (pr_result)
    if len(list_file_dataset) != 0 : 
        return list_file_dataset


    start, logFile = log.initAction("Dataset construction")
    
    ligandInPDB = loadFile.resultLigandPDB(pr_result + "resultLigandInPDB")
 
    nbLigand = len(ligandInPDB.keys())
    listligand = ligandInPDB.keys()
     
    i = 0
    while (i < nbLigand):
        nameLigand = ligandInPDB.keys()[i]
        PDB_ref = ligandInPDB[nameLigand][0]
 
        # load ligand
        l_atom_lig_ref = loadFile.ligandInPDBConnectMatrixLigand(PDB_ref, nameLigand)
        l_interest_sub = searchPDB.interestStructure(l_atom_lig_ref) # search interest structure
        if l_interest_sub == []:
            del ligandInPDB[nameLigand]
            nbLigand = nbLigand - 1
            continue
        else : 
            # control dataset quality
            checkPDBfile.checkPDB(ligandInPDB[nameLigand], nameLigand, RX, RFree)
            if ligandInPDB[nameLigand] == []:
                del ligandInPDB[nameLigand]
                nbLigand = nbLigand - 1
                continue
            else :
                i = i + 1
                
                
    # structure and file dataset and control RX + length bond
    l_file_datast = buildStructDataset (ligandInPDB, pr_result)
    log.endAction("Dataset construction", start, logFile)      
    return l_file_datast 
                
                
def buildStructDataset (d_lig_PDB, pr_init) :                
    
    pr_quality = repertory.parsingDataset(pr_init)     
    pr_bond = repertory.lengthBond (pr_init)
    ld_cn = []
    ld_coplar = []
    l_RX = []
    l_RFree = []           
    struct_dataset = structure.resolutionFilter()
            
                
    for lig in d_lig_PDB.keys () :
        for PDB in d_lig_PDB[lig]:
            l_atom_lig = loadFile.ligandInPDB(PDB, lig)
            controlLenCNBond(l_atom_lig, ld_cn)  # only one PDB by ligands
            controlCoplanarTertiaryAmine(l_atom_lig, ld_coplar)  # only one PDB by ligands
            
            appendStruct(PDB, lig, struct_dataset)
                
            l_quality = parsing.resolution(PDB)
            l_RX.append (l_quality[0])
            l_RFree.append (l_quality[1])
    
    
    
    # bond length control
    p_CN = writeFile.listFloat(ld_cn, pr_bond + "lengthCNallLigand")
    p_Cop = writeFile.listFloat(ld_coplar, pr_bond + "distanceCoplanarNter")
    
    # quality
    p_rx = writeFile.listFloat(l_RX, pr_quality + "RX")
    p_rf = writeFile.listFloat(l_RFree, pr_quality + "RFree")
    
    runScriptR.histDistance(p_CN, "CN")
    runScriptR.histDistance(p_Cop, "coplar")
    
    runScriptR.histDistance(p_rx, "RX")
    runScriptR.histDistance(p_rf, "RFree")
    
    return writeFile.resultFilterLigandPDB(struct_dataset, pr_init )
    
    


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

    type = parsing.header(filePDB)
    RX = parsing.resolution(filePDB)[0]

    if type == "SOLUTION":
        try:
            structDataset["NMR"][nameLigand].append(filePDB)
        except:
            structDataset["NMR"][nameLigand] = []
            structDataset["NMR"][nameLigand].append(filePDB)

    else:
        if RX > 3.00:
            try:
                structDataset["OUT"][nameLigand].append(filePDB)
            except:
                structDataset["OUT"][nameLigand] = []
                structDataset["OUT"][nameLigand].append(filePDB)

        if RX <= 3.00:
            try:
                structDataset["3.00"][nameLigand].append(filePDB)
            except:
                structDataset["3.00"][nameLigand] = []
                structDataset["3.00"][nameLigand].append(filePDB)

        if RX <= 2.50:
            try:
                structDataset["2.50"][nameLigand].append(filePDB)
            except:
                structDataset["2.50"][nameLigand] = []
                structDataset["2.50"][nameLigand].append(filePDB)


        if RX <= 2.00:
            try:
                structDataset["2.00"][nameLigand].append(filePDB)
            except:
                structDataset["2.00"][nameLigand] = []
                structDataset["2.00"][nameLigand].append(filePDB)
        
        if RX <= 1.50:
            try:
                structDataset["1.50"][nameLigand].append(filePDB)
            except:
                structDataset["1.50"][nameLigand] = []
                structDataset["1.50"][nameLigand].append(filePDB)
                
