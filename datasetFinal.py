import checkPDBfile
import loadFile
import log
import searchPDB
import structure
import writeFile
import pathManage
import retrieveAtom
import calcul
import toolSubstructure
import runScriptR
import parsing


def Builder(name_database, RX = 3.00, RFree = 0.25, one_PDB_by_lig = 0, debug = 1):
    """
    Dataset Builder
    in : - open file result of filter ligand PDB
    out : - log file
          - dataset file -> ligand with associated PDB
    """
    
    if one_PDB_by_lig == 0 : 
        name_dataset = name_database + "/" + str (RX) + "_" + str (RFree) + "_multiPDB"
    else : 
        name_dataset = name_database + "/" + str (RX) + "_" + str (RFree) + "_uniquePDB"
    
    pr_database = pathManage.result(name_database)
    pr_result = pathManage.result(name_dataset)
    if debug : print "== Path result " + pr_result + "==\n"
    
    # check dataSet exist !!!!!!
    # short cut
    l_file_dataset = pathManage.retriveDataSetFile (pr_result)
    if len(l_file_dataset) != 0 : 
        return l_file_dataset


    #start, logFile = log.initAction("Dataset Builder")
    
    d_lig_PDB = loadFile.LigandInPDB(pr_database + "resultLigandInPDB")
 
    nb_lig = len(d_lig_PDB.keys())
     
    i = 0
    while (i < nb_lig):
        name_lig = d_lig_PDB.keys()[i]
        PDB_ref = d_lig_PDB[name_lig][0]
 
        # load ligand
        l_atom_lig_ref = loadFile.ligandInPDBConnectMatrixLigand(PDB_ref, name_lig)
        l_interest_sub = searchPDB.interestStructure(l_atom_lig_ref) # search interest structure
        if l_interest_sub == []:
            del d_lig_PDB[name_lig]
            nb_lig = nb_lig - 1
            continue
        else : 
            # control dataset quality
            d_lig_PDB[name_lig] = checkPDBfile.CheckComplexQuality(d_lig_PDB[name_lig], name_lig, RX, RFree, one_PDB_by_lig)
            if d_lig_PDB[name_lig] == []:
                del d_lig_PDB[name_lig]
                nb_lig = nb_lig - 1
                continue
            else :
                i = i + 1
                
                
    # structure and file dataset and control RX + length bond
    WriteDataset (d_lig_PDB, pr_result, one_PDB_by_lig)
    #log.endAction("Dataset Builder", start, logFile)   
    
       
    return  Builder(name_database, RX , RFree , one_PDB_by_lig , debug = 1) 
                
                
def WriteDataset (d_lig_PDB, pr_init) :                
    """
    Write folder with dataset 
    -> run quality criteria analysis
    
    """
    
    pr_quality = pathManage.parsingDataset(pr_init)     
    pr_bond = pathManage.lengthBond (pr_init)
    ld_cn = []
    ld_coplar = []
    l_RX = []
    l_RFree = []           
    d_dataset = structure.resolutionFilter()
            
                
    for lig in d_lig_PDB.keys () :
        for PDB in d_lig_PDB[lig]:
            l_atom_lig = loadFile.ligandInPDB(PDB, lig)
            controlLenCNBond(l_atom_lig, ld_cn)  # only one PDB by ligands
            controlCoplanarTertiaryAmine(l_atom_lig, ld_coplar)  # only one PDB by ligands
            
            appendStruct(PDB, lig, d_dataset)
                
            l_quality = parsing.Quality(PDB)
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
    
    return writeFile.resultFilterLigandPDB(d_dataset, pr_init )
    
    


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

    method = parsing.header(filePDB)
    RX = parsing.Quality(filePDB)[0]

    if method == "SOLUTION":
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
                
