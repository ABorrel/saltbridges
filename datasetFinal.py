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


    # load structure    
    d_lig_PDB = loadFile.LigandInPDB(pr_database + "resultLigandInPDB")
    
    # list of ligand
    l_lig = d_lig_PDB.keys()
 
    nb_lig = len(d_lig_PDB.keys())
    print d_lig_PDB.keys()
    
    i = 0
    while (i < nb_lig):
        name_lig = l_lig[i]
        PDB_ref = d_lig_PDB[name_lig][0]
        if debug : print PDB_ref, name_lig, i
 
        # load ligand
        try : l_atom_lig_ref = loadFile.ligandInPDBConnectMatrixLigand(PDB_ref, name_lig)
        except : 
            i = i + 1
            continue
        
        # search substructure interest
        l_interest_sub = searchPDB.interestStructure(l_atom_lig_ref) # search interest structure
        if debug : print "Interest substructure", l_interest_sub
        if l_interest_sub == []:
            del d_lig_PDB[name_lig]
        else : 
            # control dataset quality
            l_PDB = checkPDBfile.CheckComplexQuality(d_lig_PDB[name_lig], name_lig, RX, RFree, one_PDB_by_lig)
            d_lig_PDB[name_lig] = l_PDB
            if d_lig_PDB[name_lig] == []:
                del d_lig_PDB[name_lig]
        i = i + 1
        
        
    if debug == 1 : print "struct ligand =>", d_lig_PDB
                
    # structure and file dataset and control RX + length bond
    WriteDataset (d_lig_PDB, pr_result)
    #log.endAction("Dataset Builder", start, logFile)   
    
       
    return  Builder(name_database, RX , RFree , one_PDB_by_lig , debug = 1) 





                
def WriteDataset (d_lig_PDB, pr_init, debug = 1) :                
    """
    Write folder with dataset 
    -> run quality criteria analysis
    
    """
    
    pr_quality = pathManage.parsingDataset(pr_init)     
    l_RX = []
    l_RFree = []           
    d_dataset = structure.resolutionFilter()
            
                
    for lig in d_lig_PDB.keys () :
        if debug == 1 : 
            print "===Write dataset-1==="
            print lig, "ligand"
            print d_lig_PDB[lig], "struct data"
        for PDB in d_lig_PDB[lig]:
            if debug == 1 :
                print "===Write dataset-2===" 
                print PDB, lig, "PDB + lig"
                print d_dataset, "append struct"
            appendStruct(PDB, lig, d_dataset)
                
            l_quality = parsing.Quality(PDB)
            l_RX.append (l_quality[0])
            l_RFree.append (l_quality[1])
    
    
    # quality
    p_rx = writeFile.listFloat(l_RX, pr_quality + "RX")
    p_rf = writeFile.listFloat(l_RFree, pr_quality + "RFree")
    
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
        
        
def appendStruct(PDB_ID, name_lig, d_dataset, debug = 0):
    """Append name PDB_ID in structure dataset
    in: PDB_ID, name ligand, structure save dataset
    out: append new ligand in structure dataset"""

    method = parsing.methodStructure(PDB_ID)
    RX = parsing.Quality(PDB_ID)[0]
    
    if debug : 
        print "**RX**", RX
        print d_dataset, "input -> append structure"
        print "**method**", method

    if method == "XRAY":
        if float(RX) > 3.00:
            if debug: print "in OUT"
            if not name_lig in d_dataset["OUT"].keys () : 
                d_dataset["OUT"][name_lig] = []
            d_dataset["OUT"][name_lig].append(PDB_ID)
            
        if float(RX) <= 1.50:
            if debug: print "in 1.5"
            if not name_lig in d_dataset["1.50"].keys () : 
                d_dataset["1.50"][name_lig] = []
            d_dataset["1.50"][name_lig].append(PDB_ID)  
            
        if float(RX) <= 2.00:
            if debug: print "in 2.0"
            if not name_lig in d_dataset["2.00"].keys () : 
                d_dataset["2.00"][name_lig] = []
            d_dataset["2.00"][name_lig].append(PDB_ID)   

        if float(RX) <= 2.50:
            if debug: print "in 2.5"
            if not name_lig in d_dataset["2.50"].keys () : 
                d_dataset["2.50"][name_lig] = []
            d_dataset["2.50"][name_lig].append(PDB_ID)
            
        if float(RX) <= 3.00:
            if debug: print "in 3.0"
            if not name_lig in d_dataset["3.00"].keys () : 
                d_dataset["3.00"][name_lig] = []
            d_dataset["3.00"][name_lig].append(PDB_ID)   

