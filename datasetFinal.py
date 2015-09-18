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
    
    nb_lig = len(d_lig_PDB.keys())
    print "NB ligand included database:", nb_lig
    
    # print d_lig_PDB.keys().index("HSO") -> search index ligand
    
    i = 0
    while i < nb_lig:
        name_lig = d_lig_PDB.keys()[i]
        
        #######################################
        # step 1 search chemical substructure #
        #######################################
        PDB_ref = d_lig_PDB[name_lig][0]
        if debug : print PDB_ref, name_lig, i, nb_lig
        # if not possible to load the ligand -> remove lig
        try : l_atom_lig_ref = loadFile.ligandInPDBConnectMatrixLigand(PDB_ref, name_lig)
        except : 
            if debug == 1 : print "Exit => load ligand-l59"
            del d_lig_PDB[name_lig]
            nb_lig = nb_lig - 1
            continue
        
        # search substructure interest
        l_interest_sub = searchPDB.interestStructure(l_atom_lig_ref) # search interest structure
        if debug : print "Interest substructures in " + str(name_lig) + "-" + str (PDB_ref) + " " + "-".join(l_interest_sub)
        if l_interest_sub == []:
            if debug == 1 : print "Exit => Not substructure-l68"
            del d_lig_PDB[name_lig]
            nb_lig = nb_lig - 1
            continue
        
        #######################################################
        # Step 2 Control quality of PDB + ligand hooked + option one #
        #######################################################
        else : 
            # control dataset quality
            if debug : print "List PDB checked -> ", d_lig_PDB[name_lig]
            l_PDB = checkPDBfile.CheckComplexQuality(d_lig_PDB[name_lig], name_lig, RX, RFree, one_PDB_by_lig)
            # remove the entrance key with the ligand
            if l_PDB == []:
                if debug == 1 : print "Exit => Not No PDB selected-l82"
                del d_lig_PDB[name_lig]
                nb_lig = nb_lig - 1
                continue
            else : 
                d_lig_PDB[name_lig] = l_PDB
        i = i + 1
        
        
    if debug == 1 : print "Number of ligand selected =>", nb_lig
                
    # structure and file dataset and control RX + length bond
    WriteDataset (d_lig_PDB, pr_result)
    
       
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
            print d_lig_PDB[lig], "list of PDB"
        for PDB in d_lig_PDB[lig]:
            if debug == 1 :
                print "===Write dataset-2===" 
                print PDB, lig, "PDB + lig"
                print d_dataset, "append struct"
            quality = parsing.Quality(PDB)
            RX = quality[0]
            RFree = quality[1]
            BuilderDatasetDict(PDB, lig, RX, d_dataset)
            l_RX.append (RX)
            l_RFree.append (RFree)
    
    
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



def BuilderDatasetDict(PDB_ID, name_lig, RX, d_dataset, debug = 1):
    """Append name PDB_ID in structure dataset
    in: PDB_ID, name ligand, structure save dataset
    out: append new ligand in structure dataset"""

    if debug : 
        print "**RX**", RX
        print d_dataset, "input -> append structure"
        
    print d_dataset.keys ()
    for RX_build in d_dataset.keys () : 
        if float(RX) <= float(RX_build) : 
            if not name_lig in d_dataset[RX_build].keys () : 
                d_dataset[RX_build][name_lig] = []
            d_dataset[RX_build][name_lig].append(PDB_ID)
        else : 
            continue
            
