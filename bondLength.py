'''
Created on 13 mai 2015

@author: borrel
'''

import pathManage
import loadFile
import parsing
import searchPDB
import runScriptR



from os import path


def BondCNAndCO (name_database, RX_thresold = 1.5):
    
    # directory
    pr_result = pathManage.result(name_database + "/CNbound" + str (RX_thresold))
    pr_database = pathManage.result(name_database)
    
    # filout with distance
    p_CN = pr_result + "distanceCN"
    p_CO = pr_result + "distanceCO"
    filout_CN = open (pr_result + "distanceCN", "w")
    filout_CO = open (pr_result + "distanceCO", "w")
    # load PDB with logand
    if not path.exists(pr_database + "resultLigandInPDB") : 
        print "ERROR => file with ligand and PDB does not exist"
        return
    else : 
        d_lig_PDB = loadFile.LigandInPDB(pr_database + "resultLigandInPDB")
    
 
    nb_lig = len(d_lig_PDB.keys())
    print d_lig_PDB.keys()
    
    i = 0
    while (i < nb_lig):
        name_lig = d_lig_PDB.keys()[i]
        
        l_PDB = d_lig_PDB[name_lig]
        
        for PDB in l_PDB : 
            # controle RX
            RX = parsing.Quality(PDB)[0]
#             print RX
            
            if RX <= RX_thresold : 
                l_atom_lig = loadFile.ligandInPDBConnectMatrixLigand(PDB, name_lig)
                
                l_distanceCN = searchPDB.BondLengthCN (l_atom_lig)
                l_distanceCO = searchPDB.BondLengthCO (l_atom_lig)
                if l_distanceCN != [] : 
                    filout_CN.write ("\n".join (l_distanceCN) + "\n")
                
                if l_distanceCO != [] : 
                    filout_CO.write ("\n".join (l_distanceCO) + "\n")    
                    
                
                # take only one PDB by ligand not more
                i = i + 1
                continue
            i = i + 1
        


    filout_CO.close ()
    filout_CN.close ()
    runScriptR.histDistance(p_CN, "CN")
    
    runScriptR.histDistance(p_CO, "CO") 



