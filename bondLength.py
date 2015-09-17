'''
Created on 13 mai 2015

@author: borrel
'''

import pathManage
import loadFile
import parsing
import runScriptR
import retrieveAtom
import calcul
import searchPDB

from os import path


def GlobalBondLength (name_database, RX_thresold = 1.5):
    
    # directory
    pr_result = pathManage.result(name_database + "/CXbound" + str (RX_thresold))
    pr_database = pathManage.result(name_database)
    
    # filout with distance
    p_CN = pr_result + "distanceCN"
    p_CO = pr_result + "distanceCO"
    p_CC = pr_result + "distanceCO"
    filout_CN = open (pr_result + "distanceCN", "w")
    filout_CO = open (pr_result + "distanceCO", "w")
    filout_CC = open (pr_result + "distanceCC", "w")
    
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
                
                l_distCN = BondLengthCandX (l_atom_lig, "N")
                l_distCO = BondLengthCandX (l_atom_lig, "O")
                l_distCC = BondLengthCandX (l_atom_lig, "C")
                if l_distCN != [] : 
                    filout_CN.write ("\n".join (l_distCN) + "\n")
                
                if l_distCO != [] : 
                    filout_CO.write ("\n".join (l_distCO) + "\n")    
                
                if l_distCC != [] : 
                    filout_CC.write ("\n".join (l_distCC) + "\n")                
                    
                
                # take only one PDB by ligand not more
                i = i + 1
                continue
            i = i + 1
        


    filout_CO.close ()
    filout_CN.close ()
    filout_CC.close ()
    
    runScriptR.histDistance(p_CN, "CN")
    runScriptR.histDistance(p_CO, "CO") 
    runScriptR.histDistance(p_CC, "CC") 



def BondLengthCandX (l_atom_lig, X_element) : 
    
    l_dist = []
    l_serialX = searchPDB.ListSerialElement(l_atom_lig, str(X_element))
    
    for serialX in l_serialX : 
#         print serialX
        l_atom_connectX, connect = retrieveAtom.atomConnect(l_atom_lig, serialX)
        
        i = 0
        nb_connect = len (connect)
        while i < nb_connect : 
            if len(l_atom_connectX) >= 2 : 
                dist = calcul.distanceTwoatoms(l_atom_connectX[0], l_atom_connectX[1])
#                 print dist
                if dist != 100 : 
                    l_dist.append (str(dist))
            i = i + 1
    
    return l_dist
                
        
 
