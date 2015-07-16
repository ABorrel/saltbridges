'''
Created on 12 mai 2015

@author: borrel
'''

# personnal
import loadFile
import parsing
import searchPDB
import writeFile
import writePDBfile

# 
from os import listdir
from copy import deepcopy
 
 
 
 
def SearchSaltBridges (pr_data, pr_result, dist_thresold = 5.0, debug = 1): 
    
    
    # do same analysis that PDB 
    l_PDB_file = listdir(pr_data)
    
    pr_summary = pr_result + "Sum/"
    p_file_lig = pr_result + "ligInPDB"
    file_lig = open (p_file_lig, "w")
    

    # Write file ligand  
    
    l_ligand = listdir(pr_data)
    for ligand in l_ligand : 
        print ligand
    
        l_file = listdir(pr_data + ligand + "/")
    
    
        for PDB_file in l_file :
            print PDB_file
        
            l_lig_PDB = parsing.retrieveListLigand(pr_data + ligand + "/" + PDB_file)
            print l_lig_PDB
            if len (l_lig_PDB) == 0 : 
                continue
            
            # write temp control
            file_lig.write (l_lig_PDB[0] + "\t" + pr_data + str (ligand) + "/" + PDB_file + "\n")
    file_lig.close ()
     
     
        
    l_lig = loadFile.resultFilterPDBLigand(p_file_lig)
    
    print "==>", l_lig
    
    nb_lig = len(l_lig)
    
    # ##Count Structure
    d_neighbor = {}
    l_neighbor_global = []
        
        
    # ##Write summary file
    d_files_summary = writeFile.openFileSummary(pr_summary)# sumary result
    
    # Creation of summary file 
    i = 0
    while i < nb_lig :
        if debug: print "Ligand: " + str(l_lig[i]["name"]) + " " + str(i) + " " + str(nb_lig)
        nb_PDB = len(l_lig[i]["PDB"])
        
        # take only one PDB by ligand
            
        j = 0
        while j < nb_PDB : 
            name_PDB = l_lig[i]["PDB"][j]
            print name_PDB
            print l_lig[i]["name"]
            
            l_atom_ligand = loadFile.ligandInPDB(l_lig[i]["PDB"][j], l_lig[i]["name"])
            
            # search neighbor for every atom in ligand selected
            searchPDB.globalNeighbors(dist_thresold, l_atom_ligand, name_PDB, l_neighbor_global)
            # search neighbor for interest 
            searchPDB.interestGroup(dist_thresold, l_atom_ligand, name_PDB, d_neighbor)
            
            j = j + 1
        i = i + 1
            
    
    
    writeFile.neighborStruct(d_neighbor, l_neighbor_global, d_files_summary)
    writeFile.closeFileSummary(d_files_summary)
    
    # case where load directly substructure => why do not load directly in dictionnary
    d_neighbor["global"] = l_neighbor_global
    return d_neighbor
    
 
 
 
def ControlPDBFormat (p_PDB):
    
    
    l_res = ["ILE", "LEU", "LYS", "PHE", "TYR", "VAL", "SER", "MET", "ARG", "TRP", "PRO", "GLY", "GLU", "ASN", "HIS", "ALA", "ASP", "GLN", "THR", "CYS"]
    l_atom_PDB = parsing.loadCoordSectionPDB(p_PDB)
    
    filout = open (p_PDB, "w")
    
    l_atom_het = []
    for atom_PDB in l_atom_PDB : 
        if atom_PDB["resName"] in l_res : 
            writePDBfile.coordinateSection(filout, [atom_PDB], "ATOM", header = 0)
        else :
            l_atom_het.append (deepcopy(atom_PDB)) 
    
    writePDBfile.coordinateSection(filout, l_atom_het, recorder = "HETATM", header = 0, connect_matrix = 1)
    
    filout.close ()
            
    

def SearchChemicalSubstruct (pr_data, pr_result, control = 0):
    
    # substructure search 
    
    p_filout = pr_result + "findSruct"
    
    filout = open (p_filout, "w")
    l_ligand = listdir(pr_data)
    
    print l_ligand
    
    for ligand in l_ligand : 
        print ligand
    
        l_file = listdir(pr_data + ligand + "/")
        
        for f in l_file : 
            if ligand == "ZM241385" : 
                group = f
            else : 
                group = f.split ("_")[-2]
            p_file_PDB = pr_data + ligand + "/" + f
            print p_file_PDB
            
            if control == 1 : 
                ControlPDBFormat(p_file_PDB)
            
            l_atom = parsing.loadCoordSectionPDB(p_file_PDB, remove_H = 1)
            if ligand == "ZM241385" : 
                ll_atom_lig = parsing.retrieveLigand(l_atom, "ZMA")
            else : 
                ll_atom_lig = parsing.retrieveLigand(l_atom, "RES")
    
            for l_atom_lig in ll_atom_lig : 
                
                l_subs = searchPDB.interestStructure(l_atom_lig, more_flex = 0)
                
                filout.write (str (ligand) + "\t" + str (f) + "\t" + str (group) + "\t" + " ".join(l_subs) + "\n")
    
    filout.close ()
       

# best group ergomine -> 2128 (2013) && 6882 (2013)
# best group taladegid -> 7527 (2013) 
# best group eticlopride -> 1285 (2010)
# ZM241385 => mod7msp
    