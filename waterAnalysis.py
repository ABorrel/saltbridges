"""
BORREL Alexandre
04-2013
"""
from os import listdir, path
from re import search
import parsing
import parseNACCESS
import repertory
import structure
import loadFile


def resolutionWater (list_PDB, path_folder_result, limit_acc = 20.0):
    
    path_folder_database = repertory.pathDitrectoryPDB()
    path_filout = path_folder_result + "statwater_" + str (limit_acc) + ".dat"
    if path.isfile(path_filout) and path.getsize(path_filout) > 0 : 
        pass
    else : 
        filout = open (path_filout, "w")
        
        for PDB_ID in list_PDB : 
            path_file_PDB = path_folder_database + PDB_ID + ".pdb"
            path_file_asa = path_folder_database + PDB_ID + ".asa"
            path_file_rsa = path_folder_database + PDB_ID + ".rsa"
                    
            if not path.isfile(path_file_asa) or not path.isfile(path_file_rsa) : 
                print "Error NACCESS", path_file_asa
                continue
                    
            resolution = parsing.resolution(PDB_ID)
            if resolution == 1000.0 : 
                continue
            number_residue_exposed, number_residue = parseNACCESS.numberResExposed(path_file_rsa, limit_acc)
            number_water = parsing.countH2O (path_file_PDB) 
                    
            filout.write (str (PDB_ID) + "\t" + str (resolution) + "\t" + str (number_residue_exposed) + "\t" + str (number_residue) + "\t" + str (number_water) + "\n")
        filout.close ()
    return path_filout




def resolutionByStructure (name_dataset) :
    
    l_structure = structure.listStructure()
    l_path = []
    
    for strut in l_structure :
        l_path.append (repertory.result(name_dataset) + "water_" + strut + ".dat")
        filout = open (repertory.result(name_dataset) + "water_" + strut + ".dat", "w") 
        l_file_summary = repertory.retrieveSummaryFile (strut, name_dataset)
        
        list_global = []
        for path_summary in l_file_summary : 
            print path_summary
            list_interest_atom = loadFile.loadSummary (path_summary)
            for interest_atom in list_interest_atom : 
                if not interest_atom in list_global : 
                    list_global.append (interest_atom)
                    
        for atom_interest in list_global : 
            rx, i_nb_atom, s_PDB = searchCountH2O(atom_interest)
            filout.write ("%s\t%s\t%s\n"%(s_PDB, rx, i_nb_atom))
        filout.close ()
    return l_path


def searchCountH2O (atom_interest):
    
    PDB = atom_interest["PDB"]
    rX = parsing.resolution(PDB)
    i_nb_H2O = nbH2O(atom_interest["neighbors"])
    
    return rX, i_nb_H2O, PDB    
    
    
        
def nbH2O (neigbor):
    
    out = 0
    for d_neighbor in neigbor : 
        if d_neighbor["resName"] == "HOH" : 
            out = out + 1
    return out
    
    
        
        
        
  
        
    
    
    
