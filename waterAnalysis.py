"""
BORREL Alexandre
04-2013
"""
from os import listdir, path
from re import search
import parsing
import parseNACCESS
import pathManage
import structure
import loadFile


def resolutionWater (l_PDB, pr_result, limit_acc = 00.0):
    
    pr_PDB = pathManage.pathDitrectoryPDB()
    p_filout = pr_result + "statwater_" + str (limit_acc) + ".dat"
    #if path.isfile(p_filout) and path.getsize(p_filout) > 0 : 
    #    return p_filout
    
    filout = open (p_filout, "w")
    filout.write("PDB ID\tResolution\tNumber of exposed residues\tNumber of residue\tNumber of water\n")
        
    for PDB_ID in l_PDB : 
        print PDB_ID
        d_PDB = loadFile.ExtractInfoPDBID(PDB_ID)
        if d_PDB == {} : 
            continue
            
        if limit_acc == 0.0 : 
            number_residue_exposed = 0
        else : 
            p_fileasa = pr_PDB + PDB_ID + ".asa"
            p_filersa = pr_PDB + PDB_ID + ".rsa"
                
            if not path.isfile(p_fileasa) or not path.isfile(p_filersa) : 
                print "Error NACCESS", p_fileasa
                continue
            else : 
                number_residue_exposed, number_residue = parseNACCESS.numberResExposed(p_filersa, limit_acc)
               
        RX = d_PDB["RX"]
        # case where resolution is not presented in the PDB file
        if RX == 100.0 : 
            continue
        if not "HOH" in d_PDB.keys () : 
            continue
            
        number_residue = len(d_PDB["protein"])
        number_water = len(d_PDB["HOH"])
            
               
        filout.write (str (PDB_ID) + "\t" + str (RX) + "\t" + str (number_residue_exposed) + "\t" + str (number_residue) + "\t" + str (number_water) + "\n")
    filout.close ()
        
    return p_filout




def resolutionByStructure (name_dataset) :
    
    l_structure = structure.ListSub()
    l_path = []
    
    for strut in l_structure :
        l_path.append (pathManage.result(name_dataset) + "water_" + strut + ".dat")
        filout = open (pathManage.result(name_dataset) + "water_" + strut + ".dat", "w") 
        l_file_summary = pathManage.retrieveSummaryFile (strut, name_dataset)
        
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
    rX = parsing.Quality(PDB)
    i_nb_H2O = nbH2O(atom_interest["neighbors"])
    
    return rX, i_nb_H2O, PDB    
    
    
        
def nbH2O (l_atom):
    
    out = 0
    for d_neighbor in l_atom : 
        if d_neighbor["resName"] == "HOH" : 
            out = out + 1
    return out
    
    
        
        
        
  
        
    
    
    
