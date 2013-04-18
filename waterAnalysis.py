"""
BORREL Alexandre
04-2013
"""
from os import listdir, path
from re import search
import parsing
import parseNACCESS
import repertory


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
