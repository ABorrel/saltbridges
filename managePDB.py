from os import system, listdir, makedirs
from re import search, compile, findall
from os.path import isfile
from urllib import urlretrieve
from shutil import copy

import loadFile
import repertory
import formatCharacter

def formatFilePDB(dir_PDB):
    """Manage global PDB database
    manage PDB database, decompress the file and move and rename file with .pdb extension"""
    listRepertory = listdir(dir_PDB)
    
    for subRepertory in listRepertory :
        filesPDB = listdir(dir_PDB + subRepertory)   
        for file in filesPDB :
            pathFile = dir_PDB + subRepertory + '/' + file
            cmdDecompress = "gunzip " + pathFile
            system(cmdDecompress)  ####run decompress
            pathFile = pathFile[0:-3]
            namePDBFile = pathFile[-8:-4] + '.pdb'
            pathFileOut = pathFile[0:-14] + namePDBFile
            cmdMove = "mv " + pathFile + ' ' + pathFileOut
            system(cmdMove)  ####run move file
    
    for repertoryPDB in listRepertory:
        repertoryPDB = dir_PDB + repertoryPDB
        cmdRemove = "rm -r " + repertoryPDB
        system(cmdRemove)  ####run remove repertory



def appendExtention(rep):
    """Rename global pdb file in repertory with .pdb extention
    in: repertory"""
    
    listFile = listdir(rep)
    for file in listFile : 
        cmd = "mv " + rep + file + " " + rep + file + ".pdb"
        system(cmd)

def retrievePDBFile(fileDataset, dir_PDB):
    ligandWithPDB = loadFile.resultFilterPDBLigand(fileDataset)
    
    for ligandPDB in ligandWithPDB :
        for pdbName in ligandPDB["PDB"] :
            pdbFilePath = dir_PDB + pdbName + ".pdb"
            if not isfile(pdbFilePath) :
                adressePDB = "http://www.pdb.org/pdb/files/" + pdbName + ".pdb"
                try:
                    pathFilePDB = urlretrieve(adressePDB)
                    print pathFilePDB
        
                    cmd = "mv " + pathFilePDB[0] + " " + pdbFilePath
                    print cmd
                    system (cmd)
        
                except:
                    print "Impossible retrieve PDB file"

  

def retrievePDB (path_folder_database, path_file_withPDBid):
    
    filin_PDB = open (path_file_withPDBid, "r")
    file_read = filin_PDB.read()
    filin_PDB.close ()
    
    path_new_folder_database = path_file_withPDBid.split (".")[0] + "/"
    
    try: makedirs(path_folder_database, mode=0777)
    except: pass
    
    regex = compile ("[0-9A-Za-z]{4}")
    list_pdb = regex.findall (file_read)
    list_pdb = list(set(list_pdb)) 
    
    for PDB_ID in list_pdb : 
        path_file_PDB = path_folder_database + PDB_ID.lower() + ".pdb"
        path_new = path_new_folder_database + PDB_ID.lower() + ".pdb"
        print path_file_PDB, path_new
        try : copy(path_file_PDB, path_new)
        except : pass
        
        
        


def retriveListPDB (name_database):
    
    if name_database == "PDB" : 
        directory_PDB = repertory.repInit + "PDB/"
        list_files = listdir(directory_PDB)
        
        list_pdb = []
        for file_PDB in list_files : 
            if file_PDB[-4:] == ".pdb" : 
                list_pdb.append (file_PDB[-8:-4])
        
        return formatCharacter.lowerList(list_pdb)
    
    else : 
        path_file_database = repertory.repInit + name_database + ".dat"
        file_database = open (path_file_database, "r")
        file_read = file_database.read()
        file_database.close ()
        
        regex = compile ("[0-9A-Za-z]{4}")
        list_pdb = regex.findall (file_read)
        list_pdb = list(set(list_pdb))
        
        return formatCharacter.lowerList(list_pdb) 
    
    
           
    






  

    
