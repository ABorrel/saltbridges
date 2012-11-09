from os import system
from re import search
from os import listdir
from os import makedirs
from os.path import isfile
from urllib import urlretrieve

import repertory
import loadFile


def formatFilePDB():
    """Manage global PDB database
    manage PDB database, decompress the file and move and rename file with .pdb extension"""

    rep = repertory.openPdbFile()
    listRepertory = listdir(rep)
    
    for subRepertory in listRepertory :
        filesPDB = listdir(rep + subRepertory)   
        for file in filesPDB :
            pathFile = rep + subRepertory +'/'+ file
            cmdDecompress = "gunzip " + pathFile
            system(cmdDecompress)####run decompress
            pathFile = pathFile[0:-3]
            namePDBFile = pathFile[-8:-4] + '.pdb'
            pathFileOut = pathFile[0:-14] + namePDBFile
            cmdMove = "mv " + pathFile + ' ' + pathFileOut
            system(cmdMove)####run move file
    
    for repertoryPDB in listRepertory:
        repertoryPDB = rep + repertoryPDB
        cmdRemove = "rm -r " + repertoryPDB
        system(cmdRemove)####run remove repertory



def appendExtention(rep):
    """Rename global pdb file in repertory with .pdb extention
    in: repertory"""
    
    listFile = listdir(rep)
    for file in listFile : 
        cmd = "mv " + rep + file + " " + rep + file + ".pdb"
        system(cmd)

def retrievePDBFile(fileDataset):
    ligandWithPDB = loadFile.resultFilterPDBLigand(fileDataset)
    repPDBFile = repertory.openPdbFile()
    
    for ligandPDB in ligandWithPDB :
        for pdbName in ligandPDB["PDB"] :
            pdbFilePath = repPDBFile + pdbName + ".pdb"
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

  
  

    