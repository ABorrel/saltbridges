import os.path
import loadFile
import repertory
import structure
from os import makedirs
from os import system
from os import listdir
from os import path
from re import search
from copy import deepcopy

import statistic


def transfomAA (aa_in):
    list_aa = {"ILE":"I", "LEU":"L", "LYS":"K", "PHE":"F", "TYR":"Y", "VAL":"V", "SER":"S", "MET":"M", "ARG":"R", "TRP":"W" , "PRO":"P", "GLY":"G", "GLU":"E", "ASN":"N", "HIS":"H", "ALA":"A", "ASP":"D", "GLN":"Q", "THR":"T", "CYS":"C" }

    if len (aa_in) == 3 : 
        return list_aa[aa_in]
    elif len (aa_in) == 1 : 
        for aa in list_aa.keys ():
            if list_aa[aa] == aa_in : 
                return aa
    
    return "ERROR"
        
    



#########file management#########

def manageligandsliste(lines, position):
    """reads lines and parse names
    returns the name"""
    
    line = lines[position]
    name = line.split(".")[0]
    return name




def retrievePositionList (list, keyName):
    """retieve in list ligand, dictionnary (keys->name ligand), the index of name ligand
    in: list ligand
    out: index (integers), position in list 
    """
    
    lenList = len(list)
    for position in range(0, lenList):
        if list[position]["name"] == keyName:
            return position


def countNumberPDB(file):
    """Count the number of PDB file that contains the ligand
    in: file result search ligand in PDB files
    out: number of PDB files"""
    
    ligandPDB = loadFile.resultLigandPDB(file)

    nbPDB = 0
    for ligand in ligandPDB:
        nbPDB = nbPDB + len(ligand["PDB"])

    return nbPDB





def existFile(pathFile):
    """test if file exist and if size is superior at 0
    in: path file
    out: boolean"""
    if os.path.isfile(pathFile):
        return 1
    else:
        return 0



def atomInLigands (atomSerial, listAtomLigand):
    """Check if atom is in listAtomLigand return position in listAtomLigand
    in: serial atom checked
    out: integer or -1 if atom is not in list atoms"""

    lenListligand = len(listAtomLigand)
    for i in range(0, lenListligand):
        if listAtomLigand[i]["number"] == atomSerial:
            return i
    return -1



def dellH(listElement):
    """Deleted in list of element the hydrogen
    in: list of element
    out: modification list"""
    
    lenList = len(listElement)

    i = 0
    while i < lenList:
        if listElement[i] == "H":
            del listElement[i]
            lenList = lenList - 1
        else:
            i = i + 1


def atomInList(listAtom, atom):
    """Check if atom is in list atoms
    in: list of atoms, atom cheked
    out: boolean"""
    
    for element in listAtom:
        if element["serial"] == atom["serial"]:
            return 1

    return 0


def checkListResult(l):
    
    if len(l) == 0 or len(l) == 1 : 
        return
    
    listGlobal = []
    
    for listSubStructure in l : 
        listSerial = []
        for atom in listSubStructure : 
            listSerial.append(atom["serial"])
        listSerial.sort()
        listGlobal.append(listSerial)
    
    nbList = len(listGlobal)
    
    
    
    i = 0
    while i < nbList :
#        print nbList 
        j = i + 1
        while j < nbList : 
            if listGlobal[i] == listGlobal[j]:

                del listGlobal[j]
                del l[j]
                nbList = nbList - 1
            else : 
                j = j + 1
            
        i = i + 1



def searchLigandInDataSetFile (datasetFile, nameLigand):
    """Search ligand in dataset file (out dataset construction)
    in: dataset file, name of ligand
    out: index of ligand or 0"""
    
    listLigandsInPDB = loadFile.resultFilterPDBLigand(datasetFile)
    nbLigand = len(listLigandsInPDB)
    i = 0
    while i < nbLigand :
        if listLigandsInPDB[i]["name"] == nameLigand : 
            return i
        
        i = i + 1
        
    return 0
        
        
def searchElementInList(list_element, stringSearch):        
    """Search element in list of string
    in: list, string search
    out: index or 0"""
    
    nbLigand = len(list_element)
    
    i = 0
    while i < nbLigand : 
        if list[i] == stringSearch : 
            return i
        i = i + 1 
    
    return 0
         
     
    
def checkFileEmpty(pathFile):
    """Check if file is empty
    in: path file
    out: boolean"""
    
    if path.isfile(pathFile) != True : 
        return 0
    
    if path.getsize(pathFile) == 0 : 
        return 1
    else : 
        return 0
    

def sumDict(dict_count) : 
    out = 0.0001
    for key in dict_count.keys() :
        if  type (dict_count[key]) is int or type (dict_count[key]) is float :  
            out = out + dict_count[key]
    return out
    
    


def colorAtomType (l_superimpose) : 
    
    for atom in l_superimpose : 
        c_atom = structure.classificationATOM(atom)
        if c_atom == "OxAcid" : 
            atom ["tempFactor"] = 100 - 8.33
        elif c_atom == "ODonAcc" : 
            atom ["tempFactor"] = 100 - (2* 8.33)
        elif c_atom ==  "OxAccept": 
            atom ["tempFactor"] = 100 - (3* 8.33) 
        elif c_atom == "OxPep"  : 
            atom ["tempFactor"] = 100 - (4* 8.33)  
        elif c_atom == "Sulfur"  : 
            atom ["tempFactor"] = 100 - (5* 8.33)
        elif c_atom == "Carom"  : 
            atom ["tempFactor"] = 100 - (6* 8.33)
        elif c_atom == "Carom"  : 
            atom ["tempFactor"] = 100 - (7* 8.33)  
        elif c_atom == "Ndonnor"  : 
            atom ["tempFactor"] = 100 - (8* 8.33)
        elif c_atom == "Nhis"  : 
            atom ["tempFactor"] = 100 - (9* 8.33)
        elif c_atom == "Nbasic"  : 
            atom ["tempFactor"] = 100 - (9* 8.33)   
        elif c_atom == "CPep"  : 
            atom ["tempFactor"] = 100 - (10* 8.33) 
        elif c_atom == "others"  : 
            atom ["tempFactor"] = 100 - (11* 8.33) 
        elif c_atom == "H2O"  : 
            atom ["tempFactor"] = 100 - (12* 8.33)    

    
       
    
