import os.path
import loadFile
import repertory
import structure
from os import makedirs
from os import system
from os import listdir
from os import path
from re import search



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


# -> move in structure
# def atLeastOneClassification(atom):
#     """Classification for at least one sudies
#     in: structure atoms
#     out: classification (string)"""
#     
#     if atom["resName"] == "GLU" or atom["resName"] == "ASP":
#         if atom["name"] == "OE1" or atom["name"] == "OE2" or atom["name"] == "OD1" or atom["name"] == "OD2":
#             return "counterIon"
#         
#         # Amphiprotic
#     if atom["resName"] == "TYR":
#         if atom["name"] == "OH":
#             return "amphiprotic"
# 
#     if atom["resName"] == "CYS":
#         if atom["name"] == "SG":
#             return "amphiprotic"
# 
# 
#     if atom["resName"] == "THR":
#         if atom["name"] == "OG1":
#             return "amphiprotic"
# 
# 
#     if atom["resName"] == "SER":
#         if atom["name"] == "OG":
#             return "amphiprotic"
# 
#     
#     if atom["resName"] == "PHE" or atom["resName"] == "TYR": 
#         if atom["name"] != "CA" : 
#             if atom["name"] != "C" :
#                 return "Carom"
#             
#     if atom["resName"] == "TRP" : 
#         if atom["name"] != "CA" : 
#             if atom["name"] != "CB" : 
#                 if atom["name"] != "CG" :
#                     return "Carom"
#     
#     
#     if atom["resName"] == "HOH":
#         if atom["name"] == "O" : 
#             return "H2O"  
#     
#     return "others"
# 
# 
# def classification (atom): # Milletti 2010
#     """Classification atoms 
#     in: atom
#     out: classification (string)"""
#     
#     listAminoAcid = ["ILE", "LEU", "LYS", "PHE", "TYR", "VAL", "SER", "MET", "ARG", "TRP", "PRO", "GLY", "GLU", "ASN", "HIS", "ALA", "ASP", "GLN", "THR", "CYS"]
#     # Oxygen acid
#     if atom["resName"] == "GLU" or atom["resName"] == "ASP":
#         if atom["name"] == "OE1" or atom["name"] == "OE2" or atom["name"] == "OD1" or atom["name"] == "OD2":
#             return "OxAcid"
# 
#     # Oxygen Donnor/acceptor
#     if atom["resName"] == "TYR":
#         if atom["name"] == "OH":
#             return "Donnor/acceptor"
# 
#     if atom["resName"] == "THR":
#         if atom["name"] == "OG1":
#             return "Donnor/acceptor"
# 
# 
#     if atom["resName"] == "SER":
#         if atom["name"] == "OG":
#             return "Donnor/acceptor"
# 
# 
#     # Sulfure
#     if atom["resName"] == "CYS":
#         if atom["name"] == "SG":
#             return "Sulfur"
# 
# 
#     if atom["resName"] == "THR":
#         if atom["name"] == "OG1":
#             return "amphiprotic"
# 
# 
#     if atom["resName"] == "SER":
#         if atom["name"] == "OG":
#             return "amphiprotic"
#     
#     # Nitrogen basic
#     if atom["resName"] == "HIS" : 
#         if atom["name"] == "NE2" or atom["name"] == "ND1" : 
#             return "Nbasic"
#         
#     if atom["resName"] == "LYS" : 
#         if atom["name"] == "NZ" : 
#             return "Nbasic"
#         
#     if atom["resName"] == "ARG" : 
#         if atom["name"] != "NH1" or  atom["name"] != "NH2" or atom["name"] != "NHE": 
#             return "Nbasic"
#         
#     if atom["resName"] in listAminoAcid : 
#         if atom["name"] == "NXT" : 
#             return "Nbasic"
# 
#     # Nitrogen donnor
#     if atom["resName"] == "ASN" : 
#         if atom["name"] == "ND2" : 
#             return "Ndonnor"
#         
#     if atom["resName"] == "GLN" : 
#         if atom["name"] == "NE2" : 
#             return "Ndonnor"
#             
#     if atom["resName"] in listAminoAcid : 
#         if atom["name"] == "N" : 
#             return "Ndonnor"
#     
#     # Caromatic
#     if atom["resName"] == "PHE" or atom["resName"] == "TYR": 
#         if atom["name"] != "CA" : 
#             if atom["name"] != "C" :
#                 return "Carom"
#             
#     if atom["resName"] == "TRP" : 
#         if atom["name"] != "CA" : 
#             if atom["name"] != "CB" : 
#                 if atom["name"] != "CG" :
#                     return "Carom"
#     
#     if atom["resName"] in listAminoAcid :
#         if atom["name"] == "O" :
#             return "OxAccept" 
#     
#     
#     if atom["resName"] == "ASN" or atom["resName"] == "GLN":
#         if atom["name"] == "OD1" or atom["name"] == "OE1" :
#             return "OxAccept"
#         
#         
#     if atom["resName"] == "HOH":
#         if atom["name"] == "O" : 
#             return "H2O" 
#     
#     
#     return "others"
#     
     
    

##################Serine Protease###########################

def checkListResult(list):
    
    if len(list) == 0 or len(list) == 1 : 
        return
    
    listGlobal = []
    
    for listSubStructure in list : 
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
                del list[j]
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
    


    
    
