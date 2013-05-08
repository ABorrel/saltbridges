# classif MILLETTI
# dico_atom_C = {"R":["CB", "CD", "CG"], "M":["CB"], "F":["CB"], "L:":["CB","CD1", "CD2", "CG"], "W":["CB"], "D":["CB"], "K":["CB", "CE"], "H":["CB"], "V":["CB", "CG1", "CG2"], "Q":["CB", "CG", "CD"], "A":["CB"], "E":["CB"], "P":["CB", "CD", "CG"], "C":["CB"], "Y":["CB"], "N":["CB"], "I":["CB", "CD1", "CG1", "CG2"]}
# dico_atom_Car = {"F":["CD1", "CD2", "CE1", "CE2", "CG", "CZ"], "W":["CD1", "CD2", "CE2", "CE3", "CG", "CH2", "CZ2", "CZ3"], "H":["CD2", "CE1"], "Y":["CD1", "CD2", "CE1", "CE2", "CG", "CZ"]}
# dico_atom_Carg = {"S":["CB"], "T":["CB"]}
# dico_atom_N = {"A":["N"], "C":["N"],"D":["N"],"E":["N"],"F":["N"],"G":["N"],"H":["N"],"I":["N"],"K":["N"],"L":["N"],"M":["N"],"N":["N", "ND2"],"P":["N"],"Q":["N", "NE2"],"R":["N"],"S":["N"],"T":["N"],"V":["N"],"W":["N"],"Y":["N"]}
# dico_atom_ND1 = {"H":["ND1"]}
# dico_atom_NE2 = {"H":["NE2"]}
# dico_atom_Nlys = {"K":["NZ"]}
# dico_atom_Ntrp = {"W":["NE1"]}
# dico_atom_O = {"A":["O"], "C":["O"],"D":["O"],"E":["O"],"F":["O"],"G":["O"],"H":["O"],"I":["O"],"K":["O"],"L":["O"],"M":["O"],"N":["O"],"P":["O"],"Q":["O"],"R":["O"],"S":["O"],"T":["O"],"V":["O"],"W":["O"],"Y":["O"]}
# dico_atom_Ocoo = {"E":["CG"]}
# dico_atom_Ooh = {"S":["CB"], "T":["CB"]}
# dico_atom_Otyr = {"Y":["OH"]}
# dico_atom_S = {"C":["SG"]}
# dico_atom_Ccoo = {"D":["CB"], "E":["CG"]}
# dico_atom_Cgln = {"Q":["CG"], "N":["CB"]}
# dico_atom_Hyd = {"R":["CZ"], "M":["SD"], "F":["CG", "CZ"], "L":["CG"], "W":["CE3", "CG", "CZ2"], "H":["CG"], "V":["CB"], "P":["CG"], "C":["SG"], "Y":["CG", "CZ"]}


def classificationATOM (atom, out_list = 0):
    """Classification atoms 
    in: atom
    out: classification (string)"""
    
    list_classif = ["OxAcid", "ODonAcc", "Sulfur", "Nhis", "Nbasic", "Carom", "OxPep", "Ndonnor","OxAccept" , "H2O", "CPep", "others"]
    if out_list : 
        return list_classif
    
    
    listAminoAcid = ["ILE", "LEU", "LYS", "PHE", "TYR", "VAL", "SER", "MET", "ARG", "TRP", "PRO", "GLY", "GLU", "ASN", "HIS", "ALA", "ASP", "GLN", "THR", "CYS"]
    # Oxygen acid
    if atom["resName"] == "GLU" or atom["resName"] == "ASP":
        if atom["name"] == "OE1" or atom["name"] == "OE2" or atom["name"] == "OD1" or atom["name"] == "OD2":
            return "OxAcid"

    # Oxygen Donnor/acceptor
    if atom["resName"] == "TYR":
        if atom["name"] == "OH":
            return "ODonAcc"

    if atom["resName"] == "THR":
        if atom["name"] == "OG1":
            return "ODonAcc"


    if atom["resName"] == "SER":
        if atom["name"] == "OG":
            return "ODonAcc"


    # Sulfure
    if atom["resName"] == "CYS":
        if atom["name"] == "SG":
            return "Sulfur"


    # Nitrogen histidine
    if atom["resName"] == "HIS" : 
        if atom["name"] == "NE2" or atom["name"] == "ND1" : 
            return "Nhis"
        
    # Nitrogen basic        
    if atom["resName"] == "LYS" : 
        if atom["name"] == "NZ" : 
            return "Nbasic"
        
    if atom["resName"] == "ARG" : 
        if atom["name"] != "NH1" or  atom["name"] != "NH2" or atom["name"] != "NHE": 
            return "Nbasic"
        
#     if atom["resName"] in listAminoAcid :  ?????
#         if atom["name"] == "NXT" : 
#             return "Nbasic"

    # Nitrogen donnor
    if atom["resName"] == "ASN" : 
        if atom["name"] == "ND2" : 
            return "Ndonnor"
        
    if atom["resName"] == "GLN" : 
        if atom["name"] == "NE2" : 
            return "Ndonnor"
            
    if atom["resName"] in listAminoAcid : 
        if atom["name"] == "N" : 
            return "Ndonnor"
    
    # Caromatic
    if atom["resName"] == "PHE" or atom["resName"] == "TYR": 
        if atom["name"] != "CA" : 
            if atom["name"] != "C" :
                return "Carom"
            
    if atom["resName"] == "TRP" : 
        if atom["name"] != "CA" : 
            if atom["name"] != "CB" : 
                return "Carom"
    
    # O peptitique
    if atom["resName"] in listAminoAcid :
        if atom["name"] == "O" :
            return "OxPep" 
    
    # O acid carboxylic
    if atom["resName"] == "ASN" or atom["resName"] == "GLN":
        if atom["name"] == "OD1" or atom["name"] == "OE1" :
            return "OxAccept"
        
    # water
    if atom["resName"] == "HOH":
        if atom["name"] == "O" : 
            return "H2O" 
    
    # C peptidique
    if atom["resName"] in listAminoAcid :
        if atom["name"] == "C" :
            return "CPep"     
    
    
    return "others"




    

def ligandPDB():
    struct = {}
    struct["name"] = ""
    struct["ligands"] = []

    return struct


def countInstanceDataSet():

    count = {}
    count["name"] = ""
    count["Number PDB"] = 0

    return count


def amine():

    listAmine = {}
    
    for element in listStructure(): 
        listAmine[element] = []

    return listAmine


def countAngle () :

    classification = classificationATOM("", out_list=1)
    count = {}
    for element in listStructure () : 
        count[element] = {}
        for classe in classification :
            count[element][classe] = {}
            count[element][classe]["angles"] = []
            count[element][classe]["distance"] = []

    return count


def countOx ():

    count = {}
    for element in listStructure() : 
        count[element] = []
    
    return count


def countAtom():

    count = {}
    for element in listStructure() : 
        count[element] = {}
    
    return count


def countResidue():

    listAminoAcid = ["ILE", "LEU", "LYS", "PHE", "TYR", "VAL", "SER", "MET", "ARG", "TRP", "PRO", "GLY", "GLU", "ASN", "HIS", "ALA", "ASP", "GLN", "THR", "CYS", "HOH"]
    
    count = {}
    for struct in listStructure() : 
        count[struct] = {}
        for aminoAcid in listAminoAcid:
            count[struct][aminoAcid] = {}
            count[struct][aminoAcid]["main"] = 0
            count[struct][aminoAcid]["side"] = 0
    return count

def countLigand():

    listS = listStructure()
    count = {}
    for struct in listS :
        count[struct] = {} 

    return count

def countbyAminoAcid():

    listAminoAcid = ["ILE", "LEU", "LYS", "PHE", "TYR", "VAL", "SER", "MET", "ARG", "TRP", "PRO", "GLY", "GLU", "ASN", "HIS", "ALA", "ASP", "GLN", "THR", "CYS"]
    listStruct = listStructure()
    count = {}
    for element in listStruct : 
        count[element] = {}
    count["global"] = {}

    for key in count.keys():
        for aminoAcid in listAminoAcid:
            count[key][aminoAcid] = {}
    return count


def countGlobalAmine(distanceMax):

    listDis = listDistance(distanceMax)
    count = {}
    
    for distance in listDis : 
        count[distance] = countType()
    return count


def countType():

    count = {}
    count["ligand"] = countLigand()
    count["angle"] = countAngle()
    count["atom"] = countAtom()
    count["residue"] = countResidue()
    count["byAA"] = countbyAminoAcid()
    count["distanceOx"] = countOx()
    count["atLeastOne"] = countAtLeastOne()
    count["proportionAtom"] = countProportionAtom()
    count["proportionType"] = countProportionType()
    count["ResidueAllAtom"] = countResidueGlobal()
    return count


def countProportionAtom ():
    
    count = {}
    listStruct = listStructure()
    for element in listStruct :
        count[element] = {}
    count["Global"] = {}
    count["GlobalAmine"] = {}
    return count
    
  

def countProportionType ():

    listS = listStructure()
    count = {}
    for struct in listS : 
        count[struct] = {}
        count[struct]["allNumberNeighbors"] = countClassificationAtoms()
        
    count["Global"] = {}
    count["Global"]["allNumberNeighbors"] = countClassificationAtoms()
    count["GlobalAmine"] = {}
    return count


def countElements():

    count = {}
    count["O"] = 0
    count["C"] = 0
    count["N"] = 0
    count["S"] = 0
    count["others"] = 0
    return count

def countClassificationAtoms():

    count = {}
    l_classe = classificationATOM("", out_list = 1)
    for classe in l_classe : 
        count[classe] = 0
    
    return count


def listAminoAcidCheck():
    listOut = {}
    listAminoAcid = ["ILE", "LEU", "LYS", "PHE", "TYR", "VAL", "SER", "MET", "ARG", "TRP", "PRO", "GLY", "GLU", "ASN", "HIS", "ALA", "ASP", "GLN", "THR", "CYS", "HOH"]

    for aminoAcid in listAminoAcid :
        listOut[aminoAcid] = {}
        listOut[aminoAcid]["side"] = 0
        listOut[aminoAcid]["main"] = 0

    return listOut



def countAtLeastOneGlobalStruct(distanceMax):
    
    listDis = listDistance(distanceMax)
    count = {}

    for distance in listDis : 
        count[distance] = {}
        count[distance]["atLeastOne"] = {}
    
        for atLeastOne in  classificationATOM("", out_list=1) : 
            count[distance]["atLeastOne"][atLeastOne] = {}
            count[distance]["atLeastOne"][atLeastOne][atLeastOne] = 0
            count[distance]["atLeastOne"][atLeastOne]["others"] = 0
 
    return count   
    
    
    
def countAtLeastOne():

    listStudy = listStructure()
    listAtLeastOne = classificationATOM("", out_list=1)
    count = {}

    for atLeastOne in listAtLeastOne : 
        count[atLeastOne] = {}
        for studyStruct in listStudy : 
            count[atLeastOne][studyStruct] = {}
            count[atLeastOne][studyStruct][atLeastOne] = 0
            count[atLeastOne][studyStruct]["others"] = 0
        
    return count



def countResidueGlobal():

    listAminoAcid = ["ILE", "LEU", "LYS", "PHE", "TYR", "VAL", "SER", "MET", "ARG", "TRP", "PRO", "GLY", "GLU", "ASN", "HIS", "ALA", "ASP", "GLN", "THR", "CYS", "HOH"]

    count = {}
    for aa in listAminoAcid :
        count[aa] = {}
        count[aa]["main"] = 0
        count[aa]["side"] = 0
    return count


def resolutionFilter():
    """
    Result of resolution filter
    """
    
    struct = {}
    struct["1.50"] = {}
    struct["2.00"] = {}
    struct["2.50"] = {}
    struct["3.00"] = {}
    struct["OUT"] = {}
    struct["NMR"] = {} # remove

    return struct


def countGroupDataset() : 

    listS = listStructure()    
    struct = {}
    for element in listS : 
        struct[element] = 0
    
    return struct


def countAngleType(count):
    
    listAngle = range(0, 190, 10)
    countAngle = {}
    print count["2.0"].keys()
    for type_strucutre in count["2.0"]["angle"].keys():
        countAngle[type_strucutre] = {}
        for distance in count.keys():       
            countAngle[type_strucutre][str(distance)] = {}
            angleTemp = 0
            for keyAngle in listAngle : 
                countAngle[type_strucutre][str(distance)][str(keyAngle)] = {}
            
                for classe in  count[distance]["angle"][type_strucutre].keys():    
                    countAngle[type_strucutre][str(distance)][str(keyAngle)][classe] = 0
                    nbAngles = len(count[distance]["angle"][type_strucutre][classe]["angles"])
                    for i in range(0, nbAngles) : 
                    # restrint angle position
                        for angleVector in count[str(distance)]["angle"][type_strucutre][classe]["angles"][i] :
                            if angleVector <= keyAngle :
                                if angleVector >= angleTemp :
                                    countAngle[type_strucutre][str(distance)][str(keyAngle)][classe] = countAngle[type_strucutre][str(distance)][str(keyAngle)][classe] + 1
                angleTemp = keyAngle

    
    for type_strucutre in countAngle.keys():
        for keyAngle in countAngle[type_strucutre]["2.0"].keys() : 
            for classe in countAngle[type_strucutre]["2.0"][keyAngle].keys() : 
                try : countAngle[type_strucutre]["2.5"][str(keyAngle)][classe] = countAngle[type_strucutre]["2.5"][str(keyAngle)][classe] - countAngle[type_strucutre]["2.0"][str(keyAngle)][classe]
                except : pass
                try : countAngle[type_strucutre]["3.0"][str(keyAngle)][classe] = countAngle[type_strucutre]["3.0"][str(keyAngle)][classe] - countAngle[type_strucutre]["2.5"][str(keyAngle)][classe] - countAngle[type_strucutre]["2.0"][str(keyAngle)][classe]
                except : pass
                try : countAngle[type_strucutre]["3.5"][str(keyAngle)][classe] = countAngle[type_strucutre]["3.5"][str(keyAngle)][classe] - countAngle[type_strucutre]["3.0"][str(keyAngle)][classe] - countAngle[type_strucutre]["2.5"][str(keyAngle)][classe] - countAngle[type_strucutre]["2.0"][str(keyAngle)][classe]
                except : pass
                try : countAngle[type_strucutre]["4.0"][str(keyAngle)][classe] = countAngle[type_strucutre]["4.0"][str(keyAngle)][classe] - countAngle[type_strucutre]["3.5"][str(keyAngle)][classe] - countAngle[type_strucutre]["3.0"][str(keyAngle)][classe] - countAngle[type_strucutre]["2.5"][str(keyAngle)][classe] - countAngle[type_strucutre]["2.0"][str(keyAngle)][classe]
                except : pass
                try : countAngle[type_strucutre]["4.5"][str(keyAngle)][classe] = countAngle[type_strucutre]["4.5"][str(keyAngle)][classe] - countAngle[type_strucutre]["4.0"][str(keyAngle)][classe] - countAngle[type_strucutre]["3.5"][str(keyAngle)][classe] - countAngle[type_strucutre]["3.0"][str(keyAngle)][classe] - countAngle[type_strucutre]["2.5"][str(keyAngle)][classe] - countAngle[type_strucutre]["2.0"][str(keyAngle)][classe]
                except : pass
                try : countAngle[type_strucutre]["5.0"][str(keyAngle)][classe] = countAngle[type_strucutre]["5.0"][str(keyAngle)][classe] - countAngle[type_strucutre]["4.5"][str(keyAngle)][classe] - countAngle[type_strucutre]["4.0"][str(keyAngle)][classe] - countAngle[type_strucutre]["3.5"][str(keyAngle)][classe] - countAngle[type_strucutre]["3.0"][str(keyAngle)][classe] - countAngle[type_strucutre]["2.5"][str(keyAngle)][classe] - countAngle[type_strucutre]["2.0"][str(keyAngle)][classe]
                except : pass
                try : countAngle[type_strucutre]["5.5"][str(keyAngle)][classe] = countAngle[type_strucutre]["5.5"][str(keyAngle)][classe] - countAngle[type_strucutre]["5.0"][str(keyAngle)][classe] - countAngle[type_strucutre]["4.5"][str(keyAngle)][classe] - countAngle[type_strucutre]["4.0"][str(keyAngle)][classe] - countAngle[type_strucutre]["3.5"][str(keyAngle)][classe] - countAngle[type_strucutre]["3.0"][str(keyAngle)][classe] - countAngle[type_strucutre]["2.5"][str(keyAngle)][classe] - countAngle[type_strucutre]["2.0"][str(keyAngle)][classe]
                except : pass
                try : countAngle[type_strucutre]["6.0"][str(keyAngle)][classe] = countAngle[type_strucutre]["6.0"][str(keyAngle)][classe] - countAngle[type_strucutre]["5.5"][str(keyAngle)][classe] - countAngle[type_strucutre]["5.0"][str(keyAngle)][classe] - countAngle[type_strucutre]["4.5"][str(keyAngle)][classe] - countAngle[type_strucutre]["4.0"][str(keyAngle)][classe] - countAngle[type_strucutre]["3.5"][str(keyAngle)][classe] - countAngle[type_strucutre]["3.0"][str(keyAngle)][classe] - countAngle[type_strucutre]["2.5"][str(keyAngle)][classe] - countAngle[type_strucutre]["2.0"][str(keyAngle)][classe]
                except : pass
                
                
    return countAngle


###############################
#       list structure        #
###############################


def listDistance (distanceMax):
    
    list_out = []
    distance = 2.0
    while distance <= distanceMax : 
        strDistance = str("%.1f" % distance)
        list_out.append(strDistance)
        distance = distance + 0.5
    return list_out


def listStructure ():
    return ["Primary", "Secondary", "Tertiary", "Diamine", "Guanidium", "Imidazole", "Pyridine", "AcidCarboxylic"]
    

# def listAtLeastOneStudy(): 
#     return ["counterIon", "Carom", "amphiprotic", "H2O"]

# def listTypeStudy ():
#     return ["OxAcid", "amphiprotic", "Nbasic", "Ndonnor", "Carom"]

def listGlobalStudy():
    return ['residue', 'proportionAtom', 'angleVector', 'ligand', 'ResidueAllAtom', 'distanceOx', 'byAA', 'H2O', 'atom', 'proportionType']
    
def selectionAngle():
    """
    Impose angle in different structure
    """
    angleStruct = {}
    
    angleStruct["Primary"] = {}
    angleStruct["Primary"]["INF"] = 90
    angleStruct["Primary"]["SUP"] = 150
    
    angleStruct["Secondary"] = {}
    angleStruct["Secondary"]["INF"] = 90
    angleStruct["Secondary"]["SUP"] = 150
    
    angleStruct["Tertiary"] = {}
    angleStruct["Tertiary"]["INF"] = 90
    angleStruct["Tertiary"]["SUP"] = 150
    
    angleStruct["Imidazole"] = {}
    angleStruct["Imidazole"]["INF"] = 0
    angleStruct["Imidazole"]["SUP"] = 30
    
    angleStruct["Guanidium"] = {}
    angleStruct["Guanidium"]["INF"] = 90
    angleStruct["Guanidium"]["SUP"] = 150
    
    angleStruct["Diamine"] = {}
    angleStruct["Diamine"]["INF"] = 90
    angleStruct["Diamine"]["SUP"] = 150
    
    angleStruct["Pyridine"] = {}
    angleStruct["Pyridine"]["INF"] = 90
    angleStruct["Pyridine"]["SUP"] = 150
    
    return angleStruct
