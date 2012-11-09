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

    classification = ["Carom", "Ndonnor", "Nbasic", "amphiprotic", "OxAcid", "OxAccept", "H2O","others"]
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
    count["OxAcid"] = 0
    count["amphiprotic"] = 0
    count["Nbasic"] = 0
    count["Ndonnor"] = 0
    count["OxAccept"] = 0
    count["Carom"] = 0
    count["H2O"] = 0
    count["others"] = 0
    
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
    
        for atLeastOne in listAtLeastOneStudy() : 
            count[distance]["atLeastOne"][atLeastOne] = {}
            count[distance]["atLeastOne"][atLeastOne][atLeastOne]  = 0
            count[distance]["atLeastOne"][atLeastOne]["others"] = 0
 
    return count   
    
    
    
def countAtLeastOne():

    listStudy = listStructure()
    listAtLeastOne = listAtLeastOneStudy()
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

    struct = {}
    struct["2.00"] = {}
    struct["2.50"] = {}
    struct["3.00"] = {}
    struct["OUT"] = {}
    struct["NMR"] = {}

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
    for type in count["2.0"]["angleVector"].keys():
        countAngle[type] = {}
        for distance in count.keys():       
            countAngle[type][str(distance)] = {}
            angleTemp = 0
            for keyAngle in listAngle : 
                countAngle[type][str(distance)][str(keyAngle)] = {}
            
                for classe in  count[distance]["angleVector"][type].keys():    
                    countAngle[type][str(distance)][str(keyAngle)][classe] = 0
                    nbAngles = len(count[distance]["angleVector"][type][classe]["angles"])
                    for i in range(0, nbAngles) : 
                        
                        for angleVector in count[str(distance)]["angleVector"][type][classe]["angles"][i] :
                            if angleVector <= keyAngle :
                                if angleVector >= angleTemp :
                                    countAngle[type][str(distance)][str(keyAngle)][classe] = countAngle[type][str(distance)][str(keyAngle)][classe] + 1
                angleTemp = keyAngle

    
    for type in countAngle.keys():
        for keyAngle in countAngle[type]["2.0"].keys() : 
            for classe in countAngle[type]["2.0"][keyAngle].keys() : 
                try : countAngle[type]["2.5"][str(keyAngle)][classe] = countAngle[type]["2.5"][str(keyAngle)][classe] - countAngle[type]["2.0"][str(keyAngle)][classe]
                except : pass
                try : countAngle[type]["3.0"][str(keyAngle)][classe] = countAngle[type]["3.0"][str(keyAngle)][classe] - countAngle[type]["2.5"][str(keyAngle)][classe] - countAngle[type]["2.0"][str(keyAngle)][classe]
                except : pass
                try : countAngle[type]["3.5"][str(keyAngle)][classe] = countAngle[type]["3.5"][str(keyAngle)][classe] - countAngle[type]["3.0"][str(keyAngle)][classe] - countAngle[type]["2.5"][str(keyAngle)][classe] - countAngle[type]["2.0"][str(keyAngle)][classe]
                except : pass
                try : countAngle[type]["4.0"][str(keyAngle)][classe] = countAngle[type]["4.0"][str(keyAngle)][classe] - countAngle[type]["3.5"][str(keyAngle)][classe] - countAngle[type]["3.0"][str(keyAngle)][classe] - countAngle[type]["2.5"][str(keyAngle)][classe] - countAngle[type]["2.0"][str(keyAngle)][classe]
                except : pass
                try : countAngle[type]["4.5"][str(keyAngle)][classe] = countAngle[type]["4.5"][str(keyAngle)][classe] - countAngle[type]["4.0"][str(keyAngle)][classe] - countAngle[type]["3.5"][str(keyAngle)][classe] - countAngle[type]["3.0"][str(keyAngle)][classe] - countAngle[type]["2.5"][str(keyAngle)][classe] - countAngle[type]["2.0"][str(keyAngle)][classe]
                except : pass
                try : countAngle[type]["5.0"][str(keyAngle)][classe] = countAngle[type]["5.0"][str(keyAngle)][classe] - countAngle[type]["4.5"][str(keyAngle)][classe] - countAngle[type]["4.0"][str(keyAngle)][classe] - countAngle[type]["3.5"][str(keyAngle)][classe] - countAngle[type]["3.0"][str(keyAngle)][classe] - countAngle[type]["2.5"][str(keyAngle)][classe] - countAngle[type]["2.0"][str(keyAngle)][classe]
                except : pass
                try : countAngle[type]["5.5"][str(keyAngle)][classe] = countAngle[type]["5.5"][str(keyAngle)][classe] - countAngle[type]["5.0"][str(keyAngle)][classe] - countAngle[type]["4.5"][str(keyAngle)][classe] - countAngle[type]["4.0"][str(keyAngle)][classe] - countAngle[type]["3.5"][str(keyAngle)][classe] - countAngle[type]["3.0"][str(keyAngle)][classe] - countAngle[type]["2.5"][str(keyAngle)][classe] - countAngle[type]["2.0"][str(keyAngle)][classe]
                except : pass
                try : countAngle[type]["6.0"][str(keyAngle)][classe] = countAngle[type]["6.0"][str(keyAngle)][classe] - countAngle[type]["5.5"][str(keyAngle)][classe] - countAngle[type]["5.0"][str(keyAngle)][classe] - countAngle[type]["4.5"][str(keyAngle)][classe] - countAngle[type]["4.0"][str(keyAngle)][classe] - countAngle[type]["3.5"][str(keyAngle)][classe] - countAngle[type]["3.0"][str(keyAngle)][classe] - countAngle[type]["2.5"][str(keyAngle)][classe] - countAngle[type]["2.0"][str(keyAngle)][classe]
                except : pass
                
                
    return countAngle


###############################
#       list structure        #
###############################


def listDistance (distanceMax):
    
    list = []
    distance = 2.0
    while distance <= distanceMax : 
        strDistance = str("%.1f" % distance)
        list.append(strDistance)
        distance = distance + 0.5
    return list


def listStructure ():
    return ["Primary", "Secondary", "Tertiary", "Diamine", "Guanidium","Imidazole","Pyridine"]
    

def listAtLeastOneStudy(): 
    return ["counterIon", "Carom", "amphiprotic", "H2O"]

def listTypeStudy ():
    return ["OxAcid", "amphiprotic", "Nbasic", "Ndonnor", "Carom"]

def listGlobalStudy():
    return ['residue', 'proportionAtom', 'angleVector', 'ligand', 'ResidueAllAtom', 'distanceOx', 'byAA', 'H2O', 'atom', 'proportionType']
    
def selectionAngle():
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
    angleStruct["Imidazole"]["INF"] = 90
    angleStruct["Imidazole"]["SUP"] = 150
    
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
