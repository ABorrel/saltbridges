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

l_res = ["ALA", "ILE", "LEU", "VAL", "MET", "CYS", "PHE", "TRP", "TYR", "HIS", "THR", "SER", "ASN", "GLN", "ASP", "GLU", "ARG", "LYS", "PRO", "GLY"]


def classificationATOM (atom="", out_list = 0):
    """Classification atoms 
    in: atom
    out: classification (string)"""
    
#     l_classif = ["Oox", "Cox", "Oh", "Oc", "Ow", "Nam", "Nim", "Ngu", "Cgu", "NaI", "Car", "Xot"]
    l_classif = ["Oox", "Oh", "Oc", "Ow", "Nam", "Nim", "Ngu", "NaI", "Car", "Xot"]
    if out_list : 
        return l_classif
    
    # Oxygen acid
    if atom["resName"] == "GLU" or atom["resName"] == "ASP":
        if atom["name"] == "OE1" or atom["name"] == "OE2" or atom["name"] == "OD1" or atom["name"] == "OD2":
            return "Oox"


#     # cabone of COO
#     if atom["resName"] == "GLU" : 
#         if atom["name"] == "CD"  : 
#             return "Cox"
        
#     if atom["resName"] == "ASP" :
#         if atom["name"] == "CG"  : 
#             return "Cox"
        

    # Oxygen Donnor/acceptor
    if atom["resName"] == "TYR":
        if atom["name"] == "OH":
            return "Oh"

    if atom["resName"] == "THR":
        if atom["name"] == "OG1":
            return "Oh"


    if atom["resName"] == "SER":
        if atom["name"] == "OG":
            return "Oh"


    # Nitrogen histidine
    if atom["resName"] == "HIS" : 
        if atom["name"] == "NE2" or atom["name"] == "ND1" : 
            return "Nim"
        
    # Nitrogen basic        
    if atom["resName"] == "LYS" : 
        if atom["name"] == "NZ" : 
            return "NaI"
        
    if atom["resName"] == "ARG" : 
        if atom["name"] == "NH1" or  atom["name"] == "NH2" or atom["name"] == "NHE" or  atom["name"] == "NE": 
            return "Ngu"
      
#     if atom["resName"] == "ARG" : 
#         if atom["name"] == "CZ" : 
#             return "Cgu"
        
    # Nitrogen donor
    if atom["resName"] == "ASN" : 
        if atom["name"] == "ND2" : 
            return "Nam"
        
    if atom["resName"] == "GLN" : 
        if atom["name"] == "NE2" : 
            return "Nam"
            
    if atom["resName"] in l_res : 
        if atom["name"] == "N" : 
            return "Nam"
    
    # Caromatic
    if atom["resName"] == "PHE" or atom["resName"] == "TYR": 
        if atom["name"] != "CA" : 
            if atom["name"] != "C" :
                return "Car"
            
    if atom["resName"] == "TRP" : 
        if atom["name"] != "CA" : 
            if atom["name"] != "CB" : 
                return "Car"
    
    # O peptitique
    if atom["resName"] in l_res :
        if atom["name"] == "O" :
            return "Oc" 
    
    # O acid carboxylic
    if atom["resName"] == "ASN" or atom["resName"] == "GLN":
        if atom["name"] == "OD1" or atom["name"] == "OE1" :
            return "Oc"
        
    # water
    if atom["resName"] == "HOH":
        if atom["name"] == "O" : 
            return "Ow" 
    
    # C peptidique
    if atom["resName"] in l_res :
        if atom["name"] == "C" :
            return "Oc"     
    
    

    
#     print atom["resName"], atom["name"]
    
    return "Xot"




    

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



def dAngleType () :

    classification = classificationATOM("", out_list=1)
    count = {}
    for element in ListSub () : 
        count[element] = {}
        for classe in classification :
            count[element][classe] = {}
            count[element][classe]["angles"] = []
            count[element][classe]["distance"] = []

    return count


def countOx ():

    count = {}
    for element in ListSub() : 
        count[element] = []
    
    return count


def countAtom():

    count = {}
    for element in ListSub() : 
        count[element] = {}
    
    return count


def countResidue():

    listAminoAcid = ["ILE", "LEU", "LYS", "PHE", "TYR", "VAL", "SER", "MET", "ARG", "TRP", "PRO", "GLY", "GLU", "ASN", "HIS", "ALA", "ASP", "GLN", "THR", "CYS", "HOH"]
    
    count = {}
    for struct in ListSub() : 
        count[struct] = {}
        for aminoAcid in listAminoAcid:
            count[struct][aminoAcid] = {}
            count[struct][aminoAcid]["main"] = 0
            count[struct][aminoAcid]["side"] = 0
    return count

def countLigand():

    listS = ListSub()
    count = {}
    for struct in listS :
        count[struct] = {} 

    return count

def countbyAminoAcid():

    listAminoAcid = ["ILE", "LEU", "LYS", "PHE", "TYR", "VAL", "SER", "MET", "ARG", "TRP", "PRO", "GLY", "GLU", "ASN", "HIS", "ALA", "ASP", "GLN", "THR", "CYS"]
    listStruct = ListSub()
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
    count["angle"] = dAngleType()
    count["atom"] = countAtom()
    count["residue"] = countResidue()
    count["byAA"] = countbyAminoAcid()
    count["distanceOx"] = countOx()
    count["atLeastOne"] = {}
    count["threeAnalysis"] = countThreeNeigbors ()
    count["proportionAtom"] = countProportionAtom()
    count["proportionType"] = countProportionType()
    count["ResidueAllAtom"] = countResidueGlobal()
    count ["numberNeighbors"] = {}
    return count


def countThreeNeigbors () : 
    
    d_out = {}
    l_subs = ListSub()
    l_subs.append("global")
    
    # MAJ 29-11-2013 -> nb_neighbor considering alway 8
    for sub_struct in l_subs : 
        
        if sub_struct == "I" : 
            nb_n = 8
        elif sub_struct == "III" : 
            nb_n = 8
        elif sub_struct == "IMD" or sub_struct == "II": 
            nb_n = 8
        else : 
            nb_n = 8
            
        d_out[sub_struct] = {}
        d_out[sub_struct]["angle1_2"]  = []
        d_out[sub_struct]["angle1_3"]  = []
        d_out[sub_struct]["angle2_3"]  = []
        for i in range(1,nb_n +1) : 
            d_out[sub_struct][i] = {} 
            d_out[sub_struct][i]["distance"] = []
            d_out[sub_struct][i]["classe"] = []
        l_classe = classificationATOM("", out_list = 1)
        for classe in l_classe : 
            for i in range(1,nb_n +1) : 
                d_out[sub_struct][i][classe] = 0.0
    
    return d_out



def countProportionAtom ():
    
    count = {}
    listStruct = ListSub()
    for element in listStruct :
        count[element] = {}
    count["Global"] = {}
    count["GlobalAmine"] = {}
    return count
    
  

def countProportionType ():

    listS = ListSub()
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

def EmptyListClassifAtom():
    
    d_out = {}
    l_classe = classificationATOM("", out_list = 1)
    for classe in l_classe : 
        d_out[classe] = []
    
    return d_out
    


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
    
    return count   
    
    
    

def countResidueGlobal():

    l_aa = ["ILE", "LEU", "LYS", "PHE", "TYR", "VAL", "SER", "MET", "ARG", "TRP", "PRO", "GLY", "GLU", "ASN", "HIS", "ALA", "ASP", "GLN", "THR", "CYS", "HOH"]

    count = {}
    for aa in l_aa :
        count[aa] = {}
        count[aa]["main"] = 0
        count[aa]["side"] = 0
    return count


def resolutionFilter():
    """
    Result of resolution filter
    -> reduce to only two datasets
    - RX 3.0 and 1.5
    """
    
    d_out = {}
    d_out["1.50"] = {}
    #d_out["2.00"] = {}
    #d_out["2.50"] = {}
    d_out["3.00"] = {}
    #d_out["OUT"] = {}

    return d_out


def countGroupDataset() : 

    listS = ListSub()    
    struct = {}
    for element in listS : 
        struct[element] = {}
        struct[element][element] = 0
        struct[element]["PDB"] = 0
        struct[element]["ligand"] = 0
    
    return struct


def countAngleDistance(d_count_angle, distance_max):
    
    l_angle = range(0, 190, 10)
    l_dist = listDistance(distance_max)
    dAngleType = {}
    
    for type_strucutre in d_count_angle.keys():
        dAngleType[type_strucutre] = {}
        for distance in l_dist:       
            dAngleType[type_strucutre][str(distance)] = {}
            angleTemp = 0
            for keyAngle in l_angle : 
                dAngleType[type_strucutre][str(distance)][str(keyAngle)] = {}
            
                for classe in  d_count_angle[type_strucutre].keys():    
                    dAngleType[type_strucutre][str(distance)][str(keyAngle)][classe] = 0
                    nbAngles = len(d_count_angle[type_strucutre][classe]["angles"])
                    for i in range(0, nbAngles) : 
                    # restrint angle position
                        for angleVector in d_count_angle[type_strucutre][classe]["angles"][i] :
                            if angleVector <= keyAngle :
                                if angleVector >= angleTemp :
                                    if d_count_angle[type_strucutre][classe]["distance"][i] <= float(distance) : 
                                        dAngleType[type_strucutre][str(distance)][str(keyAngle)][classe] = dAngleType[type_strucutre][str(distance)][str(keyAngle)][classe] + 1
                angleTemp = keyAngle

    
    for type_strucutre in dAngleType.keys():
        for keyAngle in dAngleType[type_strucutre]["2.0"].keys() : 
            for classe in dAngleType[type_strucutre]["2.0"][keyAngle].keys() : 
                try : dAngleType[type_strucutre]["2.5"][str(keyAngle)][classe] = dAngleType[type_strucutre]["2.5"][str(keyAngle)][classe] - dAngleType[type_strucutre]["2.0"][str(keyAngle)][classe]
                except : pass
                try : dAngleType[type_strucutre]["3.0"][str(keyAngle)][classe] = dAngleType[type_strucutre]["3.0"][str(keyAngle)][classe] - dAngleType[type_strucutre]["2.5"][str(keyAngle)][classe] - dAngleType[type_strucutre]["2.0"][str(keyAngle)][classe]
                except : pass
                try : dAngleType[type_strucutre]["3.5"][str(keyAngle)][classe] = dAngleType[type_strucutre]["3.5"][str(keyAngle)][classe] - dAngleType[type_strucutre]["3.0"][str(keyAngle)][classe] - dAngleType[type_strucutre]["2.5"][str(keyAngle)][classe] - dAngleType[type_strucutre]["2.0"][str(keyAngle)][classe]
                except : pass
                try : dAngleType[type_strucutre]["4.0"][str(keyAngle)][classe] = dAngleType[type_strucutre]["4.0"][str(keyAngle)][classe] - dAngleType[type_strucutre]["3.5"][str(keyAngle)][classe] - dAngleType[type_strucutre]["3.0"][str(keyAngle)][classe] - dAngleType[type_strucutre]["2.5"][str(keyAngle)][classe] - dAngleType[type_strucutre]["2.0"][str(keyAngle)][classe]
                except : pass
                try : dAngleType[type_strucutre]["4.5"][str(keyAngle)][classe] = dAngleType[type_strucutre]["4.5"][str(keyAngle)][classe] - dAngleType[type_strucutre]["4.0"][str(keyAngle)][classe] - dAngleType[type_strucutre]["3.5"][str(keyAngle)][classe] - dAngleType[type_strucutre]["3.0"][str(keyAngle)][classe] - dAngleType[type_strucutre]["2.5"][str(keyAngle)][classe] - dAngleType[type_strucutre]["2.0"][str(keyAngle)][classe]
                except : pass
                try : dAngleType[type_strucutre]["5.0"][str(keyAngle)][classe] = dAngleType[type_strucutre]["5.0"][str(keyAngle)][classe] - dAngleType[type_strucutre]["4.5"][str(keyAngle)][classe] - dAngleType[type_strucutre]["4.0"][str(keyAngle)][classe] - dAngleType[type_strucutre]["3.5"][str(keyAngle)][classe] - dAngleType[type_strucutre]["3.0"][str(keyAngle)][classe] - dAngleType[type_strucutre]["2.5"][str(keyAngle)][classe] - dAngleType[type_strucutre]["2.0"][str(keyAngle)][classe]
                except : pass
                try : dAngleType[type_strucutre]["5.5"][str(keyAngle)][classe] = dAngleType[type_strucutre]["5.5"][str(keyAngle)][classe] - dAngleType[type_strucutre]["5.0"][str(keyAngle)][classe] - dAngleType[type_strucutre]["4.5"][str(keyAngle)][classe] - dAngleType[type_strucutre]["4.0"][str(keyAngle)][classe] - dAngleType[type_strucutre]["3.5"][str(keyAngle)][classe] - dAngleType[type_strucutre]["3.0"][str(keyAngle)][classe] - dAngleType[type_strucutre]["2.5"][str(keyAngle)][classe] - dAngleType[type_strucutre]["2.0"][str(keyAngle)][classe]
                except : pass
                try : dAngleType[type_strucutre]["6.0"][str(keyAngle)][classe] = dAngleType[type_strucutre]["6.0"][str(keyAngle)][classe] - dAngleType[type_strucutre]["5.5"][str(keyAngle)][classe] - dAngleType[type_strucutre]["5.0"][str(keyAngle)][classe] - dAngleType[type_strucutre]["4.5"][str(keyAngle)][classe] - dAngleType[type_strucutre]["4.0"][str(keyAngle)][classe] - dAngleType[type_strucutre]["3.5"][str(keyAngle)][classe] - dAngleType[type_strucutre]["3.0"][str(keyAngle)][classe] - dAngleType[type_strucutre]["2.5"][str(keyAngle)][classe] - dAngleType[type_strucutre]["2.0"][str(keyAngle)][classe]
                except : pass
                
                
    return dAngleType


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


def ListSub ():
    
    return ["I", "II", "III",  "IMD", "GAI", "COO"]
    

def listGlobalStudy():
    return ['residue', 'proportionAtom', 'angleVector', 'ligand', 'ResidueAllAtom', 'distanceOx', 'byAA', 'H2O', 'atom', 'proportionType']
    
def criteraAngle(subs = "", loose = 0):
    """
    Impose angle in different structure
    """
    d_criteria = {}
    
    if loose == 0 : 
        d_criteria["I"] = {}
        d_criteria["I"]["angle"] = [50, 180]
        d_criteria["I"]["distance"] = [2.0,4.0]
        
        d_criteria["II"] = {}
        d_criteria["II"]["angle"] = [50, 180]
        d_criteria["II"]["distance"] = [2.0,4.0]
        
        d_criteria["III"] = {}
        d_criteria["III"]["angle"] = [50, 180]
        d_criteria["III"]["distance"] = [2.0,4.0]
        
        d_criteria["IMD"] = {}
        d_criteria["IMD"]["angle"] = [0, 90]
        d_criteria["IMD"]["distance"] = [2.0,5.0]
        
        # change criterion because I change 11-08 the number of angle considered -> need more data
        d_criteria["GAI"] = {}
        d_criteria["GAI"]["angle"] = [0, 180]#[40, 80, 120, 180]
        d_criteria["GAI"]["distance"] = [2.0,5]#[3.5, 5.5, 3.5, 5.5] #need change distance also !!!!
        
        d_criteria["COO"] = {}
        d_criteria["COO"]["angle"] = [50, 180]
        d_criteria["COO"]["distance"] = [2.0,4.5]
        
        d_criteria["global"] = {}
        d_criteria["global"]["angle"] = [0, 180]
        d_criteria["global"]["distance"] = [2.0,4.0]
        
    else : 
        
        d_criteria["I"] = {}
        d_criteria["I"]["angle"] = [90,130]
        d_criteria["I"]["distance"] = [2.0,4.0]
    
        d_criteria["II"] = {}
        d_criteria["II"]["angle"] = [80,120]
        d_criteria["II"]["distance"] = [2.0,4.0]
    
        d_criteria["III"] = {}
        d_criteria["III"]["angle"] = [90,130]
        d_criteria["III"]["distance"] = [2.0,4.0]
    
        d_criteria["IMD"] = {}
        d_criteria["IMD"]["angle"] = [110,150]
        d_criteria["IMD"]["distance"] = [2.0,5.0]
    
        d_criteria["GAI"] = {}
        d_criteria["GAI"]["angle"] = [0,180]
        d_criteria["GAI"]["distance"] = [2.0,5.0]
   
        d_criteria["COO"] = {}
        d_criteria["COO"]["angle"] = [130,160]
        d_criteria["COO"]["distance"] = [2.0,4.0]
    
    
    if subs == "" : 
        return d_criteria
    else : 
        return d_criteria[subs]



def splitAreaDistance ():

    distStruct = {}
    
    distStruct["I"] = {}
    distStruct["I"]["distance"] = 3.5
    
    distStruct["II"] = {}
    distStruct["II"]["distance"] = 3.5
    
    distStruct["III"] = {}
    distStruct["III"]["distance"] = 3.5
    
    distStruct["IMD"] = {}
    distStruct["IMD"]["distance"] = 3.5
    
    distStruct["GAI"] = {}
    distStruct["GAI"]["distance"] = 4.5
    
    distStruct["COO"] = {}
    distStruct["COO"]["distance"] = 4.25
    
    distStruct["global"] = {}
    distStruct["global"]["distance"] = 3.5
    return distStruct


def substructureCoord (type_substructure):
    
    if type_substructure == "I" : 
        
        atom1 = {}
        atom1["serial"] = 1
        atom1["name"] = "C01"
        atom1["resName"] = "PRI"
        atom1["x"] = 0.000
        atom1["y"] = 0.000
        atom1["z"] = 0.000
        atom1["element"] = "C"
        atom1["charge"] = "0"
        atom1["occupancy"] = "1"
        atom1["tempFactor"] = "0"
        atom1["connect"] = []  
        atom1["char"] = ""
    
        
        atom2 = {}
        atom2["serial"] = 2
        atom2["name"] = "N01"
        atom2["resName"] = "PRI"
        atom2["x"] = 1.118
        atom2["y"] = 0.931
        atom2["z"] = -0.210
        atom2["element"] = "N"
        atom2["charge"] = "0"
        atom2["occupancy"] = "1"
        atom2["tempFactor"] = "0"
        atom2["connect"] = [] 
        atom2["char"] = ""
        return [atom1, atom2]
    
    
    elif type_substructure == "II" : 
        
        atom1 = {}
        atom1["serial"] = 1
        atom1["name"] = "C01"
        atom1["resName"] = "SEC"
        atom1["x"] = 0.000
        atom1["y"] = 0.000
        atom1["z"] = 0.000
        atom1["element"] = "C"
        atom1["charge"] = "0"
        atom1["occupancy"] = "1"
        atom1["tempFactor"] = "0"
        atom1["connect"] = []  
        atom1["char"] = ""
        
        atom2 = {}
        atom2["serial"] = 2
        atom2["name"] = "N01"
        atom2["resName"] = "SEC"
        atom2["x"] = 1.118
        atom2["y"] = 0.931
        atom2["z"] = -0.210
        atom2["element"] = "N"
        atom2["charge"] = "0"
        atom2["occupancy"] = "1"
        atom2["tempFactor"] = "0"
        atom2["connect"] = [] 
        atom2["char"] = ""
        
        atom3 = {}
        atom3["serial"] = 3
        atom3["name"] = "C02"
        atom3["resName"] = "SEC"
        atom3["x"] = 1.469
        atom3["y"] = 1.571
        atom3["z"] = 1.066
        atom3["element"] = "C"
        atom3["charge"] = "0"
        atom3["occupancy"] = "1"
        atom3["tempFactor"] = "0"
        atom3["connect"] = []   
        atom3["char"] = ""     
        
        return [atom2, atom1, atom3]    
    
    
    elif type_substructure == "III" : 
        
        
        atom1 = {}
        atom1["serial"] = 1
        atom1["name"] = "C01"
        atom1["resName"] = "TER"
        atom1["x"] = 0.451
        atom1["y"] = -1.374
        atom1["z"] = -0.266
        atom1["element"] = "C"
        atom1["charge"] = "0"
        atom1["occupancy"] = "1"
        atom1["tempFactor"] = "0"
        atom1["connect"] = []  
        atom1["char"] = ""
        
        atom2 = {}
        atom2["serial"] = 2
        atom2["name"] = "C02"
        atom2["resName"] = "TER"
        atom2["x"] = -0.468
        atom2["y"] = 0.102
        atom2["z"] = 1.390
        atom2["element"] = "C"
        atom2["charge"] = "0"
        atom2["occupancy"] = "1"
        atom2["tempFactor"] = "0"
        atom2["connect"] = [] 
        atom2["char"] = ""
        
        atom3 = {}
        atom3["serial"] = 3
        atom3["name"] = "C03"
        atom3["resName"] = "TER"
        atom3["x"] = -1.099
        atom3["y"] = 0.341
        atom3["z"] = -0.915
        atom3["element"] = "C"
        atom3["charge"] = "0"
        atom3["occupancy"] = "1"
        atom3["tempFactor"] = "0"
        atom3["connect"] = [] 
        atom3["char"] = ""       
        
        atom4 = {}
        atom4["serial"] = 4
        atom4["name"] = "N01"
        atom4["resName"] = "TER"
        atom4["x"] = 0.000
        atom4["y"] = 0.000
        atom4["z"] = 0.000
        atom4["element"] = "N"
        atom4["charge"] = "0"
        atom4["occupancy"] = "1"
        atom4["tempFactor"] = "0"
        atom4["connect"] = [] 
        atom4["char"] = ""
        
        return [atom4, atom1, atom2, atom3]    
    
    elif type_substructure == "IMD" : 
        
        atom1 = {}
        atom1["serial"] = 1
        atom1["name"] = "C2"
        atom1["resName"] = "IMD"
        atom1["x"] = 0.121
        atom1["y"] = -1.265
        atom1["z"] = 0.019
        atom1["element"] = "C"
        atom1["charge"] = "0"
        atom1["occupancy"] = "1"
        atom1["tempFactor"] = "0"
        atom1["connect"] = []  
        atom1["char"] = ""
        
        atom2 = {}
        atom2["serial"] = 2
        atom2["name"] = "N1"
        atom2["resName"] = "IMD"
        atom2["x"] = 1.123
        atom2["y"] = -0.401
        atom2["z"] = 0.073
        atom2["element"] = "N"
        atom2["charge"] = "0"
        atom2["occupancy"] = "1"
        atom2["tempFactor"] = "0"
        atom2["connect"] = [] 
        atom2["char"] = ""
        
        atom3 = {}
        atom3["serial"] = 3
        atom3["name"] = "C4"
        atom3["resName"] = "IMD"
        atom3["x"] = -0.776
        atom3["y"] = 0.744
        atom3["z"] = -0.057
        atom3["element"] = "C"
        atom3["charge"] = "0"
        atom3["occupancy"] = "1"
        atom3["tempFactor"] = "0"
        atom3["connect"] = [] 
        atom3["char"] = ""       
        
        atom4 = {}
        atom4["serial"] = 4
        atom4["name"] = "N3"
        atom4["resName"] = "IMD"
        atom4["x"] = -1.061
        atom4["y"] = -0.592
        atom4["z"] = -0.062
        atom4["element"] = "N"
        atom4["charge"] = "0"
        atom4["occupancy"] = "1"
        atom4["tempFactor"] = "0"
        atom4["connect"] = [] 
        atom4["char"] = ""
        
        atom5 = {}
        atom5["serial"] = 5
        atom5["name"] = "C5"
        atom5["resName"] = "IMD"
        atom5["x"] = 0.569
        atom5["y"] = 0.900
        atom5["z"] = 0.025
        atom5["element"] = "C"
        atom5["charge"] = "0"
        atom5["occupancy"] = "1"
        atom5["tempFactor"] = "0"
        atom5["connect"] = [] 
        atom5["char"] = ""
        
        return [atom2, atom1, atom4, atom3, atom5] 
        
        
    elif type_substructure == "GAI" : 

        atom1 = {}
        atom1["serial"] = 1
        atom1["name"] = "N01"
        atom1["resName"] = "GUA"
        atom1["x"] = -2.547
        atom1["y"] = 0.166
        atom1["z"] = 1.501
        atom1["element"] = "N"
        atom1["charge"] = "0"
        atom1["occupancy"] = "1"
        atom1["tempFactor"] = "0"
        atom1["connect"] = []  
        
        atom2 = {}
        atom2["serial"] = 2
        atom2["name"] = "C01"
        atom2["resName"] = "GUA"
        atom2["x"] = -2.102
        atom2["y"] = 1.011
        atom2["z"] = 0.362
        atom2["element"] = "C"
        atom2["charge"] = "0"
        atom2["occupancy"] = "1"
        atom2["tempFactor"] = "0"
        atom2["connect"] = [] 
        
        atom3 = {}
        atom3["serial"] = 3
        atom3["name"] = "N02"
        atom3["resName"] = "GUA"
        atom3["x"] = -0.810
        atom3["y"] = 1.366
        atom3["z"] = -0.160
        atom3["element"] = "N"
        atom3["charge"] = "0"
        atom3["occupancy"] = "1"
        atom3["tempFactor"] = "0"
        atom3["connect"] = []        
        
        atom4 = {}
        atom4["serial"] = 4
        atom4["name"] = "N03"
        atom4["resName"] = "GUA"
        atom4["x"] = -3.233
        atom4["y"] = 1.358
        atom4["z"] = -0.327
        atom4["element"] = "N"
        atom4["charge"] = "0"
        atom4["occupancy"] = "1"
        atom4["tempFactor"] = "0"
        atom4["connect"] = [] 
        
        atom5 = {}
        atom5["serial"] = 5
        atom5["name"] = "C02"
        atom5["resName"] = "GUA"
        atom5["x"] = -1.845
        atom5["y"] = -0.525
        atom5["z"] = 2.579
        atom5["element"] = "C"
        atom5["charge"] = "0"
        atom5["occupancy"] = "1"
        atom5["tempFactor"] = "0"
        atom5["connect"] = [] 
        
        return [atom5, atom2, atom1, atom3, atom4] 
    
    elif type_substructure == "COO" :
    
        atom1 = {}
        atom1["serial"] = 1
        atom1["name"] = "C01"
        atom1["resName"] = "COO"
        atom1["x"] = 0.000
        atom1["y"] = 0.000
        atom1["z"] = 0.000
        atom1["element"] = "C"
        atom1["charge"] = "0"
        atom1["occupancy"] = "1"
        atom1["tempFactor"] = "0"
        atom1["connect"] = []  
        
        atom2 = {}
        atom2["serial"] = 2
        atom2["name"] = "C02"
        atom2["resName"] = "COO"
        atom2["x"] = -1.356
        atom2["y"] = 0.460
        atom2["z"] = 0.567
        atom2["element"] = "C"
        atom2["charge"] = "0"
        atom2["occupancy"] = "1"
        atom2["tempFactor"] = "0"
        atom2["connect"] = [] 
        
        atom3 = {}
        atom3["serial"] = 3
        atom3["name"] = "O01"
        atom3["resName"] = "COO"
        atom3["x"] = -2.372
        atom3["y"] = -0.009
        atom3["z"] = 0.133
        atom3["element"] = "O"
        atom3["charge"] = "0"
        atom3["occupancy"] = "1"
        atom3["tempFactor"] = "0"
        atom3["connect"] = [] 
        
        atom4 = {}
        atom4["serial"] = 4
        atom4["name"] = "O02"
        atom4["resName"] = "COO"
        atom4["x"] = -1.387
        atom4["y"] = 1.300
        atom4["z"] = 1.457
        atom4["element"] = "O"
        atom4["charge"] = "0"
        atom4["occupancy"] = "1"
        atom4["tempFactor"] = "0"
        atom4["connect"] = [] 
        
       
        return [atom1, atom3, atom2, atom4] 
        
        
def DFile2K (pr_result):
    
    d_file_open = {}
    l_interest = ListSub ()
    l_type_neighbor = classificationATOM (out_list = 1)
    
    for interest in l_interest : 
        d_file_open[interest] = {}
        for type_neighbor in l_type_neighbor : 
            d_file_open[interest][type_neighbor] = open (pr_result + str(interest) + "_" + str (type_neighbor), "w")
    
    
        d_file_open[interest]["density"] = open (pr_result + str(interest) + "_density", "w")
    
    return d_file_open



def closeDFile2K (d_in) :
     
    l_path = []
    for k1 in d_in.keys () : 
        for k2 in d_in[k1].keys() :
            l_path.append (d_in[k1][k2].name) 
            d_in[k1][k2].close ()
    
    return l_path
    
    

        
def nbNeighbor () : 
     
    nbNeighborStruct = {}
    
    nbNeighborStruct["I"] = 4
    nbNeighborStruct["II"] = 3
    nbNeighborStruct["III"] = 3
    nbNeighborStruct["IMD"] = 5
    nbNeighborStruct["GAI"] = 8
    nbNeighborStruct["COO"] = 5
    nbNeighborStruct["global"] = 8

    return nbNeighborStruct
        
        
        
 
 
def genericAtom (x, y, z, serialAtom = 2, nameAtom = "O", char = "", resName = "XXX", ChainID = "A", resSeq = 2, iCode = "", occupancy = 1.0, tempFactor = 10.0, element = "O", charge = 0):
     
    
    out = {}
    out["serial"] = int(serialAtom)
    out["name"] = nameAtom
    out["char"] = char
    out["resName"] = resName
    out["chainID"] = ChainID
    out["resSeq"] = int (resSeq)
    out["iCode"] = str(iCode)
    out["x"] = float (x)
    out["y"] = float (y)
    out["z"] = float (z)
    out["element"] = element
    
    out["charge"] = charge
    out["occupancy"] = occupancy
    out["tempFactor"] = tempFactor
    
    out["connect"] = []
    
    return out
    
           
def Pka (subs = ""):
    
    d_pka = {}
    
    d_pka["I"] = [7.75, 10.63]
    d_pka["II"] = [9.29, 11.01]
    d_pka["III"] = [8.31,10.65]
    d_pka["IMD"] = [5.1, 7.75]
    d_pka["GAI"] = [8.33, 13.71]
    d_pka["COO"] = [1.84, 4.40]
    d_pka["global"] = [0.0, 0.1]
    
    if subs != "" : 
        return d_pka[subs]
    else : 
        return d_pka
    
    
