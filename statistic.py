import calcul
import loadFile
import log
import proportionAtoms
import proportionType
import retrieveAtom
import runScriptR
import searchPDB
import structure
import writeFile
import tool
import repertory
import runScriptR

from os import path 
from copy import deepcopy
from numpy import sum


def parseDataSet(path_file_dataset, one_PDB_by_ligand = 0):
    """Statistic of dataset, repetition ligand in PDB file or in PDB
    out : file with count of repetition in file or files"""

    # log 
    print path_file_dataset, "path"
    start, logFile = log.initAction("Parsing dataset, ligand representation " + str(path.splitext(path.basename(path_file_dataset))[0]))
    
    dataSetGlobal = loadFile.resultFilterPDBLigand(path_file_dataset)
    print dataSetGlobal, "data"

    countAmine = structure.countGroupDataset()
    listCount = []
    listPDB = []
    for element in dataSetGlobal:
        count = structure.countInstanceDataSet()
        # print element["name"]
        logFile.write(element["name"] + "\n")
        count["name"] = element["name"]
        if one_PDB_by_ligand == 1 : 
            count["Number PDB"] = 1
            listPDB.append (element["PDB"][0])
        else : 
            count["Number PDB"] = len(element["PDB"])
        
            for PDB in element["PDB"] : 
                if not PDB in listPDB : 
                    listPDB.append(PDB)
        
        l_at_ligand = loadFile.ligandInPDB(element["PDB"][0], element["name"])
        list_struct = searchPDB.interestStructure(l_at_ligand)
        
        flag_guanidium = 0
        flag_diamine = 0
        flag_imidazole = 0
        flag_pyridine = 0
        flag_acidcarboxylic = 0
        for struct in list_struct:
            if struct == "Guanidium" :
                flag_guanidium = flag_guanidium + 1 
            elif struct == "Diamine" :
                flag_diamine = flag_diamine + 1 
            elif struct == "Pyridine" :
                flag_pyridine = flag_pyridine + 1 
            elif struct == "Imidazole" :
                flag_imidazole = flag_imidazole + 1 
            elif struct == "AcidCarboxylic" :
                flag_acidcarboxylic = flag_acidcarboxylic + 1     
                
            else :
                # 
                countAmine[struct] = countAmine[struct] + count["Number PDB"]
            
        countAmine["Imidazole"] = countAmine["Imidazole"] + int(flag_imidazole / 2) * count["Number PDB"]
        countAmine["Pyridine"] = countAmine["Pyridine"] + int(flag_pyridine / 2) * count["Number PDB"]
        countAmine["Diamine"] = countAmine["Diamine"] + int(flag_diamine / 2) * count["Number PDB"] 
        countAmine["Guanidium"] = countAmine["Guanidium"] + (int(flag_guanidium / 3) * count["Number PDB"])
        # divise by 2 because equivalent oxygen
        countAmine["AcidCarboxylic"] = countAmine["AcidCarboxylic"] + (int(flag_acidcarboxylic / 2) * count["Number PDB"])
        listCount.append(count)
        
    numberPDB = len(listPDB)
    if one_PDB_by_ligand == 1 : 
        path_file_dataset = path_file_dataset + "OnePDB"
    writeFile.parsingDataSet(listCount, countAmine, numberPDB, path_file_dataset)
    log.endAction("Parsing dataset, ligand representation", start, logFile)


def distanceAnalysisOxygen(amine, count):


    for interestGroup in amine.keys():
        for atom in amine[interestGroup]:
            if atom["neighbors"] == []:
                continue
            for neighbor in atom["neighbors"]:
                countOxygen(neighbor, count[interestGroup])



def countOxygen(neighbor, count):

    if neighbor["name"] == "OD1" or neighbor["name"] == "OD2" or neighbor["name"] == "OE1" or neighbor["name"] == "OE2":
        count.append(neighbor["distance"])


def atom(amine, count):
    
    for typeAmine in amine.keys(): 
        for atomAzote in amine[typeAmine]:
            if atomAzote["neighbors"] == []:
                continue
            for neighbor in atomAzote["neighbors"]: 
                if neighbor["name"] in count[typeAmine].keys():
                    count[typeAmine][neighbor["name"]] = count[typeAmine][neighbor["name"]] + 1 
                else: 
                    count[typeAmine][neighbor["name"]] = 1



def residue(amine, count):

    atomMajorChain = ["C", "O", "CA", "N", "OXT", "NXT"]

    for typeAmine in amine.keys():
        for atomAzote in amine[typeAmine]:
            if atomAzote["neighbors"] == []:
                continue
            listCheck = structure.listAminoAcidCheck()
            for neighbor in atomAzote["neighbors"]:
                if neighbor["resName"] in count[typeAmine].keys():
                    if neighbor["name"] in atomMajorChain:
                        if listCheck[neighbor["resName"]]["main"] == 1:
                            continue
                        else:
                            count[typeAmine][neighbor["resName"]]["main"] = count[typeAmine][neighbor["resName"]]["main"] + 1
                            listCheck[neighbor["resName"]]["main"] = 1
                    else:
                        if listCheck[neighbor["resName"]]["side"] == 1:
                            continue
                        else:
                            count[typeAmine][neighbor["resName"]]["side"] = count[typeAmine][neighbor["resName"]]["side"] + 1
                            listCheck[neighbor["resName"]]["side"] = 1



def ligand(amine, count):

    listAminoAcid = ["ILE", "LEU", "LYS", "PHE", "TYR", "VAL", "SER", "MET", "ARG", "TRP", "PRO", "GLY", "GLU", "ASN", "HIS", "ALA", "ASP", "GLN", "THR", "CYS"]
    listCheck = []

    for typeAmine in amine.keys():
        for atomAzote in amine[typeAmine]:
            for neighbor in atomAzote["neighbors"]:
                if not neighbor["resName"] in listAminoAcid:
                    if not neighbor["resSeq"] in listCheck:
                        if neighbor["resName"] in count[typeAmine].keys():
                            count[typeAmine][neighbor["resName"]] = count[typeAmine][neighbor["resName"]] + 1
                            listCheck.append(neighbor["resSeq"])
                        else:
                            count[typeAmine][neighbor["resName"]] = 1
                            listCheck.append(neighbor["resSeq"])


def atomByAa(amine, count):

    for typeAmine in amine.keys():
        if typeAmine != "global":
            for atomAzote in amine[typeAmine]:
                if atomAzote["neighbors"] == []:
                    continue
                for neighbor in atomAzote["neighbors"]:
                    if neighbor["resName"] in count[typeAmine].keys():
                        if not neighbor["name"] in count[typeAmine][neighbor["resName"]].keys():
                            count[typeAmine][neighbor["resName"]][neighbor["name"]] = {}
                            count[typeAmine][neighbor["resName"]][neighbor["name"]]["3.5"] = 0
                            count[typeAmine][neighbor["resName"]][neighbor["name"]]["4.5"] = 0
                        if not neighbor["name"] in count["global"][neighbor["resName"]].keys():
                            count["global"][neighbor["resName"]][neighbor["name"]] = {}
                            count["global"][neighbor["resName"]][neighbor["name"]]["3.5"] = 0
                            count["global"][neighbor["resName"]][neighbor["name"]]["4.5"] = 0

                        if neighbor["distance"] < 3.5:
                            count[typeAmine][neighbor["resName"]][neighbor["name"]]["3.5"] = count[typeAmine][neighbor["resName"]][neighbor["name"]]["3.5"] + 1
                            count["global"][neighbor["resName"]][neighbor["name"]]["3.5"] = count["global"][neighbor["resName"]][neighbor["name"]]["3.5"] + 1
                        else:
                            count[typeAmine][neighbor["resName"]][neighbor["name"]]["4.5"] = count[typeAmine][neighbor["resName"]][neighbor["name"]]["4.5"] + 1
                            count["global"][neighbor["resName"]][neighbor["name"]]["4.5"] = count["global"][neighbor["resName"]][neighbor["name"]]["4.5"] + 1


def countAtLeastOne (struct_neighbor, struct_count, l_classifSearch):
    
    
    # check key structure count
    k_struct_count = "_".join(l_classifSearch)
    if not k_struct_count in struct_count.keys () : 
        struct_count[k_struct_count] = {}

        # in struct
        struct_count = struct_count[k_struct_count]
    else :
        struct_count = struct_count[k_struct_count]
    
    
    if type (struct_neighbor) is list : 
        struct_count["other"] = 0.0001 # -> proportion calcul
        struct_count[k_struct_count] = 0.0 # -> proportion calcul
        implementAtLeastOne(struct_neighbor, struct_count, k_struct_count,l_classifSearch)
    else :        
        for type_subsearch in struct_neighbor.keys():
            struct_count[type_subsearch] = {}
            struct_count[type_subsearch]["other"] = 0.0001 # -> proportion calcul
            struct_count[type_subsearch][k_struct_count] = 0.0 # -> proportion calcul
            implementAtLeastOne(struct_neighbor[type_subsearch], struct_count[type_subsearch], k_struct_count, l_classifSearch,)


def implementAtLeastOne (l_central_atom, struct_count_in, key_atleastOne_type, l_type_atLeastOne):
    
    for atom_central in l_central_atom:
        struct_count_in["other"] = struct_count_in["other"] + 1
        for neighbor in atom_central["neighbors"]:
            classif_neighbor = structure.classificationATOM(neighbor)
            if classif_neighbor in l_type_atLeastOne:
                struct_count_in[key_atleastOne_type] = struct_count_in[key_atleastOne_type] + 1
                struct_count_in["other"] = struct_count_in["other"] - 1
                break # found at least one

                    
def planarityImidazole (atom_interest_close, p_dir_result) : 
    
    l_imidazole_atom_central = atom_interest_close["Imidazole"]
    
    
    p_dir_result = repertory.coplorIMD (p_dir_result)
    p_filout = p_dir_result + "coplarRing.txt"
    filout = open (p_filout, "w")

    nb_imd = len (l_imidazole_atom_central)
    i = 0
    while i < nb_imd : 
        PDB_ID = l_imidazole_atom_central[i]["PDB"]
        serial_at_central = l_imidazole_atom_central[i]["serial"]
        name_ligand =  l_imidazole_atom_central[i]["resName"]
        
        # load ligand
        l_at_lig = loadFile.ligandInPDB(PDB_ID, name_ligand)
        
        # load structure
        l_at_subs = retrieveAtom.substructure ("Imidazole", serial_at_central, l_at_lig)
        
        # coplar
        d_coplar = calcul.coplanarPoint(l_at_subs[3], [l_at_subs[0],l_at_subs[1], l_at_subs[2]])
        filout.write (str(d_coplar) + "\n")
        
        # criterion dell no IMD -> not used
#         if d_coplar > 0.01 :
#             del l_imidazole_atom_central[i]
#             nb_imd = nb_imd - 1
#         else : 
        i = i + 1
    
        filout.write (str(d_coplar) + "\n")
    filout.close ()
    
    
def planarityGuanidium (atom_interest_close, p_dir_result) : 
    
    l_guanidium_atom_central = atom_interest_close["Guanidium"]
    
    
    p_dir_result = repertory.coplorGUA (p_dir_result)
    p_filout = p_dir_result + "coplarRing.txt"
    filout = open (p_filout, "w")

    nb_gua = len (l_guanidium_atom_central)
    i = 0
    while i < nb_gua : 
        PDB_ID = l_guanidium_atom_central[i]["PDB"]
        serial_at_central = l_guanidium_atom_central[i]["serial"]
        name_ligand =  l_guanidium_atom_central[i]["resName"]
        
        # load ligand
        l_at_lig = loadFile.ligandInPDB(PDB_ID, name_ligand)
        
        # load structure
        l_at_subs = retrieveAtom.substructure ("Guanidium", serial_at_central, l_at_lig)
#         print len (l_at_subs)
#         print l_at_subs[0]["element"], l_at_subs[1]["element"], l_at_subs[2]["element"], l_at_subs[3]["element"], l_at_subs[4]["element"]
        
        # coplar
        d_coplar = calcul.coplanarPoint(l_at_subs[0], [l_at_subs[1],l_at_subs[2], l_at_subs[3]])
#         print d_coplar
        filout.write (str(d_coplar) + "\n")
        
        # criterion dell no IMD -> not used
        if d_coplar > 0.5 :
            del l_guanidium_atom_central[i]
            nb_gua = nb_gua - 1
        else : 
            i = i + 1
    
        
        filout.write (str(d_coplar) + "\n")
    filout.close ()   
    
    


def globalRunStatistic(atom_interest_close, global_atom_close, max_distance, option_angle, path_dir_result):
    """
    Search close environment of different amines
    arg: -> distance max 
         -> file with dataset
    """
    
#     start, logFile = log.initAction("satistic")

    
    # ##Count Structure
    countStruct = structure.countGlobalAmine(max_distance)  # global structure count
    countAtLeastOneGlobal = structure.countAtLeastOneGlobalStruct(max_distance)
    
    distanceAnalysisOxygen(atom_interest_close, countStruct[str(max_distance)]["distanceOx"])
    atom(atom_interest_close, countStruct[str(max_distance)]["atom"])
    ligand(atom_interest_close, countStruct[str(max_distance)]["ligand"])
    atomByAa(atom_interest_close, countStruct[str(max_distance)]["byAA"])
    
    relationNeighbors (atom_interest_close, countStruct[str(max_distance)]["threeAnalysis"])
    relationNeighbors (global_atom_close, countStruct[str(max_distance)]["threeAnalysis"])

    # number of neighbor average
    numberNeighbor (atom_interest_close, countStruct[str(max_distance)]["numberNeighbors"])
    numberNeighbor (global_atom_close, countStruct[str(max_distance)]["numberNeighbors"])

    distance = max_distance

    while distance >= 2:
#         print distance
        # reduce structure with distance criterion
        neighborDistance(distance, max_distance, atom_interest_close)
        neighborDistanceList(distance, max_distance, global_atom_close) # analyse every distance
        

        proportionAtoms.amine(atom_interest_close, countStruct[str(distance)]["proportionAtom"])
        proportionType.amine(atom_interest_close, countStruct[str(distance)]["proportionType"])
        proportionAtoms.globalNeighbors(global_atom_close, countStruct[str(distance)]["proportionAtom"]["Global"])
        proportionType.globalNeighbors(global_atom_close, countStruct[str(distance)]["proportionType"]["Global"])
        
        residue(atom_interest_close, countStruct[str(distance)]["residue"])
#                 
#         # cumul at least one -> interest group
        countAtLeastOne(atom_interest_close, countStruct[str(distance)]["atLeastOne"], ["OxAcid"])
        countAtLeastOne(atom_interest_close, countStruct[str(distance)]["atLeastOne"], ["H2O"])
        countAtLeastOne(atom_interest_close, countStruct[str(distance)]["atLeastOne"], ["ODonAcc"])
        countAtLeastOne(atom_interest_close, countStruct[str(distance)]["atLeastOne"], ["Carom"])
        countAtLeastOne(atom_interest_close, countStruct[str(distance)]["atLeastOne"], ["OxAcid", "ODonAcc"])
        countAtLeastOne(atom_interest_close, countStruct[str(distance)]["atLeastOne"], ["OxAcid", "ODonAcc", "H2O"])
        countAtLeastOne(atom_interest_close, countStruct[str(distance)]["atLeastOne"], ["Nhis", "Nbasic"])
#        -> every atom         
        countAtLeastOne(global_atom_close, countAtLeastOneGlobal[str(distance)]["atLeastOne"], ["OxAcid"])
        countAtLeastOne(global_atom_close, countAtLeastOneGlobal[str(distance)]["atLeastOne"], ["H2O"])
        countAtLeastOne(global_atom_close, countAtLeastOneGlobal[str(distance)]["atLeastOne"], ["ODonAcc"])
        countAtLeastOne(global_atom_close, countAtLeastOneGlobal[str(distance)]["atLeastOne"], ["Carom"])
        countAtLeastOne(global_atom_close, countAtLeastOneGlobal[str(distance)]["atLeastOne"], ["OxAcid", "ODonAcc"])
        countAtLeastOne(global_atom_close, countAtLeastOneGlobal[str(distance)]["atLeastOne"], ["OxAcid", "ODonAcc", "H2O"])
        countAtLeastOne(global_atom_close, countAtLeastOneGlobal[str(distance)]["atLeastOne"], ["Nhis", "Nbasic"])
#                 
        globalAtomResidue(global_atom_close, countStruct[str(distance)]["ResidueAllAtom"])
        angle(atom_interest_close, countStruct[str(distance)]["angle"])
        
        
        
        distance = distance - 0.5
#                 
# 
    writeFile.countGlobalCount(max_distance, countStruct, path_dir_result)
    writeFile.resultAtLeastOne(countAtLeastOneGlobal, countStruct, max_distance, path_dir_result)

#     log.endAction("Statistical analysis neighbors ", start, logFile)



def numberNeighbor (struct_neighbor, count) : 
    
    if type (struct_neighbor) is list : 
        if not "Global" in count.keys () : 
            count["Global"] = []
        for atom_central in struct_neighbor : 
            count["Global"].append (len(atom_central["neighbors"]))
            
    else :
        for substruct in struct_neighbor.keys():
            for atom_central in struct_neighbor[substruct]:
                nbNeighbor = len(atom_central["neighbors"])
                if not substruct in count.keys () : 
                    count[substruct] = []
                
                count[substruct].append (nbNeighbor)




def neighborDistance(distance, distanceGlobal, struct_neighbor):

    if distance == distanceGlobal:
        return
    else:
        for substruct in struct_neighbor.keys():
            for azote in struct_neighbor[substruct]:
                nbNeighbor = len(azote["neighbors"])
                i = 0
                while i < nbNeighbor:
                    if float(azote["neighbors"][i]["distance"]) >= float(distance):
                        del azote["neighbors"][i]
                        nbNeighbor = nbNeighbor - 1
                    else:
                        i = i + 1


def neighborDistanceList(distance, distanceGlobal, listAtom):

    if distance == distanceGlobal:
        return

    else:
        for atom in listAtom:
            nbNeighbor = len(atom["neighbors"])
            i = 0
            while i < nbNeighbor:
                if float(atom["neighbors"][i]["distance"]) >= float(distance):
                    del atom["neighbors"][i]
                    nbNeighbor = nbNeighbor - 1
                else:
                    i = i + 1


def globalAtomResidue (listAtom, count):

    atomMajorChain = ["C", "O", "CA", "N", "OXT", "NXT"]

    for atom in listAtom:
        if atom["neighbors"] == []:
            continue
        listCheck = structure.listAminoAcidCheck()
        for neighbor in atom["neighbors"]:
            if neighbor["resName"] in count.keys():
                if neighbor["name"] in atomMajorChain:
                    if listCheck[neighbor["resName"]]["main"] == 1:
                        continue
                    else:
                        count[neighbor["resName"]]["main"] = count[neighbor["resName"]]["main"] + 1
                        listCheck[neighbor["resName"]]["main"] = 1
                else:
                    if listCheck[neighbor["resName"]]["side"] == 1:
                        continue
                    else:
                        count[neighbor["resName"]]["side"] = count[neighbor["resName"]]["side"] + 1
                        listCheck[neighbor["resName"]]["side"] = 1


def angle(struct_neighbor, countAngle):

    for substructure in struct_neighbor.keys():
        for nitrogen in struct_neighbor[substructure]:
            nbNeighbor = len(nitrogen["neighbors"])
            i = 0
            while i < nbNeighbor:
                classif_neighbor = structure.classificationATOM(nitrogen["neighbors"][i])
                countAngle[substructure][classif_neighbor]["distance"].append(nitrogen["neighbors"][i]["distance"])
                # for angle in nitrogen["neighbors"][i]["angle"] : 
                countAngle[substructure][classif_neighbor]["angles"].append(nitrogen["neighbors"][i]["angle"])
                i = i + 1


def relationNeighbors (struct_neighbor, countStruct) : 

    l_substruct = structure.listStructure()
    
    for substruct in l_substruct : 
        if substruct == "Primary" : 
            searchNeighbor (struct_neighbor, countStruct, substruct, 8)
        elif substruct == "Tertiary" : 
            searchNeighbor (struct_neighbor, countStruct, substruct, 8)
        elif substruct == "Imidazole" or substruct == "Secondary": 
            searchNeighbor (struct_neighbor, countStruct, substruct, 8)
        else : 
            searchNeighbor (struct_neighbor, countStruct, substruct,  8)





def searchNeighbor (struct_neighbor, countStruct, sub_struct, nb_neighbor):
    """
    Search neigbor in proximity
    """

    if type (struct_neighbor) is list : 
        sub_struct = "global"
        neig_temp = deepcopy(struct_neighbor)
        for atom_central in neig_temp : 
            l_neighbor = atom_central["neighbors"]
            if len(l_neighbor) == 0 : # no neighbor continue 
                continue
            
            dtemp_angle = {}
            # count type of classes atoms
            for i in range(1,nb_neighbor+1) : 
                classif_first, atom_close = searchMoreClose (l_neighbor)
                dtemp_angle[i] = atom_close
                if classif_first == None : 
                    continue
                countStruct[sub_struct][i]["distance"].append(str(atom_close["distance"]))
                countStruct[sub_struct][i]["classe"].append(classif_first)
                countStruct[sub_struct][i][classif_first] = countStruct[sub_struct][i][classif_first] + 1
            
            # angles
            try : 
                countStruct[sub_struct]["angle1_3"].append (calcul.angleVector(dtemp_angle[1], atom_central, dtemp_angle[3]))
            except : 
                countStruct[sub_struct]["angle1_3"].append ("NA")
                
            try :   
                countStruct[sub_struct]["angle1_2"].append (calcul.angleVector(dtemp_angle[1], atom_central, dtemp_angle[2]))   
            except : 
                countStruct[sub_struct]["angle1_2"].append ("NA") 
                
            try : 
                countStruct[sub_struct]["angle2_3"].append (calcul.angleVector(dtemp_angle[3], atom_central, dtemp_angle[2]))  
            except : 
                countStruct[sub_struct]["angle2_3"].append ("NA")  
                
        
    else :
                    neig_temp = deepcopy(struct_neighbor[sub_struct])
                    for atom_central in neig_temp : 
                        l_neighbor = atom_central["neighbors"]
                        if len(l_neighbor) == 0  : 
                            continue
                        
                        dtemp_angle = {}
                        for i in range(1,nb_neighbor+1) : 
                            classif_first, atom_close = searchMoreClose (l_neighbor)
                            dtemp_angle[i] = atom_close
                            
                            if classif_first == None : 
                                continue
                            countStruct[sub_struct][i]["distance"].append(str(atom_close["distance"]))
                            countStruct[sub_struct][i]["classe"].append(classif_first)
                            countStruct[sub_struct][i][classif_first] = countStruct[sub_struct][i][classif_first] + 1
                
                        # angle
                        try : 
                            countStruct[sub_struct]["angle1_3"].append (calcul.angleVector(dtemp_angle[1], atom_central, dtemp_angle[3]))
                        except : 
                            countStruct[sub_struct]["angle1_3"].append ("NA")
                        try :   
                            countStruct[sub_struct]["angle1_2"].append (calcul.angleVector(dtemp_angle[1], atom_central, dtemp_angle[2]))   
                        except : 
                            countStruct[sub_struct]["angle1_2"].append ("NA") 
                        try : 
                            countStruct[sub_struct]["angle2_3"].append (calcul.angleVector(dtemp_angle[3], atom_central, dtemp_angle[2]))  
                        except : 
                            countStruct[sub_struct]["angle2_3"].append ("NA")
            
        
    
def searchMoreClose (l_neighbors, option_lcopy = 0) : 
    
    if l_neighbors == [] : 
        return None, None
    
    if option_lcopy == 1 : 
        l_neighbors_use = deepcopy(l_neighbors)
    else : 
        l_neighbors_use = l_neighbors
    
    # no neighbor
    if len (l_neighbors_use) == 0 : 
        return None, None
    
    d = 10
    i_out = 0
    nb_neigbor = len (l_neighbors_use)
    i = 0
    while i < nb_neigbor : 
        if l_neighbors_use[i]["distance"] < d : 
            classe_out = structure.classificationATOM(l_neighbors_use[i])
            d = l_neighbors_use[i]["distance"]
            atom_out = deepcopy(l_neighbors_use[i])
            i_out = i
        i = i + 1 
    
    del l_neighbors_use[i_out]
    return classe_out, atom_out


def lenBondAnalysis (struct_neighbor, substruct, p_dir_result ):

    p_dir_result = repertory.bondLength (p_dir_result)

    l_distance = []
    l_first = []
    l_coplar = []
    l_first_coplar = []
    l_angle = []
    
    for at_central in struct_neighbor[substruct] : 
        PDB_ID = at_central["PDB"]
        serial_at_central = at_central["serial"]
        name_ligand =  at_central["resName"] 
        
        # all atom ligand
        l_at_lig = loadFile.ligandInPDB(PDB_ID, name_ligand)
        l_at_subs = retrieveAtom.substructure (substruct, serial_at_central, l_at_lig)
        
        # first neighbor
        if len (at_central["neighbors"]) == 0 : 
            continue
        else : 
            classif_first, atom_close = searchMoreClose (at_central["neighbors"], option_lcopy = 1)
            
            if substruct == "Tertiary" : 
                l_coplar.append (calcul.coplanar(l_at_subs[0], l_at_lig))
                l_first_coplar.append (classif_first)
            
            # atom from neighbor search
            N_ref =  l_at_subs[0]
            
            l_distance_temp = []
            # distance
            for c_atom in l_at_subs[1:] : 
                l_distance_temp.append (calcul.distanceTwoatoms(c_atom, N_ref))
            
            l_angle.append(sum(atom_close["angle"]))
            l_distance.append (sum (l_distance_temp))
            l_first.append (classif_first)
    
    
    writeFile.lenBondType (l_distance, l_first, p_dir_result + substruct + ".dat") 
    writeFile.lenBondType (l_angle, l_first, p_dir_result + substruct + "_angle.dat") 
    
    if substruct == "Tertiary" : 
        writeFile.lenBondType (l_coplar, l_first_coplar, p_dir_result + substruct + "_coplar.dat") 
        runScriptR.barplotLenBond (p_dir_result + substruct + "_coplar.dat")
         
    runScriptR.barplotLenBond (p_dir_result + substruct + ".dat")
    runScriptR.barplotLenBond (p_dir_result + substruct + "_angle.dat")        
            
            
            
        
    
    
    
    
    
    
    



