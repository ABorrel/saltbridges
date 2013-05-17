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
from os import path 



def parseDataSet(path_file_dataset):
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
        count["Number PDB"] = len(element["PDB"])
        
        for PDB in element["PDB"] : 
            if not PDB in listPDB : 
                listPDB.append(PDB)
        
        atomLigand = loadFile.ligandInPDB(element["PDB"][0], element["name"])
        list_struct = searchPDB.interestStructure(atomLigand)
        
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
                countAmine[struct] = countAmine[struct] + count["Number PDB"]
            
        countAmine["Imidazole"] = countAmine["Imidazole"] + int(flag_imidazole / 2) * count["Number PDB"]
        countAmine["Pyridine"] = countAmine["Pyridine"] + int(flag_pyridine / 2) * count["Number PDB"]
        countAmine["Diamine"] = countAmine["Diamine"] + int(flag_diamine / 2) * count["Number PDB"] 
        countAmine["Guanidium"] = countAmine["Guanidium"] + (int(flag_guanidium / 2) * count["Number PDB"])
        countAmine["AcidCarboxylic"] = countAmine["AcidCarboxylic"] + (int(flag_acidcarboxylic / 2) * count["Number PDB"])
        listCount.append(count)
        
    numberPDB = len(listPDB)
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
    

    distance = max_distance

    while distance >= 2:
        print distance
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
#        -> every atom         
        countAtLeastOne(global_atom_close, countAtLeastOneGlobal[str(distance)]["atLeastOne"], ["OxAcid"])
        countAtLeastOne(global_atom_close, countAtLeastOneGlobal[str(distance)]["atLeastOne"], ["H2O"])
        countAtLeastOne(global_atom_close, countAtLeastOneGlobal[str(distance)]["atLeastOne"], ["ODonAcc"])
        countAtLeastOne(global_atom_close, countAtLeastOneGlobal[str(distance)]["atLeastOne"], ["Carom"])
        countAtLeastOne(global_atom_close, countAtLeastOneGlobal[str(distance)]["atLeastOne"], ["OxAcid", "ODonAcc"])
        countAtLeastOne(global_atom_close, countAtLeastOneGlobal[str(distance)]["atLeastOne"], ["OxAcid", "ODonAcc", "H2O"])
#                 
        globalAtomResidue(global_atom_close, countStruct[str(distance)]["ResidueAllAtom"])
        angle(atom_interest_close, countStruct[str(distance)]["angle"])
        
        
        
        distance = distance - 0.5
#                 
# 
    writeFile.countGlobalCount(max_distance, countStruct, path_dir_result)
    writeFile.resultAtLeastOne(countAtLeastOneGlobal, countStruct, max_distance, path_dir_result)

#     log.endAction("Statistical analysis neighbors ", start, logFile)



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
