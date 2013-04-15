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
    print path_file_dataset, "pqth"
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
        
        for struct in list_struct:
            if struct == "Guanidium" :
                flag_guanidium = flag_guanidium + 1 
                
            elif struct == "Diamine" :
                flag_diamine = flag_diamine + 1 
                
            elif struct == "Pyridine" :
                flag_pyridine = flag_pyridine + 1 
                
            elif struct == "Imidazole" :
                flag_imidazole = flag_imidazole + 1 
            else :
                countAmine[struct] = countAmine[struct] + count["Number PDB"]
            
        countAmine["Imidazole"] = countAmine["Imidazole"] + int(flag_imidazole / 2) * count["Number PDB"]
        countAmine["Pyridine"] = countAmine["Pyridine"] + int(flag_pyridine / 2) * count["Number PDB"]
        countAmine["Diamine"] = countAmine["Diamine"] + int(flag_diamine / 2) * count["Number PDB"] 
        countAmine["Guanidium"] = countAmine["Guanidium"] + (int(flag_guanidium / 2) * count["Number PDB"])
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


def countAtLeastOne (amine, count, type, globalIndex):
   
    
    if globalIndex == 0 :
        count = count[type]
        for typeAmine in amine.keys():
            for nitrogen in amine[typeAmine]:
                flag = 0
                for neighbor in nitrogen["neighbors"]:
                    if neighbor["classificationAtLeastOne"] == type:
                        flag = 1
                        break
                if flag == 0: 
                    count[typeAmine]["others"] = count[typeAmine]["others"] + 1
                else:
                    count[typeAmine][type] = count[typeAmine][type] + 1
                    
    else : 
        count = count[type]
        for atom in amine :
            flag = 0
            for neighbor in atom["neighbors"]:
                if neighbor["classificationAtLeastOne"] == type :
                    flag = 1
                    break
            if flag == 0:
                count["others"] = count["others"] + 1
            else:
                count[type] = count[type] + 1
                
                
                
def countAtLeastOneListStudy (amine, count, ListType, globalIndex):
   
    keyName = ListType[0]
    for i in range(1, len(ListType)) : 
        keyName = keyName + "_" + str(ListType[i])

    if globalIndex == 0 :
        listKey = count.keys()
        if not keyName in listKey : 
            count [keyName] = {}
            listStructureStudy = structure.listStructure()
            for typeAmine in listStructureStudy:
                count[keyName][typeAmine] = {}
                count[keyName][typeAmine][keyName] = 0
                count[keyName][typeAmine]["others"] = 0
        else : 
            for typeAmine in amine.keys():
                for nitrogen in amine[typeAmine]:
                    flag = 0
                    for neighbor in nitrogen["neighbors"]:
                        if neighbor["classificationAtLeastOne"] in ListType:
                            flag = 1
                            break
                    if flag == 0: 
                        count[keyName][typeAmine]["others"] = count[keyName][typeAmine]["others"] + 1
                    else:
                        count[keyName][typeAmine][keyName] = count[keyName][typeAmine][keyName] + 1

    else :
        listKey = count.keys()
        if not keyName in listKey : 
            count[keyName] = {}
            count[keyName][keyName] = 0
            count[keyName]["others"] = 0
        
        for atom in amine :
            flag = 0
            for neighbor in atom["neighbors"]:
                if neighbor["classificationAtLeastOne"] in ListType :
                    flag = 1
                    break
            if flag == 0:
                count[keyName]["others"] = count[keyName]["others"] + 1
            else:
                count[keyName][keyName] = count[keyName][keyName] + 1


def neighborsAmine(distanceMax, path_dataset_file, one_ligandby_complexe, angleOption, dir_result):
    """
    Search close environment of different amines
    arg: -> distance max 
         -> file with dataset
    """
    
    start, logFile = log.initAction("search neighbors in " +path.basename(path_dataset_file))
    list_ligands_in_PDB = loadFile.resultFilterPDBLigand(path_dataset_file)
    nbLigand = len(list_ligands_in_PDB)
    #print list_ligands_in_PDB
    
    # ##Count Structure
    countStruct = structure.countGlobalAmine(distanceMax)  # global structure count
    countAtLeastOneGlobal = structure.countAtLeastOneGlobalStruct(distanceMax)
    
    # ##Write summary file
    filesAmine = writeFile.openFileAmine(dir_result)
    filesWithoutAtLeastOne = writeFile.openFilesWithoutSummary(distanceMax, dir_result)
    
    # pass for check different value
    #i = 1274
    #nbLigand = 1275
    

    
#     i = 0
#     while i < nbLigand : 
#         if list_ligands_in_PDB[i]["name"] == "IMD" : 
#             
#             print i
#             i = nbLigand
#         else :
#             i = i + 1
#         
# 
#     return 

    # inialization    
    i = 139
    while i < nbLigand :
        print "Ligand: " + str(list_ligands_in_PDB[i]["name"]) + " " + str(i) + " " + str(nbLigand)
        logFile.write("Ligand: " + str(list_ligands_in_PDB[i]["name"]) + " " + str(i) + "\n")
        nbPDB = len(list_ligands_in_PDB[i]["PDB"])
        if one_ligandby_complexe == 1 : 
            nbPDB = 1 # take first complexe without selection
            
        j = 0
        while j < nbPDB : 
            print "PDB", list_ligands_in_PDB[i]["PDB"][j], j
            list_atom_ligand = loadFile.ligandInPDB(list_ligands_in_PDB[i]["PDB"][j], list_ligands_in_PDB[i]["name"])
            #print list_atom_ligand
            globalAtom = searchPDB.globalNeighbors(distanceMax, list_atom_ligand, list_ligands_in_PDB[i]["PDB"][j])
            amine = searchPDB.interestGroup(distanceMax, list_atom_ligand, list_ligands_in_PDB[i]["PDB"][j], angleOption)

            distanceAnalysisOxygen(amine, countStruct[str(distanceMax)]["distanceOx"])
            atom(amine, countStruct[str(distanceMax)]["atom"])
            ligand(amine, countStruct[str(distanceMax)]["ligand"])
            atomByAa(amine, countStruct[str(distanceMax)]["byAA"])
            
            writeFile.amine(amine, filesAmine)
            

            distance = distanceMax

            while distance >= 2:
                writeFile.withoutAtLeastOneSummary(amine, filesWithoutAtLeastOne, distance)
                neighborDistance(distance, distanceMax, amine)
                neighborDistanceList(distance, distanceMax, globalAtom)
                proportionAtoms.amine(amine, countStruct[str(distance)]["proportionAtom"])
                proportionType.amine(amine, countStruct[str(distance)]["proportionType"])
                proportionAtoms.globalNeighbors(globalAtom, countStruct[str(distance)]["proportionAtom"]["Global"])
                proportionType.globalNeighbors(globalAtom, countStruct[str(distance)]["proportionType"]["Global"])
                residue(amine, countStruct[str(distance)]["residue"])
                
                
                countAtLeastOne(amine, countStruct[str(distance)]["atLeastOne"], "counterIon", 0)
                countAtLeastOne(amine, countStruct[str(distance)]["atLeastOne"], "H2O", 0)
                countAtLeastOne(amine, countStruct[str(distance)]["atLeastOne"], "amphiprotic", 0)
                countAtLeastOne(amine, countStruct[str(distance)]["atLeastOne"], "Carom", 0)
                countAtLeastOneListStudy(amine, countStruct[str(distance)]["atLeastOne"], ["counterIon", "amphiprotic"], 0)
                countAtLeastOneListStudy(amine, countStruct[str(distance)]["atLeastOne"], ["counterIon", "amphiprotic", "H2O"], 0)
                
                countAtLeastOne(globalAtom, countAtLeastOneGlobal[str(distance)]["atLeastOne"], "counterIon", 1)
                countAtLeastOne(globalAtom, countAtLeastOneGlobal[str(distance)]["atLeastOne"], "H2O", 1)
                countAtLeastOne(globalAtom, countAtLeastOneGlobal[str(distance)]["atLeastOne"], "amphiprotic", 1)
                countAtLeastOne(globalAtom, countAtLeastOneGlobal[str(distance)]["atLeastOne"], "Carom", 1)
                countAtLeastOneListStudy(globalAtom, countAtLeastOneGlobal[str(distance)]["atLeastOne"], ["counterIon", "amphiprotic"], 1)
                countAtLeastOneListStudy(globalAtom, countAtLeastOneGlobal[str(distance)]["atLeastOne"], ["counterIon", "amphiprotic", "H2O"], 1)
                
                globalAtomResidue(globalAtom, countStruct[str(distance)]["ResidueAllAtom"])
                angle(amine, countStruct[str(distance)]["angle"])
                distance = distance - 0.5
                
            j = j + 1
        i = i + 1

    writeFile.closeFileAmine(filesAmine)
    writeFile.closeFilesWithoutSummary(filesWithoutAtLeastOne)
    writeFile.countGlobalAmine(distanceMax, countStruct, dir_result)
    writeFile.resultAtLeastOneGlobal(countAtLeastOneGlobal, distanceMax, dir_result)



    log.endAction("Statistical analysis neighbors " + str(path_dataset_file), start, logFile)



def neighborDistance(distance, distanceGlobal, amine):

    if distance == distanceGlobal:
        return
    else:
        for type in amine.keys():
            for azote in amine[type]:
                nbNeighbor = len(azote["neighbors"])
                i = 0
                while i < nbNeighbor:
                    if azote["neighbors"][i]["distance"] >= distance:
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
                if atom["neighbors"][i]["distance"] >= distance:
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


def angle(amine, countAngle):

    for type in amine.keys():
        for nitrogen in amine[type]:
            nbNeighbor = len(nitrogen["neighbors"])
            i = 0
            while i < nbNeighbor:
                countAngle[type][nitrogen["neighbors"][i]["classification"]]["distance"].append(nitrogen["neighbors"][i]["distance"])
                # for angle in nitrogen["neighbors"][i]["angle"] : 
                countAngle[type][nitrogen["neighbors"][i]["classification"]]["angles"].append(nitrogen["neighbors"][i]["angle"])
                i = i + 1
