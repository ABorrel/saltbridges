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
import pathManage

from os import path 
from copy import deepcopy
from numpy import sum
from re import search


def ParseDataSet(p_dataset, one_PDB_by_ligand = 0):
    """Statistic of dataset, repetition ligand in PDB file or in PDB
    out : file with count of repetition in file or files"""

    # log 
    print p_dataset, "path"
    start, logFile = log.initAction("Parsing dataset, ligand representation " + str(path.splitext(path.basename(p_dataset))[0]))
    
    d_dataset = loadFile.resultFilterPDBLigand(p_dataset)
    print d_dataset, "data"

    countAmine = structure.countGroupDataset()
    listCount = []
    listPDB = []
    for element in d_dataset:
        print '**', element
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
        countAmine["Guanidium"] = countAmine["Guanidium"] + (int(flag_guanidium / 3) * count["Number PDB"])

        # divise by 2 because equivalent oxygen
        countAmine["AcidCarboxylic"] = countAmine["AcidCarboxylic"] + (int(flag_acidcarboxylic / 2) * count["Number PDB"])
        listCount.append(count)
        
    numberPDB = len(listPDB)
    if one_PDB_by_ligand == 1 : 
        p_dataset = p_dataset + "OnePDB"
    writeFile.parsingDataSet(listCount, countAmine, numberPDB, p_dataset)
    log.endAction("Parsing dataset, ligand representation", start, logFile)


def distanceAnalysis(stAtm, dir_out, logfile):

    d_file = structure.DFile2K(dir_out)

    # implement structure
    for interestGroup in stAtm.keys():
        for atom in stAtm[interestGroup]:
            if atom["neighbors"] == []:
                continue
            for neighbor in atom["neighbors"]:
                type_atom = structure.classificationATOM(neighbor)
                d_file[interestGroup][type_atom].write (str (neighbor["distance"]) + "\t" + interestGroup + "\t" + atom["resName"] + "\n")
                d_file[interestGroup]["density"].write (str (neighbor["distance"]) + "\t" + type_atom + "\t" + atom["resName"] + "\n")

    l_p_file = structure.closeDFile2K(d_file)
        
    # Run R
    for p_file in l_p_file : 
        if search("density", p_file) : 
            runScriptR.plotDistanceDensity(p_file, logfile)
        else : 
            runScriptR.plotDistance(p_file, logfile)


def angle(struct_neighbor, pr_result, d_max, log_file):

    d_count_global = {}
    log_file.write ("[ANGLE] - count structure implementation\n")

    for substructure in struct_neighbor.keys():
        d_count_global[substructure] = {}
        for central_atom in struct_neighbor[substructure]:
            nbNeighbor = len(central_atom["neighbors"])
            i = 0
            while i < nbNeighbor:
                classif_neighbor = structure.classificationATOM(central_atom["neighbors"][i])
                if not classif_neighbor in d_count_global[substructure].keys () : 
                    d_count_global[substructure][classif_neighbor] = {}
                    d_count_global[substructure][classif_neighbor]["distance"] = [central_atom["neighbors"][i]["distance"]]
                    d_count_global[substructure][classif_neighbor]["angles"] = [central_atom["neighbors"][i]["angle"]]
                else : 
                    d_count_global[substructure][classif_neighbor]["distance"].append(central_atom["neighbors"][i]["distance"])
                    d_count_global[substructure][classif_neighbor]["angles"].append(central_atom["neighbors"][i]["angle"])
                i = i + 1
    
    # write global
    l_p_angle = writeFile.resultAngle(d_count_global, pr_result)
    
    # count for barplot
    d_count_angle = structure.countAngleDistance(d_count_global, d_max)
    
    writeFile.dAngleType(d_count_angle, pr_result)
    
    # plot angle
    runScriptR.plotAngle(l_p_angle, log_file)
    

def atomProx(stAtom, pr_result, max_distance, logFile):
    
    # structure count
    stCount = {}
    
    
    for interestGroup in stAtom.keys(): 
        stCount[interestGroup] = {}
        for atom_central in stAtom[interestGroup]:
            if atom_central["neighbors"] == []:
                continue
            for neighbor in atom_central["neighbors"]: 
                if neighbor["name"] in stCount[interestGroup].keys():
                    stCount[interestGroup][neighbor["name"]] = stCount[interestGroup][neighbor["name"]] + 1 
                else: 
                    stCount[interestGroup][neighbor["name"]] = 1
    
    
    l_files_result = writeFile.resultCount(stCount, "ATM", pr_result)
    for file_result in l_files_result : 
        runScriptR.barplotQuantity(max_distance, "Atoms", file_result, logFile)


def ligandProx(stAtom, pr_result, max_distance, logFile):
    
    # structure count
    stCount = {}

    #water in residue list also, because ligand check only the ligand prox
    l_amino_acid = ["ILE", "LEU", "LYS", "PHE", "TYR", "VAL", "SER", "MET", "ARG", "TRP", "PRO", "GLY", "GLU", "ASN", "HIS", "ALA", "ASP", "GLN", "THR", "CYS", "HOH"]

    for interestGroup in stAtom.keys():
        stCount[interestGroup] = {}
        for atomCentral in stAtom[interestGroup]:
            l_check = []
            for neighbor in atomCentral["neighbors"]:
                if not neighbor["resName"] in l_amino_acid:
                    if not neighbor["resName"] in l_check:
                        if neighbor["resName"] in stCount[interestGroup].keys():
                            stCount[interestGroup][neighbor["resName"]] = stCount[interestGroup][neighbor["resName"]] + 1
                            l_check.append(neighbor["resName"])
                        else:
                            stCount[interestGroup][neighbor["resName"]] = 1
                            l_check.append(neighbor["resName"])
                            
    l_files_result = writeFile.resultCount(stCount, "HET", pr_result)
    for file_result in l_files_result : 
        runScriptR.barplotQuantity(max_distance, "Ligands", file_result, logFile)
    

def resProx(stAtom, pr_result, max_distance, logFile):

    l_atom_mainchain = ["C", "O", "CA", "N", "OXT", "NXT"]
    l_amino_acid = ["ILE", "LEU", "LYS", "PHE", "TYR", "VAL", "SER", "MET", "ARG", "TRP", "PRO", "GLY", "GLU", "ASN", "HIS", "ALA", "ASP", "GLN", "THR", "CYS", "HOH"]

    # structure count
    stCount = {}
    
    

    for interestSubstruct in stAtom.keys():
        stCount[interestSubstruct] = {}
        for central_atom in stAtom[interestSubstruct]:
            if central_atom["neighbors"] == []:
                continue
            for neighbor in central_atom["neighbors"]:
                print neighbor["resSeq"], neighbor["resName"]
                
                if not neighbor["resName"] in l_amino_acid : 
                    continue # case ligand closed
                elif not neighbor["resName"] in stCount[interestSubstruct].keys() :
                    stCount[interestSubstruct][neighbor["resName"]] = {}
                    stCount[interestSubstruct][neighbor["resName"]]["main"] = 0
                    stCount[interestSubstruct][neighbor["resName"]]["side"] = 0
                else : 
                    if neighbor["name"] in l_atom_mainchain:
                        stCount[interestSubstruct][neighbor["resName"]]["main"] = stCount[interestSubstruct][neighbor["resName"]]["main"] + 1
                    else : 
                        stCount[interestSubstruct][neighbor["resName"]]["side"] = stCount[interestSubstruct][neighbor["resName"]]["side"] + 1 
        for aa in l_amino_acid : 
            if not aa in stCount[interestSubstruct].keys () :
                stCount[interestSubstruct][aa] = {}
                stCount[interestSubstruct][aa]["main"] = 0
                stCount[interestSubstruct][aa]["side"] = 0
        
    l_files_result = writeFile.resultCountAA(stCount, pr_result)
    for file_result in l_files_result : 
        runScriptR.barplotQuantity(max_distance, "Residues", file_result, logFile)


def classifResProx(stAtom, pr_result, max_distance, logFile):
    
    # structure count
    stCount = {}
    
    # variable
    l_distance = structure.listDistance(max_distance)
    l_amino_acid = ["ILE", "LEU", "LYS", "PHE", "TYR", "VAL", "SER", "MET", "ARG", "TRP", "PRO", "GLY", "GLU", "ASN", "HIS", "ALA", "ASP", "GLN", "THR", "CYS", "HOH"]
    
    for interestGroup in stAtom.keys(): 
        stCount[interestGroup] = {}
        for atom_central in stAtom[interestGroup]:
            if atom_central["neighbors"] == []:
                continue
            l_check = []
            for neighbor in atom_central["neighbors"]: 
                for distance in l_distance : 
                    if not distance in stCount[interestGroup].keys (): 
                        stCount[interestGroup][distance]={}
                        for aa in l_amino_acid : 
                            stCount[interestGroup][distance][aa]=0
                    if neighbor["distance"]< float(distance) and not neighbor["resSeq"] in l_check: 
                        res = neighbor["resName"]
                        if not res in l_amino_acid : 
                            continue
                        else : 
                            stCount[interestGroup][distance][res] = stCount[interestGroup][distance][res] + 1
                            l_check.append (neighbor["resSeq"])
    
    print stCount.keys ()
    # write file
    l_files_result = writeFile.resultResProx(stCount, max_distance ,pr_result)
    for file_result in l_files_result : 
        runScriptR.barplotResDist(file_result, logFile)
        
    
def atomByAa(stAtom, pr_result ,max_distance, logFile ):


    stCount = {}
    l_amino_acid = ["ILE", "LEU", "LYS", "PHE", "TYR", "VAL", "SER", "MET", "ARG", "TRP", "PRO", "GLY", "GLU", "ASN", "HIS", "ALA", "ASP", "GLN", "THR", "CYS"]
     

    for substruct in stAtom.keys():
        stCount[substruct] = {}
        for atom_central in stAtom[substruct]:
            if atom_central["neighbors"] == []:
                continue
            for neighbor in atom_central["neighbors"]:
                res = neighbor["resName"]
                atom = neighbor["name"]
                if not res in l_amino_acid : 
                    continue
                if not res in stCount[substruct].keys () : 
                    stCount[substruct][res] = {}
                    stCount[substruct][res][">3.5"] = {}
                    stCount[substruct][res]["<3.5"] = {}
                if not atom in stCount[substruct][res]["<3.5"].keys () : 
                    stCount[substruct][res][">3.5"][atom] =  0
                    stCount[substruct][res]["<3.5"][atom] =  0
                    
                if neighbor["distance"] >= 3.5 : 
                    stCount[substruct][res][">3.5"][atom] = stCount[substruct][res][">3.5"][atom] + 1
                else : 
                    print stCount[substruct][res]
                    stCount[substruct][res]["<3.5"][atom] = stCount[substruct][res]["<3.5"][atom] + 1 
                

    l_files_result = writeFile.resultByAA(stCount, max_distance ,pr_result)
    for p_file in l_files_result : 
        runScriptR.barplotQuantityByAA(max_distance, p_file.split ("/")[-1], p_file, logFile)


def numberNeighbor (stAtom, pr_result, max_distance, logFile) : 
    
    
    stCount = {}
    
    # different hist with different thresold
    l_distance = structure.listDistance(max_distance)
    
    for substruct in stAtom.keys () : 
        stCount[substruct] = {}
        for atom_central in stAtom[substruct] :
            if atom_central["neighbors"] == []:
                continue
            
            for distance in l_distance : 
                count = 0
                for neighbor in atom_central["neighbors"] : 
                        if float(neighbor["distance"]) < float(distance) : 
                            count = count + 1
                
                if not distance in stCount[substruct] :
                    stCount[substruct][distance] = []
                
                stCount[substruct][distance].append (count)    
                    
                    
    l_p_filout = writeFile.disributionNumberNeighbor (stCount, pr_result)
    for p_filin in l_p_filout : 
        runScriptR.plotNbNeighbor(p_filin, logFile)


def neighborAtomComposition(stAtom, pr_result, max_distance, logFile) : 
    
    stCount = {}
    for substruct in stAtom.keys () : 
        stCount[substruct] = {}
        searchNeighbor (stAtom, stCount, substruct)
    
    # write files
    l_files_result = writeFile.proportionByPositionNeighbors(stCount, pr_result)
    
    print stCount
    for file_result in l_files_result : 
        runScriptR.proportionAtomClassNeighbor (file_result, logFile)


def firstNeighbor (stAtom, stAtomGlobal, pr_result, logFile):


    stCount = {}
    for substruct in stAtom.keys () : 
        stCount[substruct] = {}
        searchNeighbor (stAtom, stCount, substruct)
    stCount["global"] = {}
    searchNeighbor (stAtomGlobal, stCount, "global")
    
    
    # write files
    l_files_count = writeFile.countFirstNeighbor(stCount, pr_result)
    l_files_dist = writeFile.distanceCountStruct(stCount, pr_result)
    
    
    for file_result in l_files_count : 
        runScriptR.AFCPieFirstNeighbor (file_result, logFile)      
    
    for files_dist in l_files_dist : 
        runScriptR.multiHist(files_dist)  


def allNeighbors (stAtom, stAtomGlobal, pr_result, logFile):

    stCount = {}
    for substruct in stAtom.keys () : 
        stCount[substruct] = {}
        searchNeighbor (stAtom, stCount, substruct)
    
    # global
    stCount["global"] = {}
    searchNeighbor (stAtomGlobal, stCount, "global")
    # write files
    l_files_result = writeFile.countNeighborsAll(stCount, pr_result)
    
    for file_result in l_files_result : 
        runScriptR.AFCPieFirstNeighbor (file_result, logFile)        
    
    
def searchNeighbor (stAtom, stCount, sub_struct):
    """
    Search neigbor in proximity
    """
 
    l_type_atom = structure.classificationATOM("", out_list = 1)
    stCount[sub_struct]["angle1_3"] = []
    stCount[sub_struct]["angle1_2"] = []
    stCount[sub_struct]["angle2_3"] = []
 
 
    neig_temp = deepcopy(stAtom[sub_struct])
    for atom_central in neig_temp : 
        l_neighbor = atom_central["neighbors"]
        nb_neighbor =  len(l_neighbor)   
        if nb_neighbor == 0  : 
            continue
        
                    
        dtemp_angle = {}
        for i in range(1,nb_neighbor+1) :
            classif_first, atom_first = searchMoreClose (l_neighbor) # remove the atom closer
            dtemp_angle[i] = atom_first
                             
            if classif_first == None : 
                continue
            if not i in stCount[sub_struct].keys () :
                stCount[sub_struct][i] = {}
                stCount[sub_struct][i]["distance"] = []
                stCount[sub_struct][i]["classe"] = []
                for type_atom in l_type_atom : 
                    stCount[sub_struct][i][type_atom] = 0
        
            stCount[sub_struct][i]["distance"].append(str(atom_first["distance"]))
            stCount[sub_struct][i]["classe"].append(classif_first)
            stCount[sub_struct][i][classif_first] = stCount[sub_struct][i][classif_first] + 1
        
                 
        # angle
        try : 
            stCount[sub_struct]["angle1_3"].append (calcul.angleVector(dtemp_angle[1], atom_central, dtemp_angle[3]))
        except : 
            stCount[sub_struct]["angle1_3"].append ("NA")
        try :   
            stCount[sub_struct]["angle1_2"].append (calcul.angleVector(dtemp_angle[1], atom_central, dtemp_angle[2]))   
        except : 
            stCount[sub_struct]["angle1_2"].append ("NA") 
        try : 
            stCount[sub_struct]["angle2_3"].append (calcul.angleVector(dtemp_angle[3], atom_central, dtemp_angle[2]))  
        except : 
            stCount[sub_struct]["angle2_3"].append ("NA")
                
    
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
     

def globalRunStatistic(struct_atom_close, global_atom_close, max_distance, pr_result):
    """
    Search close environment of different amines
    arg: -> distance max 
         -> file with dataset
    """
    
    start, logFile = log.initAction("RUN Statistic")

#    # proportion salt bridges
    saltBridges (struct_atom_close, pathManage.resultSaltBridges(pr_result), logFile)

# #    distribution distance interest group and type atoms -> distance type
#     distanceAnalysis(struct_atom_close, repertory.resultDistance(pr_result), logFile)
#       
# #     angle -> directory angles
#     angle(struct_atom_close, pr_result, max_distance, logFile)
#       
# #     global analysis proximity -1 atom ligand // -2 aa type // -3 atom classification
#     ligandProx(struct_atom_close, repertory.countGlobalProx (pr_result, name_in = "hetProx"), max_distance, logFile)
#     atomProx(struct_atom_close, repertory.countGlobalProx (pr_result, name_in = "atmProx"), max_distance, logFile)
#     resProx(struct_atom_close, repertory.countGlobalProx (pr_result, name_in = "resProx"), max_distance, logFile)
#     classifResProx(struct_atom_close, repertory.countGlobalProx (pr_result, name_in = "classifAtmProx"), max_distance, logFile)
#     atomByAa(struct_atom_close, repertory.countGlobalProx (pr_result, name_in = "byAA") ,max_distance, logFile )
#       
#       
# #     analyse number of neighbors -> number of atom type (C, O, N)
#     numberNeighbor (struct_atom_close, repertory.countNeighbor(pr_result, "numberHist"), max_distance, logFile)
#     neighborAtomComposition(struct_atom_close, repertory.countNeighbor(pr_result, "propotionPosition"), max_distance, logFile)
#     firstNeighbor (struct_atom_close, global_atom_close, repertory.countNeighbor(pr_result, "firstNeighbor"), logFile)
#     allNeighbors (struct_atom_close, global_atom_close, repertory.countNeighbor(pr_result, "allNeighbor"), logFile)
    
#   # with two area defintion
#     d_area1, d_area2 = splitTwoArea (struct_atom_close)
#     d_global1, d_global2 = splitTwoArea (global_atom_close)
#     allNeighbors (d_area1, d_global1, repertory.twoArea(pr_result, "neighborArea1"), logFile)
#     allNeighbors (d_area2, d_global2, repertory.twoArea(pr_result, "neighborArea2"), logFile)

#    # combination
    combinationNeighbors (struct_atom_close, pathManage.combination(pr_result), logFile)
    
    
    
    
    
    
# # # # # # # #     atom(struct_atom_close, countStruct[str(max_distance)]["atom"])
# # # # # # # #     ligand(struct_atom_close, countStruct[str(max_distance)]["ligand"])
# # # # # # # #     atomByAa(struct_atom_close, countStruct[str(max_distance)]["byAA"])
# # # # # # # #     
# # # # # # # #     relationNeighbors (struct_atom_close, countStruct[str(max_distance)]["threeAnalysis"])
# # # # # # # #     relationNeighbors (global_atom_close, countStruct[str(max_distance)]["threeAnalysis"])
# # # # # # # # 
# # # # # # # #     # number of neighbor average
# # # # # # # #     numberNeighbor (struct_atom_close, countStruct[str(max_distance)]["numberNeighbors"])
# # # # # # # #     numberNeighbor (global_atom_close, countStruct[str(max_distance)]["numberNeighbors"])
# # # # # # # # 
# # # # # # # #     distance = max_distance
# # # # # # # # 
# # # # # # # #     while distance >= 2:
# # # # # # # # #         print distance
# # # # # # # #         # reduce structure with distance criterion
# # # # # # # #         neighborDistance(distance, max_distance, struct_atom_close)
# # # # # # # #         neighborDistanceList(distance, max_distance, global_atom_close) # analyse every distance
# # # # # # # #         
# # # # # # # # 
# # # # # # # #         proportionAtoms.stAtom(struct_atom_close, countStruct[str(distance)]["proportionAtom"])
# # # # # # # #         proportionType.stAtom(struct_atom_close, countStruct[str(distance)]["proportionType"])
# # # # # # # #         proportionAtoms.globalNeighbors(global_atom_close, countStruct[str(distance)]["proportionAtom"]["Global"])
# # # # # # # #         proportionType.globalNeighbors(global_atom_close, countStruct[str(distance)]["proportionType"]["Global"])
# # # # # # # #         
# # # # # # # #         residue(struct_atom_close, countStruct[str(distance)]["residue"])
# # # # # # # # #                 
# # # # # # # # #         # cumul at least one -> interest group
# # # # # # # #         countAtLeastOne(struct_atom_close, countStruct[str(distance)]["atLeastOne"], ["OxAcid"])
# # # # # # # #         countAtLeastOne(struct_atom_close, countStruct[str(distance)]["atLeastOne"], ["H2O"])
# # # # # # # #         countAtLeastOne(struct_atom_close, countStruct[str(distance)]["atLeastOne"], ["ODonAcc"])
# # # # # # # #         countAtLeastOne(struct_atom_close, countStruct[str(distance)]["atLeastOne"], ["Carom"])
# # # # # # # #         countAtLeastOne(struct_atom_close, countStruct[str(distance)]["atLeastOne"], ["OxAcid", "ODonAcc"])
# # # # # # # #         countAtLeastOne(struct_atom_close, countStruct[str(distance)]["atLeastOne"], ["OxAcid", "ODonAcc", "H2O"])
# # # # # # # #         countAtLeastOne(struct_atom_close, countStruct[str(distance)]["atLeastOne"], ["Nhis", "Nbasic"])
# # # # # # # # #        -> every atom         
# # # # # # # #         countAtLeastOne(global_atom_close, countAtLeastOneGlobal[str(distance)]["atLeastOne"], ["OxAcid"])
# # # # # # # #         countAtLeastOne(global_atom_close, countAtLeastOneGlobal[str(distance)]["atLeastOne"], ["H2O"])
# # # # # # # #         countAtLeastOne(global_atom_close, countAtLeastOneGlobal[str(distance)]["atLeastOne"], ["ODonAcc"])
# # # # # # # #         countAtLeastOne(global_atom_close, countAtLeastOneGlobal[str(distance)]["atLeastOne"], ["Carom"])
# # # # # # # #         countAtLeastOne(global_atom_close, countAtLeastOneGlobal[str(distance)]["atLeastOne"], ["OxAcid", "ODonAcc"])
# # # # # # # #         countAtLeastOne(global_atom_close, countAtLeastOneGlobal[str(distance)]["atLeastOne"], ["OxAcid", "ODonAcc", "H2O"])
# # # # # # # #         countAtLeastOne(global_atom_close, countAtLeastOneGlobal[str(distance)]["atLeastOne"], ["Nhis", "Nbasic"])
# # # # # # # # #                 
# # # # # # # #         globalAtomResidue(global_atom_close, countStruct[str(distance)]["ResidueAllAtom"])
# # # # # # # #         angle(struct_atom_close, countStruct[str(distance)]["angle"])
# # # # # # # #         
# # # # # # # #         
# # # # # # # #         
# # # # # # # #         distance = distance - 0.5
# # # # # # # # #                 
# # # # # # # # # 
# # # # # # # #     writeFile.countGlobalCount(max_distance, countStruct, pr_result)
# # # # # # # #     writeFile.resultAtLeastOne(countAtLeastOneGlobal, countStruct, max_distance, pr_result)

    log.endAction("END Staistic run !!!", start, logFile)







def saltBridges (st_atom, pr_result, logFile):
    
    filout = open (pr_result + "proportionSaltBridges", "w")
    filout_sum = open (pr_result + "proportionSaltBridges.sum", "w")
    st_count = {}
    l_interaction = ["salt-bridges","H-bond","water","other"]
    
    for type_subs in st_atom.keys ():
        
        if not type_subs in st_count.keys ():
            st_count[type_subs] = {}
            st_count[type_subs]["salt-bridges"] = 0
            st_count[type_subs]["H-bond"] = 0
            st_count[type_subs]["water"] = 0
            st_count[type_subs]["other"] = 0
            st_count[type_subs]["out_distance"] = 0
            st_count[type_subs]["out_angle"] = 0
            
        
        for atom_central in st_atom[type_subs] : 
            type_stabilisation = retrieveInteraction (atom_central["neighbors"], type_subs)
            st_count[type_subs][type_stabilisation] = st_count[type_subs][type_stabilisation] + 1
            if type_stabilisation == "out_angle" or type_stabilisation == "out_distance": 
                st_count[type_subs]["other"] = st_count[type_subs]["other"] + 1
    
    filout.write ("\t".join (l_interaction) + "\n")
    for sub in st_count.keys () :
        filout_sum.write ("== " + str (sub) + " ==\n")
        filout.write (sub)
        for interaction in l_interaction : 
            filout.write ("\t" + str(st_count[sub][interaction]))
            filout_sum.write (interaction + ": " + str (st_count[sub][interaction]) + "\n")
        
        filout_sum.write (">-< Out distance: " + str (st_count[sub]["out_distance"]) + "\n")
        filout_sum.write (">-< Out angle: " + str (st_count[sub]["out_angle"]) + "\n")
        filout_sum.write ("Global: " + str (st_count[sub]["salt-bridges"] + st_count[sub]["H-bond"] + st_count[sub]["water"] + st_count[sub]["other"]) + "\n")
        filout.write ("\n")
    
    filout.close ()
    filout_sum.close ()
    # write summary
    
    runScriptR.saltBridgesProportion(pr_result + "proportionSaltBridges")
      
        
def retrieveInteraction (l_atoms, subs) : 
     
    st_angle = structure.criteraAngle()
     
    flag_water = 0
    flag_ox = 0
    flag_hbond = 0
    flag_out_distance = 0
    flag_out_angle = 0
     
    for atom in l_atoms : 
        type_atom = structure.classificationATOM(atom)
        #print atom.keys ()
        
        
        if atom["distance"] >= st_angle[subs]["distance"][0] and atom["distance"] <= st_angle[subs]["distance"][1] : 
            if atom["angle"] != [] and atom["angle"][0] >= st_angle[subs]["angle"][0] and atom["angle"][0] <= st_angle[subs]["angle"][1] : 
                
                print atom["angle"], atom["distance"], "****----***** OK", subs
                
                if subs == "AcidCarboxylic" : 
                    if type_atom == "Nbasic" : 
                        flag_ox = 1
                    elif type_atom == "Ndonnor" or type_atom == "Nhis" : 
                        flag_hbond = 1
                    elif  type_atom == "H2O" : 
                        flag_water = 1
                else : 
                
                    if type_atom == "OxAcid" : 
                        flag_ox = 1
                    elif type_atom == "OxPep" or type_atom == "OxAccept" or type_atom == "ODonAcc" : 
                        flag_hbond = 1
                    elif  type_atom == "H2O" : 
                        flag_water = 1
            else : 
                if subs == "AcidCarboxylic" :
                    if type_atom == "Nbasic" : 
                        flag_out_angle = 1
                else : 
                    if type_atom == "OxAcid" : 
                        flag_out_angle = 1
        else : 
            if subs == "AcidCarboxylic" :
                if type_atom == "Nbasic" : 
                    flag_out_distance = 1
            else : 
                if type_atom == "OxAcid" : 
                    flag_out_distance = 1
                        
                        
                    
#         print atom["angle"], atom["distance"], "****----***** NO-OK"
    
    if flag_ox == 1 : 
        return "salt-bridges"
    if flag_hbond == 1 : 
        return "H-bond"
    if flag_water == 1 : 
        return "water"
    if flag_out_angle == 1 :
        return "out_angle"
    if flag_out_distance == 1 : 
        return "out_distance"
    
    return "other"
     
     
def planarityImidazole (atom_interest_close, p_dir_result) : 
     
    l_imidazole_atom_central = atom_interest_close["Imidazole"]
     
     
    p_dir_result = pathManage.coplorIMD (p_dir_result)
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
     
     
    p_dir_result = pathManage.coplorGUA (p_dir_result)
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
     
     


def combinationNeighbors (st_atom, pr_result, logFile):

    d_nb_neighbor = structure.nbNeighbor ()

    d_count = {}
    
    for struct in st_atom.keys () : 
        if not struct in d_count.keys () : 
            d_count[struct] = {}
            
        for atom_central in st_atom[struct] : 
            nb_ind = d_nb_neighbor[struct]
            l_neighbor = deepcopy(atom_central["neighbors"])
            if len (l_neighbor) < nb_ind : 
                continue
            else : 
                l_combination = []
                for i_neighbors in range (1,nb_ind + 1) : 
                    atom_class, d_atom = searchMoreClose (l_neighbor) 
                    l_combination.append (atom_class)
                
                l_combination.sort ()
                k = "_".join (l_combination)
                if not k in d_count[struct] : 
                    d_count[struct][k] = 0
                
                d_count[struct][k] =  d_count[struct][k]  + 1   
    
    for struct in d_count.keys () : 
        filout = open (pr_result + struct + "_combi", "w")
        for combi in d_count[struct].keys () : 
            filout.write (combi + "\t" + str (d_count[struct][combi]) + "\n")
        
        filout.close ()
        runScriptR.barplotCombination (pr_result + struct + "_combi", logFile)
    
                    
                    
                    
                    
                    
                    
                    
                    
                     
            
            
            
        nb_neighbors = d_nb_neighbor[struct]
        
        

    
    
                
# def neighborDistance(distance, distanceGlobal, struct_neighbor):
# 
#     if distance == distanceGlobal:
#         return
#     else:
#         for substruct in struct_neighbor.keys():
#             for azote in struct_neighbor[substruct]:
#                 nbNeighbor = len(azote["neighbors"])
#                 i = 0
#                 while i < nbNeighbor:
#                     if float(azote["neighbors"][i]["distance"]) >= float(distance):
#                         del azote["neighbors"][i]
#                         nbNeighbor = nbNeighbor - 1
#                     else:
#                         i = i + 1
# 
# 
# def neighborDistanceList(distance, distanceGlobal, listAtom):
# 
#     if distance == distanceGlobal:
#         return
# 
#     else:
#         for atom in listAtom:
#             nbNeighbor = len(atom["neighbors"])
#             i = 0
#             while i < nbNeighbor:
#                 if float(atom["neighbors"][i]["distance"]) >= float(distance):
#                     del atom["neighbors"][i]
#                     nbNeighbor = nbNeighbor - 1
#                 else:
#                     i = i + 1


# def globalAtomResidue (listAtom, count):
# 
#     atomMajorChain = ["C", "O", "CA", "N", "OXT", "NXT"]
# 
#     for atom in listAtom:
#         if atom["neighbors"] == []:
#             continue
#         listCheck = structure.listAminoAcidCheck()
#         for neighbor in atom["neighbors"]:
#             if neighbor["resName"] in count.keys():
#                 if neighbor["name"] in atomMajorChain:
#                     if listCheck[neighbor["resName"]]["main"] == 1:
#                         continue
#                     else:
#                         count[neighbor["resName"]]["main"] = count[neighbor["resName"]]["main"] + 1
#                         listCheck[neighbor["resName"]]["main"] = 1
#                 else:
#                     if listCheck[neighbor["resName"]]["side"] == 1:
#                         continue
#                     else:
#                         count[neighbor["resName"]]["side"] = count[neighbor["resName"]]["side"] + 1
#                         listCheck[neighbor["resName"]]["side"] = 1

  


# def relationNeighbors (struct_neighbor, countStruct) : 
# 
#     l_substruct = structure.listStructure()
#     
#     for substruct in l_substruct : 
#         if substruct == "Primary" : 
#             searchNeighbor (struct_neighbor, countStruct, substruct, 8)
#         elif substruct == "Tertiary" : 
#             searchNeighbor (struct_neighbor, countStruct, substruct, 8)
#         elif substruct == "Imidazole" or substruct == "Secondary": 
#             searchNeighbor (struct_neighbor, countStruct, substruct, 8)
#         else : 
#             searchNeighbor (struct_neighbor, countStruct, substruct,  8)








# def lenBondAnalysis (struct_neighbor, substruct, p_dir_result ):
# 
#     p_dir_result = repertory.bondLength (p_dir_result)
# 
#     l_distance = []
#     l_first = []
#     l_coplar = []
#     l_first_coplar = []
#     l_angle = []
#     
#     for at_central in struct_neighbor[substruct] : 
#         PDB_ID = at_central["PDB"]
#         serial_at_central = at_central["serial"]
#         name_ligand =  at_central["resName"] 
#         
#         # all atom ligand
#         l_at_lig = loadFile.ligandInPDB(PDB_ID, name_ligand)
#         l_at_subs = retrieveAtom.substructure (substruct, serial_at_central, l_at_lig)
#         
#         # first neighbor
#         if len (at_central["neighbors"]) == 0 : 
#             continue
#         else : 
#             classif_first, atom_close = searchMoreClose (at_central["neighbors"], option_lcopy = 1)
#             
#             if substruct == "Tertiary" : 
#                 l_coplar.append (calcul.coplanar(l_at_subs[0], l_at_lig))
#                 l_first_coplar.append (classif_first)
#             
#             # atom from neighbor search
#             N_ref =  l_at_subs[0]
#             
#             l_distance_temp = []
#             # distance
#             for c_atom in l_at_subs[1:] : 
#                 l_distance_temp.append (calcul.distanceTwoatoms(c_atom, N_ref))
#             
#             l_angle.append(sum(atom_close["angle"]))
#             l_distance.append (sum (l_distance_temp))
#             l_first.append (classif_first)
#     
#     
#     writeFile.lenBondType (l_distance, l_first, p_dir_result + substruct + ".dat") 
#     writeFile.lenBondType (l_angle, l_first, p_dir_result + substruct + "_angle.dat") 
#     
#     if substruct == "Tertiary" : 
#         writeFile.lenBondType (l_coplar, l_first_coplar, p_dir_result + substruct + "_coplar.dat") 
#         runScriptR.barplotLenBond (p_dir_result + substruct + "_coplar.dat")
#          
#     runScriptR.barplotLenBond (p_dir_result + substruct + ".dat")
#     runScriptR.barplotLenBond (p_dir_result + substruct + "_angle.dat")        
#             
#             
            

def splitTwoArea (st_atom_sub) : 
    
    st_division = structure.splitAreaDistance()
    
    d_area1 = {}
    d_area2 = {}
    
    for sub in st_atom_sub.keys () :
        
        d_area1[sub] = []
        d_area2[sub] = []
        
        dist_split = st_division[sub]["distance"]
        #print dist_split, "DISTANCE"
        
        for atom_central in st_atom_sub[sub] : 
            if atom_central["neighbors"] == [] : 
                continue
            else :
#                 print 'l892', len (atom_central["neighbors"]) 
                atom_central1 =  deepcopy(atom_central)
                atom_central1["neighbors"] = []
                atom_central2 =  deepcopy(atom_central)
                atom_central2["neighbors"] = []
               
               
            for neighbor in atom_central["neighbors"] : 
                if neighbor["distance"] < dist_split : 
                    atom_central1["neighbors"].append (neighbor)
                else : 
                    atom_central2["neighbors"].append (neighbor)
    
    
            if atom_central1["neighbors"] != [] : 
                d_area1[sub].append (deepcopy(atom_central1))
            if atom_central2["neighbors"] != [] : 
                d_area2[sub].append (deepcopy(atom_central2))
            
            #print "l908", len (atom_central1["neighbors"]), len (atom_central2["neighbors"])
    
    return d_area1, d_area2
            
