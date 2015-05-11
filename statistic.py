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
import superimpose

from os import path 
from copy import deepcopy
from numpy import sum
from re import search


def ParseDataSet(p_dataset, debug = 0):
    """Statistic of dataset, repetition ligand in PDB file or in PDB
    out : file with count of repetition in file or files"""

    # log 
    if debug : print p_dataset, "==path dataset=="
    start, logFile = log.initAction("Parsing dataset, ligand representation " + str(path.splitext(path.basename(p_dataset))[0]))
    
    l_lig_dataset = loadFile.resultFilterPDBLigand(p_dataset)
    if debug : print l_lig_dataset, "dico loaded"

    d_count = structure.countGroupDataset()
    l_count = []
    l_PDB = []
    for lig_dataset in l_lig_dataset:
        if debug : print '**', lig_dataset
        count = structure.countInstanceDataSet()
        # print lig_dataset["name"]
        logFile.write(lig_dataset["name"] + "\n")
        count["name"] = lig_dataset["name"]
        count["Number PDB"] = len(lig_dataset["PDB"])
        
        for PDB in lig_dataset["PDB"] : 
            if not PDB in l_PDB : 
                l_PDB.append(PDB)
        
        l_at_ligand = loadFile.ligandInPDB(lig_dataset["PDB"][0], lig_dataset["name"])
        l_sub_found = searchPDB.interestStructure(l_at_ligand)
        if debug : print l_sub_found, "==l48-statistic=="
        
        # global count -> unique list => no unique because pb with multisub in the ligand
        l_sub_unique = sorted(set(l_sub_found),key=l_sub_found.index) 
        
        for sub_unique in l_sub_unique:
            d_count[sub_unique]["PDB"] = d_count[sub_unique]["PDB"] + count["Number PDB"]
            d_count[sub_unique]["ligand"] = d_count[sub_unique]["ligand"] + 1
        
        
        for sub_found in l_sub_found:
            d_count[sub_found][sub_found] = d_count[sub_found][sub_found] + count["Number PDB"]
        l_count.append(count)
        

    # divise number for complexe queries
    for sub in d_count.keys () : 
        if sub == "Imidazole" or sub == "AcidCarboxylic" : 
            d_count[sub][sub] = d_count[sub][sub] / 2
        elif sub == "Guanidium" : 
            d_count[sub][sub] = d_count[sub][sub] / 3
            
    n_PDB = len(l_PDB)
    writeFile.AnalysisDataSet(l_count, d_count, n_PDB, p_dataset)
    log.endAction("Parsing dataset, ligand representation", start, logFile)


def distanceAnalysis(stAtm, dir_out, logfile):

    d_file = structure.DFile2K(dir_out)

    # implement structure
    for subs in stAtm.keys():
        if subs == "global" : continue
        for atom in stAtm[subs]:
            if atom["neighbors"] == []:
                continue
            for neighbor in atom["neighbors"]:
                type_atom = structure.classificationATOM(neighbor)
                d_file[subs][type_atom].write (str (neighbor["distance"]) + "\t" + subs + "\t" + atom["resName"] + "\n")
                d_file[subs]["density"].write (str (neighbor["distance"]) + "\t" + type_atom + "\t" + atom["resName"] + "\n")

    l_p_file = structure.closeDFile2K(d_file)
        
    # Run R
    for p_file in l_p_file : 
        if search("density", p_file) : 
            runScriptR.plotDistanceDensity(p_file, logfile)
        else : 
            runScriptR.plotDistance(p_file, logfile)


def angleSubs(st_atom, pr_result, d_max, log_file):

    d_count_global = {}
    log_file.write ("[ANGLE] - count structure implementation\n")

    for subs in st_atom.keys():
        if subs == "global" : continue
        d_count_global[subs] = {}
        for central_atom in st_atom[subs]:
            nbNeighbor = len(central_atom["neighbors"])
            i = 0
            while i < nbNeighbor:
                classif_neighbor = structure.classificationATOM(central_atom["neighbors"][i])
                if not classif_neighbor in d_count_global[subs].keys () : 
                    d_count_global[subs][classif_neighbor] = {}
                    d_count_global[subs][classif_neighbor]["distance"] = [central_atom["neighbors"][i]["distance"]]
                    d_count_global[subs][classif_neighbor]["angles"] = [central_atom["neighbors"][i]["angleSubs"]]
                else : 
                    d_count_global[subs][classif_neighbor]["distance"].append(central_atom["neighbors"][i]["distance"])
                    d_count_global[subs][classif_neighbor]["angles"].append(central_atom["neighbors"][i]["angleSubs"])
                i = i + 1
    
    # write global
    l_p_angle = writeFile.resultAngle(d_count_global, pr_result)
    
    # count for barplot
    d_count_angle = structure.countAngleDistance(d_count_global, d_max)
    
    writeFile.dAngleType(d_count_angle, pr_result)
    
    # plot angleSubs
    runScriptR.plotAngle(l_p_angle, log_file)
    

def atomProx(st_atom, pr_result, max_distance, logFile):
    
    # structure count
    stCount = {}
    
    
    for subs in st_atom.keys(): 
        if subs == "global" : continue
        stCount[subs] = {}
        for atom_central in st_atom[subs]:
            if atom_central["neighbors"] == []:
                continue
            for neighbor in atom_central["neighbors"]: 
                if neighbor["name"] in stCount[subs].keys():
                    stCount[subs][neighbor["name"]] = stCount[subs][neighbor["name"]] + 1 
                else: 
                    stCount[subs][neighbor["name"]] = 1
    
    
    l_files_result = writeFile.resultCount(stCount, "ATM", pr_result)
    for file_result in l_files_result : 
        runScriptR.barplotQuantity(max_distance, "Atoms", file_result, logFile)


def ligandProx(st_atom, pr_result, max_distance, logFile):
    
    # structure count
    stCount = {}

    #water in residue list also, because ligand check only the ligand prox
    l_amino_acid = ["ILE", "LEU", "LYS", "PHE", "TYR", "VAL", "SER", "MET", "ARG", "TRP", "PRO", "GLY", "GLU", "ASN", "HIS", "ALA", "ASP", "GLN", "THR", "CYS", "HOH"]

    for subs in st_atom.keys():
        if subs == "global" : continue
        stCount[subs] = {}
        for atomCentral in st_atom[subs]:
            l_check = []
            for neighbor in atomCentral["neighbors"]:
                if not neighbor["resName"] in l_amino_acid:
                    if not neighbor["resName"] in l_check:
                        if neighbor["resName"] in stCount[subs].keys():
                            stCount[subs][neighbor["resName"]] = stCount[subs][neighbor["resName"]] + 1
                            l_check.append(neighbor["resName"])
                        else:
                            stCount[subs][neighbor["resName"]] = 1
                            l_check.append(neighbor["resName"])
                            
    l_files_result = writeFile.resultCount(stCount, "HET", pr_result)
    for file_result in l_files_result : 
        runScriptR.barplotQuantity(max_distance, "Ligands", file_result, logFile)
    

def resProx(st_atom, pr_result, max_distance, logFile):

    l_atom_mainchain = ["C", "O", "CA", "N", "OXT", "NXT"]
    l_amino_acid = ["ILE", "LEU", "LYS", "PHE", "TYR", "VAL", "SER", "MET", "ARG", "TRP", "PRO", "GLY", "GLU", "ASN", "HIS", "ALA", "ASP", "GLN", "THR", "CYS", "HOH"]

    # structure count
    stCount = {}
    
    

    for subs in st_atom.keys():
        if subs == "global" : continue
        stCount[subs] = {}
        for central_atom in st_atom[subs]:
            if central_atom["neighbors"] == []:
                continue
            for neighbor in central_atom["neighbors"]:
                print neighbor["resSeq"], neighbor["resName"]
                
                if not neighbor["resName"] in l_amino_acid : 
                    continue # case ligand closed
                elif not neighbor["resName"] in stCount[subs].keys() :
                    stCount[subs][neighbor["resName"]] = {}
                    stCount[subs][neighbor["resName"]]["main"] = 0
                    stCount[subs][neighbor["resName"]]["side"] = 0
                else : 
                    if neighbor["name"] in l_atom_mainchain:
                        stCount[subs][neighbor["resName"]]["main"] = stCount[subs][neighbor["resName"]]["main"] + 1
                    else : 
                        stCount[subs][neighbor["resName"]]["side"] = stCount[subs][neighbor["resName"]]["side"] + 1 
        for aa in l_amino_acid : 
            if not aa in stCount[subs].keys () :
                stCount[subs][aa] = {}
                stCount[subs][aa]["main"] = 0
                stCount[subs][aa]["side"] = 0
        
    l_files_result = writeFile.resultCountAA(stCount, pr_result)
    for file_result in l_files_result : 
        runScriptR.barplotQuantity(max_distance, "Residues", file_result, logFile)


def classifResProx(st_atom, pr_result, max_distance, logFile):
    
    # structure count
    stCount = {}
    
    # variable
    l_distance = structure.listDistance(max_distance)
    l_amino_acid = ["ILE", "LEU", "LYS", "PHE", "TYR", "VAL", "SER", "MET", "ARG", "TRP", "PRO", "GLY", "GLU", "ASN", "HIS", "ALA", "ASP", "GLN", "THR", "CYS", "HOH"]
    
    for subs in st_atom.keys(): 
        if subs == "global" : continue
        stCount[subs] = {}
        for atom_central in st_atom[subs]:
            if atom_central["neighbors"] == []:
                continue
            l_check = []
            for neighbor in atom_central["neighbors"]: 
                for distance in l_distance : 
                    if not distance in stCount[subs].keys (): 
                        stCount[subs][distance]={}
                        for aa in l_amino_acid : 
                            stCount[subs][distance][aa]=0
                    if neighbor["distance"]< float(distance) and not neighbor["resSeq"] in l_check: 
                        res = neighbor["resName"]
                        if not res in l_amino_acid : 
                            continue
                        else : 
                            stCount[subs][distance][res] = stCount[subs][distance][res] + 1
                            l_check.append (neighbor["resSeq"])
    
    print stCount.keys ()
    # write file
    l_files_result = writeFile.resultResProx(stCount, max_distance ,pr_result)
    for file_result in l_files_result : 
        runScriptR.barplotResDist(file_result, logFile)
        
    
def atomByAa(st_atom, pr_result ,max_distance, logFile ):


    stCount = {}
    l_amino_acid = ["ILE", "LEU", "LYS", "PHE", "TYR", "VAL", "SER", "MET", "ARG", "TRP", "PRO", "GLY", "GLU", "ASN", "HIS", "ALA", "ASP", "GLN", "THR", "CYS"]
     

    for subs in st_atom.keys():
        if subs == "global" : continue
        stCount[subs] = {}
        for atom_central in st_atom[subs]:
            if atom_central["neighbors"] == []:
                continue
            for neighbor in atom_central["neighbors"]:
                res = neighbor["resName"]
                atom = neighbor["name"]
                if not res in l_amino_acid : 
                    continue
                if not res in stCount[subs].keys () : 
                    stCount[subs][res] = {}
                    stCount[subs][res][">3.5"] = {}
                    stCount[subs][res]["<3.5"] = {}
                if not atom in stCount[subs][res]["<3.5"].keys () : 
                    stCount[subs][res][">3.5"][atom] =  0
                    stCount[subs][res]["<3.5"][atom] =  0
                    
                if neighbor["distance"] >= 3.5 : 
                    stCount[subs][res][">3.5"][atom] = stCount[subs][res][">3.5"][atom] + 1
                else : 
                    print stCount[subs][res]
                    stCount[subs][res]["<3.5"][atom] = stCount[subs][res]["<3.5"][atom] + 1 
                

    l_files_result = writeFile.resultByAA(stCount, max_distance ,pr_result)
    for p_file in l_files_result : 
        runScriptR.barplotQuantityByAA(max_distance, p_file.split ("/")[-1], p_file, logFile)


def numberNeighbor (st_atom, pr_result, max_distance, logFile) : 
    
    
    stCount = {}
    
    # different hist with different thresold
    l_distance = structure.listDistance(max_distance)
    
    for subs in st_atom.keys () : 
        if subs == "global" : continue
        stCount[subs] = {}
        for atom_central in st_atom[subs] :
            if atom_central["neighbors"] == []:
                continue
            
            for distance in l_distance : 
                count = 0
                for neighbor in atom_central["neighbors"] : 
                        if float(neighbor["distance"]) < float(distance) : 
                            count = count + 1
                
                if not distance in stCount[subs] :
                    stCount[subs][distance] = []
                
                stCount[subs][distance].append (count)    
                    
                    
    l_p_filout = writeFile.disributionNumberNeighbor (stCount, pr_result)
    for p_filin in l_p_filout : 
        runScriptR.plotNbNeighbor(p_filin, logFile)


def neighborAtomComposition(st_atom, pr_result, max_distance, logFile) : 
    
    stCount = {}
    for subs in st_atom.keys () : 
        if subs == "global" : continue
        stCount[subs] = {}
        searchNeighbor (st_atom, stCount, subs)
    
    # write files
    l_files_result = writeFile.proportionByPositionNeighbors(stCount, pr_result)
    
    for file_result in l_files_result : 
        runScriptR.proportionAtomClassNeighbor (file_result, logFile)


def firstNeighbor (st_atom, pr_result, logFile):


    d_count = {}
    for subs in st_atom.keys () : 
        #if subs == "global" : continue
        d_count[subs] = {}
        searchNeighbor (st_atom, d_count, subs)
    #d_count["global"] = {}
    #searchNeighbor (st_atom["global"], d_count, "global")
    
    
    # write files
    l_files_count = writeFile.countFirstNeighbor(d_count, pr_result)
    l_files_dist = writeFile.distanceCountStruct(d_count, pr_result)
    
    
    for file_result in l_files_count : 
        runScriptR.AFCPieFirstNeighbor (file_result, logFile)      
    
    for files_dist in l_files_dist : 
        runScriptR.multiHist(files_dist)  


def allNeighbors (st_atom, pr_result, logFile):

    st_count = {}
    for subs in st_atom.keys () : 
        st_count[subs] = {}
        searchNeighbor (st_atom, st_count, subs)
    
    # write files
    l_files_result = writeFile.countNeighborsAll(st_count, pr_result)
    
    for file_result in l_files_result : 
        runScriptR.AFCPieFirstNeighbor (file_result, logFile)        
    
    
def searchNeighbor (st_atom, d_stock, subs):
    """
    Search neigbor in proximity
    """
 
    l_type_atom = structure.classificationATOM("", out_list = 1)
    d_stock[subs]["angle1_3"] = []
    d_stock[subs]["angle1_2"] = []
    d_stock[subs]["angle2_3"] = []
    
    d_stock[subs]["distance1_3"] = []
    d_stock[subs]["distance1_2"] = []
    d_stock[subs]["distance2_3"] = []
 
 
    neig_temp = deepcopy(st_atom[subs])
    for atom_central in neig_temp : 
        l_neighbor = atom_central["neighbors"]
        nb_neighbor =  len(l_neighbor)   
        if nb_neighbor == 0  : 
            continue
        
                    
        dtemp_angle = {}
        for i in range(1, nb_neighbor+1) :
            classif_first, atom_first = searchMoreClose (l_neighbor) # remove the atom closer
            dtemp_angle[i] = atom_first
                             
            if classif_first == None : 
                continue
            if not i in d_stock[subs].keys () :
                d_stock[subs][i] = {}
                d_stock[subs][i]["distance"] = []
                d_stock[subs][i]["classe"] = []
                for type_atom in l_type_atom : 
                    d_stock[subs][i][type_atom] = 0
        
            d_stock[subs][i]["distance"].append(str(atom_first["distance"]))
            d_stock[subs][i]["classe"].append(classif_first)
            d_stock[subs][i][classif_first] = d_stock[subs][i][classif_first] + 1
        
                 
        # angleSubs
        try : 
            d_stock[subs]["angle1_3"].append (calcul.angleVector(dtemp_angle[1], atom_central, dtemp_angle[3]))
        except : 
            d_stock[subs]["angle1_3"].append ("NA")
        try :   
            d_stock[subs]["angle1_2"].append (calcul.angleVector(dtemp_angle[1], atom_central, dtemp_angle[2]))   
        except : 
            d_stock[subs]["angle1_2"].append ("NA") 
        try : 
            d_stock[subs]["angle2_3"].append (calcul.angleVector(dtemp_angle[3], atom_central, dtemp_angle[2]))  
        except : 
            d_stock[subs]["angle2_3"].append ("NA")
        
        # distance
        try : 
            d_stock[subs]["distance1_3"].append (calcul.distanceTwoatoms(dtemp_angle[1],  dtemp_angle[3]))
        except : 
            d_stock[subs]["distance1_3"].append ("NA")
        try :   
            d_stock[subs]["distance1_2"].append (calcul.distanceTwoatoms(dtemp_angle[1],  dtemp_angle[2]))   
        except : 
            d_stock[subs]["distance1_2"].append ("NA") 
        try : 
            d_stock[subs]["distance2_3"].append (calcul.distanceTwoatoms(dtemp_angle[2],  dtemp_angle[3])) 
        except : 
            d_stock[subs]["distance2_3"].append ("NA")        
                
    
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
     

def globalRunStatistic(st_atom, max_distance, pr_result):
    """
    Search close environment of different amines
    arg: -> distance max 
         -> file with dataset
    """
    
    
    
    start, logFile = log.initAction("RUN Statistic")
# 
# #    # proportion salt bridges
    saltBridges (st_atom, pathManage.resultSaltBridges(pr_result), logFile)
#  
#     # distribution distance interest group and type atoms -> distance type
    distanceAnalysis(st_atom, pathManage.resultDistance(pr_result), logFile)
#         
#     # angleSubs -> directory angles
    angleSubs(st_atom, pr_result, max_distance, logFile)
#         
#     # global analysis proximity -1 atom ligand // -2 aa type // -3 atom classification
    ligandProx(st_atom, pathManage.countGlobalProx (pr_result, name_in = "hetProx"), max_distance, logFile)
    atomProx(st_atom, pathManage.countGlobalProx (pr_result, name_in = "atmProx"), max_distance, logFile)
    resProx(st_atom, pathManage.countGlobalProx (pr_result, name_in = "resProx"), max_distance, logFile)
    classifResProx(st_atom, pathManage.countGlobalProx (pr_result, name_in = "classifAtmProx"), max_distance, logFile)
    atomByAa(st_atom, pathManage.countGlobalProx (pr_result, name_in = "byAA") ,max_distance, logFile )
#         
#         
#     # analyse number of neighbors -> number of atom type (C, O, N)
    numberNeighbor (st_atom, pathManage.countNeighbor(pr_result, "numberHist"), max_distance, logFile)
    neighborAtomComposition(st_atom, pathManage.countNeighbor(pr_result, "propotionPosition"), max_distance, logFile)
    firstNeighbor (st_atom, pathManage.countNeighbor(pr_result, "firstNeighbor"), logFile)
    allNeighbors (st_atom, pathManage.countNeighbor(pr_result, "allNeighbor"), logFile)
#      
#     # with two area defintion
    d_area1, d_area2 = splitTwoArea (st_atom)
    allNeighbors (d_area1, pathManage.twoArea(pr_result, "neighborArea1"), logFile)
    allNeighbors (d_area2, pathManage.twoArea(pr_result, "neighborArea2"), logFile)
# 
# #    # combination
    combinationNeighbors (st_atom, pathManage.combination(pr_result), logFile)
    combinationNeighborsAngle (st_atom, pathManage.combination(pr_result, "angleSubs"))
    superimpose.SuperimposeFirstNeighbors (st_atom, pathManage.combination(pr_result, "superimposed"))
#     
    
    
    
    
    
    
# # # # # # # #     atom(st_atom, countStruct[str(max_distance)]["atom"])
# # # # # # # #     ligand(st_atom, countStruct[str(max_distance)]["ligand"])
# # # # # # # #     atomByAa(st_atom, countStruct[str(max_distance)]["byAA"])
# # # # # # # #     
# # # # # # # #     relationNeighbors (st_atom, countStruct[str(max_distance)]["threeAnalysis"])
# # # # # # # #     relationNeighbors (st_atom_global, countStruct[str(max_distance)]["threeAnalysis"])
# # # # # # # # 
# # # # # # # #     # number of neighbor average
# # # # # # # #     numberNeighbor (st_atom, countStruct[str(max_distance)]["numberNeighbors"])
# # # # # # # #     numberNeighbor (st_atom_global, countStruct[str(max_distance)]["numberNeighbors"])
# # # # # # # # 
# # # # # # # #     distance = max_distance
# # # # # # # # 
# # # # # # # #     while distance >= 2:
# # # # # # # # #         print distance
# # # # # # # #         # reduce structure with distance criterion
# # # # # # # #         neighborDistance(distance, max_distance, st_atom)
# # # # # # # #         neighborDistanceList(distance, max_distance, st_atom_global) # analyse every distance
# # # # # # # #         
# # # # # # # # 
# # # # # # # #         proportionAtoms.stAtom(st_atom, countStruct[str(distance)]["proportionAtom"])
# # # # # # # #         proportionType.stAtom(st_atom, countStruct[str(distance)]["proportionType"])
# # # # # # # #         proportionAtoms.globalNeighbors(st_atom_global, countStruct[str(distance)]["proportionAtom"]["Global"])
# # # # # # # #         proportionType.globalNeighbors(st_atom_global, countStruct[str(distance)]["proportionType"]["Global"])
# # # # # # # #         
# # # # # # # #         residue(st_atom, countStruct[str(distance)]["residue"])
# # # # # # # # #                 
# # # # # # # # #         # cumul at least one -> interest group
# # # # # # # #         countAtLeastOne(st_atom, countStruct[str(distance)]["atLeastOne"], ["OxAcid"])
# # # # # # # #         countAtLeastOne(st_atom, countStruct[str(distance)]["atLeastOne"], ["H2O"])
# # # # # # # #         countAtLeastOne(st_atom, countStruct[str(distance)]["atLeastOne"], ["ODonAcc"])
# # # # # # # #         countAtLeastOne(st_atom, countStruct[str(distance)]["atLeastOne"], ["Carom"])
# # # # # # # #         countAtLeastOne(st_atom, countStruct[str(distance)]["atLeastOne"], ["OxAcid", "ODonAcc"])
# # # # # # # #         countAtLeastOne(st_atom, countStruct[str(distance)]["atLeastOne"], ["OxAcid", "ODonAcc", "H2O"])
# # # # # # # #         countAtLeastOne(st_atom, countStruct[str(distance)]["atLeastOne"], ["Nhis", "Nbasic"])
# # # # # # # # #        -> every atom         
# # # # # # # #         countAtLeastOne(st_atom_global, countAtLeastOneGlobal[str(distance)]["atLeastOne"], ["OxAcid"])
# # # # # # # #         countAtLeastOne(st_atom_global, countAtLeastOneGlobal[str(distance)]["atLeastOne"], ["H2O"])
# # # # # # # #         countAtLeastOne(st_atom_global, countAtLeastOneGlobal[str(distance)]["atLeastOne"], ["ODonAcc"])
# # # # # # # #         countAtLeastOne(st_atom_global, countAtLeastOneGlobal[str(distance)]["atLeastOne"], ["Carom"])
# # # # # # # #         countAtLeastOne(st_atom_global, countAtLeastOneGlobal[str(distance)]["atLeastOne"], ["OxAcid", "ODonAcc"])
# # # # # # # #         countAtLeastOne(st_atom_global, countAtLeastOneGlobal[str(distance)]["atLeastOne"], ["OxAcid", "ODonAcc", "H2O"])
# # # # # # # #         countAtLeastOne(st_atom_global, countAtLeastOneGlobal[str(distance)]["atLeastOne"], ["Nhis", "Nbasic"])
# # # # # # # # #                 
# # # # # # # #         globalAtomResidue(st_atom_global, countStruct[str(distance)]["ResidueAllAtom"])
# # # # # # # #         angleSubs(st_atom, countStruct[str(distance)]["angleSubs"])
# # # # # # # #         
# # # # # # # #         
# # # # # # # #         
# # # # # # # #         distance = distance - 0.5
# # # # # # # # #                 
# # # # # # # # # 
# # # # # # # #     writeFile.countGlobalCount(max_distance, countStruct, pr_result)
# # # # # # # #     writeFile.resultAtLeastOne(countAtLeastOneGlobal, countStruct, max_distance, pr_result)

    log.endAction("END Staistic run !!!", start, logFile)







def saltBridges (st_atom, pr_result, logFile, debug = 1):
    
    filout = open (pr_result + "proportionSaltBridges", "w")
    filout_sum = open (pr_result + "proportionSaltBridges.sum", "w")
    st_count = {}
    l_interaction = ["salt-bridges","H-bond","water","other"]
    
    for type_subs in st_atom.keys ():
        # remove global
        if type_subs == "global" : 
            continue
        if debug : print "=> control l574 statistic.py", type_subs
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
        filout_sum.write (">-< Out angleSubs: " + str (st_count[sub]["out_angle"]) + "\n")
        filout_sum.write ("Global: " + str (st_count[sub]["salt-bridges"] + st_count[sub]["H-bond"] + st_count[sub]["water"] + st_count[sub]["other"]) + "\n")
        filout.write ("\n")
    
    filout.close ()
    filout_sum.close ()
    # write summary
    
    runScriptR.saltBridgesProportion(pr_result + "proportionSaltBridges")
      
        
def retrieveInteraction (l_atoms, subs, debug = 1) : 
     
    st_angle = structure.criteraAngle(subs)
     
    flag_water = 0
    flag_ox = 0
    flag_hbond = 0
    flag_out_distance = 0
    flag_out_angle = 0
     
    for atom in l_atoms : 
        type_atom = structure.classificationATOM(atom)
        #print atom.keys ()
        # print atom
        for angle in atom["angleSubs"] :
            if debug == 1 : 
                if subs == "Imidazole" : print angle, subs 
            if angle == "NA" : 
                continue
        if atom["distance"] >= st_angle["distance"][0] and atom["distance"] <= st_angle["distance"][1] : 
            if atom["angleSubs"] != [] and atom["angleSubs"][0] >= st_angle["angle"][0] and atom["angleSubs"][0] <= st_angle["angle"][1] : 
                
#                 if debug == 1 : print atom["angleSubs"], atom["distance"], "****----***** OK", subs
                
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
    
    for subs in st_atom.keys () : 
        if not subs in d_count.keys () : 
            d_count[subs] = {}
            
        for atom_central in st_atom[subs] : 
            # number of neighbors considered
            nb_ind = d_nb_neighbor[subs]
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
                if not k in d_count[subs] : 
                    d_count[subs][k] = 0
                
                d_count[subs][k] =  d_count[subs][k]  + 1   
    
    for subs in d_count.keys () : 
        filout = open (pr_result + subs + "_combi", "w")
        for combi in d_count[subs].keys () : 
            filout.write (combi + "\t" + str (d_count[subs][combi]) + "\n")
        
        filout.close ()
        runScriptR.barplotCombination (pr_result + subs + "_combi", logFile)
    
    
    
    
def combinationNeighborsAngle (st_atom, pr_result):

    d_relation_neighbors = {}

    for sub in st_atom.keys () : 
        
        d_relation_neighbors [sub] = {}
        searchNeighbor (st_atom, d_relation_neighbors, sub)
        
    
    l_p_filout = writeFile.RelationAngleDistNeighbors (d_relation_neighbors, pr_result)
    
    for p_filout in l_p_filout : 
        runScriptR.DistVSAngleNeighbor (p_filout)
    
        

    
def lenBondAnalysis (struct_neighbor, substruct, p_dir_result ):
 
    p_dir_result = pathManage.bondLength (p_dir_result)
 
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
             
            l_angle.append(sum(atom_close["angleSubs"]))
            l_distance.append (sum (l_distance_temp))
            l_first.append (classif_first)
     
     
    writeFile.lenBondType (l_distance, l_first, p_dir_result + substruct + ".dat") 
    writeFile.lenBondType (l_angle, l_first, p_dir_result + substruct + "_angle.dat") 
     
    if substruct == "Tertiary" : 
        writeFile.lenBondType (l_coplar, l_first_coplar, p_dir_result + substruct + "_coplar.dat") 
        runScriptR.barplotLenBond (p_dir_result + substruct + "_coplar.dat")
          
    runScriptR.barplotLenBond (p_dir_result + substruct + ".dat")
    runScriptR.barplotLenBond (p_dir_result + substruct + "_angle.dat")        
             
             
            

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
            
