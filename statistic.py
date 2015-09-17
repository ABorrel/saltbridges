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
import structure

from os import path 
from copy import deepcopy
from numpy import sum, mean, std
from re import search


def ParseDataSet(p_dataset, debug = 0):
    """Statistic of dataset, repetition ligand in PDB file or in PDB
    out : file with d_count of repetition in file or files"""

    # log 
    if debug : print p_dataset, "==path dataset=="
    
    l_d_lig = loadFile.resultFilterPDBLigand(p_dataset)
    if debug : print l_d_lig, "dico loaded"

    d_sub = structure.countGroupDataset()
    l_count = []
    l_PDB_global = []
    for d_lig in l_d_lig:
        if debug : print '**', d_lig
        d_count = {}
        # print d_lig["name"]
        d_count["name"] = d_lig["name"]
        d_count["Number PDB"] = len(d_lig["PDB"])
        
        for PDB in d_lig["PDB"] : 
            if not PDB in l_PDB_global : 
                l_PDB_global.append (PDB)
        
        l_at_ligand = loadFile.ExtractInfoPDBID(PDB)[d_lig["name"]][0]
        
        
        print d_lig["name"]
        
        l_sub_found = searchPDB.interestStructure(l_at_ligand)
        print l_sub_found
        
        
# # #         if l_sub_found == [] : 
# # #             print d_lig["PDB"], d_lig["name"]
# # #             gggg
        
        if debug : print l_sub_found, "==l48-statistic=="
        
        # global d_count -> unique list => no unique because pb with multisub in the ligand
        l_sub_unique = sorted(set(l_sub_found),key=l_sub_found.index) 
        
        for sub_unique in l_sub_unique:
            d_sub[sub_unique]["PDB"] = d_sub[sub_unique]["PDB"] + d_count["Number PDB"]
            d_sub[sub_unique]["ligand"] = d_sub[sub_unique]["ligand"] + 1
        
        for sub_found in l_sub_found:
            
            if sub_found == "IMD" : print l_sub_found, d_lig["name"], d_lig["PDB"][0]
            
            d_sub[sub_found][sub_found] = d_sub[sub_found][sub_found] + d_count["Number PDB"]
            
        l_count.append(d_count)
        

    # divise number for complexe queries
    for sub in d_sub.keys () : 
        if sub == "IMD" or sub == "COO" : 
            print d_sub[sub][sub]
            d_sub[sub][sub] = d_sub[sub][sub] / 2
        elif sub == "GAI" : 
            d_sub[sub][sub] = d_sub[sub][sub] / 3
            
    n_PDB = len(l_PDB_global)
    writeFile.AnalysisDataSet(l_count, d_sub, n_PDB, p_dataset)
    print d_sub["IMD"]



def globalRunStatistic(st_atom, max_distance, pr_result):
    """
    Search close environment of different amines
    arg: -> distance max 
         -> file with dataset
    """
    
    print st_atom.keys ()
    
    start, logFile = log.initAction("RUN Statistic")
# # 
    # proportion salt bridges
    CountInteraction (st_atom, pathManage.resultInteraction(pr_result), logFile)
    EnvironmentInteraction (st_atom, pathManage.resultInteraction (pr_result, name_in = "conditional"), logFile)

    # loose    
#     CountInteraction (st_atom, pathManage.resultInteraction(pr_result, "loose"), logFile, restrained = 0)
#     EnvironmentInteraction (st_atom, pathManage.resultInteraction (pr_result, name_in = "loose/conditional"), logFile, restrained = 0)    
    
# #  
# #     # distribution distance interest group and type atoms -> distance type
    DistTypeAtom(st_atom, pathManage.resultDistance(pr_result), logFile)
# #         
# #     # angleSubs -> directory angles
    AngleSubs(st_atom, pr_result, max_distance)
    AngleSelect (st_atom, pr_result)
# #         
# #     # global analysis proximity -1 atom ligand // -2 aa type // -3 atom classification
#     ligandProx(st_atom, pathManage.countGlobalProx (pr_result, name_in = "hetProx"), max_distance, logFile)
#     atomProx(st_atom, pathManage.countGlobalProx (pr_result, name_in = "atmProx"), max_distance, logFile)
#     resProx(st_atom, pathManage.countGlobalProx (pr_result, name_in = "resProx"), max_distance, logFile)
#     classifResProx(st_atom, pathManage.countGlobalProx (pr_result, name_in = "classifAtmProx"), max_distance, logFile)
#     atomByAa(st_atom, pathManage.countGlobalProx (pr_result, name_in = "byAA") ,max_distance, logFile )
# #         
# #         
# #     # analyse number of neighbors -> number of atom type (C, O, N)
    MeansNumberNeighbors (st_atom, pathManage.countNeighbor(pr_result, "MeansNumberNeighbor"))
    ResidueClose (st_atom, pathManage.countNeighbor(pr_result, "residuesNeighbor"), logFile)
#     numberNeighbor (st_atom, pathManage.countNeighbor(pr_result, "numberHist"), max_distance, logFile)
#     neighborAtomComposition(st_atom, pathManage.countNeighbor(pr_result, "propotionPosition"), max_distance, logFile)
#     firstNeighbor (st_atom, pathManage.countNeighbor(pr_result, "firstNeighbor"), logFile)
#     allNeighbors (st_atom, pathManage.countNeighbor(pr_result, "allNeighbor"), logFile)
    
# #      
# #     # with two area defintion
#     d_area1, d_area2 = splitTwoArea (st_atom)
#     allNeighbors (d_area1, pathManage.twoArea(pr_result, "neighborArea1"), logFile)
#     allNeighbors (d_area2, pathManage.twoArea(pr_result, "neighborArea2"), logFile)
# # 
# # #    # combination
#     combinationNeighbors (st_atom, pathManage.combination(pr_result), logFile)
#     combinationNeighborsAngle (st_atom, pathManage.combination(pr_result, "angleSubs"))
#     superimpose.SuperimposeFirstNeighbors (st_atom, pathManage.combination(pr_result, "superimposed"))
# #     

    log.endAction("END Statistic run !!!", start, logFile)


def MergeDataSet (pr_result, name_dataset1, name_dataset2):
    
    # merge proportion interaction
    pr_merged_portion = pathManage.CreatePathDir(pr_result + "MergeInteractProportion/")
    
    # path file dataset 1
    p_atleastone1 = pathManage.resultInteraction(pr_result + name_dataset1 + "/") + "interact_dependant"
    p_atleastonecum1 = pathManage.resultInteraction(pr_result + name_dataset1 + "/") + "interact_independant"
    
    # path file dataset 2
    p_atleastone2 = pathManage.resultInteraction(pr_result + name_dataset2 + "/") + "interact_dependant"
    p_atleastonecum2 = pathManage.resultInteraction(pr_result + name_dataset2 + "/") + "interact_independant"
    

    if not path.exists(p_atleastone1) or not path.exists(p_atleastonecum1) or not path.exists(p_atleastonecum2) or not path.exists(p_atleastone2) : 
        print p_atleastone1
        print p_atleastone2
        print p_atleastonecum1
        print p_atleastonecum2
        print pr_merged_portion
        print "ERROR-file, l165 statistic"
        return 

    else : 
        runScriptR.MergeProportionAndDataset(p_atleastone1, p_atleastonecum1, p_atleastone2, p_atleastonecum2, pr_merged_portion)

def DistTypeAtom(stAtm, dir_out, logfile):

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


def AngleSubs(st_atom, pr_result, d_max):

    d_count_global = {}

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
                    if subs == "IMD" :
                        angle =  central_atom["neighbors"][i]["angleSubs"][0]
                        if angle >= 90.0 : 
                            angle = 180.0 - angle
                        d_count_global[subs][classif_neighbor]["angles"].append([angle])
                    else : 
                        d_count_global[subs][classif_neighbor]["angles"].append(central_atom["neighbors"][i]["angleSubs"])
                i = i + 1
    
    # write global
    l_p_angle = writeFile.resultAngle(d_count_global, pr_result)
    
    # count for barplot
    d_count_angle = structure.countAngleDistance(d_count_global, d_max)
    
    writeFile.dAngleType(d_count_angle, pr_result)
    
    # plot angleSubs
    runScriptR.plotAngle(l_p_angle)
    

def AngleSelect (st_atom, pr_result) : 
    
    criteria = structure.criteraAngle()
    
    for sub in st_atom.keys () : 
        # do not need global
        if sub == "global" : 
            continue
        l_distlimit = criteria[sub]["distance"]
        nb_test = len (l_distlimit) / 2
        d_file = {}
        d_file[sub] = {}
        i = 0
        while i < nb_test : 
            min_dist = l_distlimit[i]
            max_dist = l_distlimit[i + 1]
            name_folder = str (sub) + str (min_dist) + "-" + str (max_dist)
            
            for atom_central in st_atom[sub] : 
                for neighbor in atom_central["neighbors"] : 
                    if neighbor ["distance"] >= min_dist and neighbor["distance"] <= max_dist : 
                        type_atom = structure.classificationATOM(neighbor)
                        if not type_atom in d_file[sub].keys () : 
                            p_filetype = pathManage.ResultAngleCriteria(pr_result, name_folder) + "angle_" + str (type_atom)
                            d_file[sub][type_atom] = open (p_filetype, "w")
                        for angle in neighbor["angleSubs"] : 
                            if sub == "IMD" and angle >= 90 : 
                                angle = 180 - angle
                            d_file[sub][type_atom].write(str (angle) + "\t" + str(neighbor["distance"]) + "\t" + str (type_atom) + "\t" + str(atom_central["resName"]) + "\t" + str (atom_central["PDB"]) + "\n")
            i = i + 2
        
        for k in d_file[sub].keys () : 
            d_file[sub][k].close ()
            
            # run plot
            runScriptR.HistAngle (d_file[sub][k].name)
    



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




def MeansNumberNeighbors (st_atom, pr_result) : 
    
    criteria = structure.criteraAngle()
    l_atom_type = structure.classificationATOM(out_list = 1)
    
    d_list = {}
    
    for subs in st_atom.keys () :
        if not subs in d_list.keys (): 
            d_list[subs] = {} 
        for atom_sub in st_atom[subs] : 
            if not "Nb" in d_list[subs].keys () : 
                d_list[subs]["Nb"] = []
                
            d_count = {} # temp count by subs
            d_count["Nb"] = 0
            for type_atom in l_atom_type : 
                if not type_atom in d_list[subs].keys () : 
                    d_list[subs][type_atom] = []
                if not type_atom in d_count : 
                    d_count[type_atom] = 0
            
            
            for atom_neighbor in atom_sub["neighbors"] :
                if atom_neighbor["distance"] >= criteria[subs]["distance"][0] and atom_neighbor["distance"] <= criteria[subs]["distance"][1] : 
                    type_atom = structure.classificationATOM(atom_neighbor)
                    d_count[type_atom] = d_count[type_atom] + 1
                    d_count["Nb"] = d_count["Nb"] + 1
                      
            for type_atom in l_atom_type : 
                d_list[subs][type_atom].append (d_count[type_atom])
            d_list[subs]["Nb"].append (d_count["Nb"])
    
    
    
    p_filout = pr_result + "meanNumberofNeighbors"
    filout = open (p_filout, "w")

    # header
    filout.write ("M\tSD")
    for type_atom in l_atom_type : 
        filout.write ("\t" + str (type_atom) + "\t" + str (type_atom) + "SD")
    filout.write ("\n")
    
    # result
    l_sub = structure.ListSub()
    l_sub.append ("global")
    for subs in l_sub : 
        filout.write (str(subs) + "\t" + str (mean (d_list[subs]["Nb"])) + "\t" + str (std (d_list[subs]["Nb"])))
        for type_atom in l_atom_type : 
            filout.write ("\t" + str (mean (d_list[subs][type_atom])) + "\t" + str (std (d_list[subs][type_atom])))
        filout.write ("\n")
    filout.close ()
        
    runScriptR.MeansNumberNeighbors (p_filout)
    
    
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
        runScriptR.AFCBarplot (file_result, logFile)      
    
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
        runScriptR.AFCBarplot (file_result, logFile)        
    


def ResidueClose (st_atom, pr_result, logFile) : 
    
    d_count = {}
    d_list = {}
    
    for subs in st_atom.keys () : 
        d_count[subs] = {}
        d_list[subs] = {}
        SearchResiduesComposition (st_atom[subs], d_count[subs], d_list[subs], subs)
    
    ################
    # global Count #
    ################
    
    #writeFile
    d_file = writeFile.CountNeighborRes (d_count, pr_result)
    #run
    l_p_count = d_file["count"]
    for p_count in l_p_count : 
        runScriptR.plotNbNeighbor(p_count, logFile)
    
    l_p_res = d_file["res"]
    for p_res in l_p_res : 
        runScriptR.AFCBarplot (p_res, logFile)
        

    ################
    # Means and SD #
    ################
    
    p_filout = pr_result + "meanNumberResNeighbors"
    filout = open (p_filout, "w")

    l_res = structure.l_res
    # header
    filout.write ("M\tSD")
    for res in l_res : 
        filout.write ("\t" + str (res) + "\t" + str (res) + "SD")
    filout.write ("\n")
    
    # result
    l_sub = structure.ListSub()
    l_sub.append ("global")
    for subs in l_sub : 
        filout.write (str(subs) + "\t" + str (mean (d_count[subs]["count"])) + "\t" + str (std (d_count[subs]["count"])))
        for res in l_res : 
            filout.write ("\t" + str (mean (d_list[subs][res])) + "\t" + str (std (d_list[subs][res])))
        filout.write ("\n")
    filout.close ()
    
    runScriptR.MeansNumberNeighbors (p_filout)
    
    
def SearchResiduesComposition (st_atom, d_count, d_list, subs):
    """
    Count the proportion in AA
    -> maybe need optimization with variable distance cut off
    """
    criteria = structure.criteraAngle()
    
    l_res = structure.l_res

    # initialize structure stock
    d_count["res"] = {}
    d_count["count"] = []
    for res in l_res : 
        d_count["res"][res] = 0
        d_list[res] = []
    
    for atom_subs in st_atom :
        l_res_temp = [] 
        d_temp = {}
        for res in l_res : 
            d_temp[res] = 0
        for neighbor in atom_subs["neighbors"] : 
            if neighbor ["distance"] >= criteria[subs]["distance"][0] and neighbor ["distance"] <= criteria[subs]["distance"][1] : 
                if not neighbor["resName"] in l_res : 
                    continue
                res_temp = str(neighbor["resName"] + "_" + str (neighbor["resSeq"]))
    #             print res_temp
                if not res_temp in l_res_temp : 
                    l_res_temp.append(res_temp)
                    d_temp[neighbor["resName"]] =  d_temp[neighbor["resName"]] + 1
        
        d_count["count"].append (len (l_res_temp))
        
        for res in l_res : 
            d_list[res].append (d_temp[res])
            d_count["res"][res] = d_count["res"][res] + d_temp[res]
        
    
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
     




def CountInteraction (st_atom, pr_result, logFile, restrained = 1, debug = 1):
    
    st_count = {}
    l_interactions_search = ["COO", "H2O", "OH", "NH", "N", "Other"]
    
    for type_subs in st_atom.keys ():
        # remove global
        if type_subs == "global" : 
            continue
        if debug : print "=> control l574 statistic.py", type_subs
        if not type_subs in st_count.keys ():
            st_count[type_subs] = {}
            st_count[type_subs]["dependant"] = {}
            st_count[type_subs]["independant"] = {}
            #initialisation
            for interact in l_interactions_search : 
                st_count[type_subs]["independant"][interact] = 0
                st_count[type_subs]["dependant"][interact] = 0
            
            print st_count
            
            
        for atom_central in st_atom[type_subs] : 
            l_interact_found = GetInteractions (atom_central["neighbors"], type_subs, restrained)
            
            CountIndependant (l_interact_found, st_count[type_subs]["independant"])
            CountDependant(l_interact_found, st_count[type_subs]["dependant"], type_subs)
            

    p_filout_dependant = pr_result + "interact_dependant"
    p_filout_independant = pr_result + "interact_independant"
    
    filout_dependant = open (p_filout_dependant, "w")
    filout_independant = open (p_filout_independant, "w")
    
    # header
    filout_dependant.write ("\t".join (l_interactions_search) + "\n")
    filout_independant.write ("\t".join (l_interactions_search[:-1]) + "\n") 
    
    for sub in st_count.keys () :
        filout_dependant.write (sub + "\t")
        filout_independant.write (sub + "\t")
        
        filout_dependant.write ("\t".join([str (st_count[sub]["dependant"][k]) for k in l_interactions_search]) + "\n")
        filout_independant.write ("\t".join([str (st_count[sub]["independant"][k]) for k in l_interactions_search[:-1]]) + "\n")                       
                                

    filout_dependant.close ()
    filout_independant.close ()
    
    
    runScriptR.InteractionProportion(p_filout_dependant)
    runScriptR.InteractionProportion(p_filout_independant)

    # merge plot 
    runScriptR.MergeProportionInteractAtLeasNotAtLeast (p_filout_dependant, p_filout_independant)
    


def CountIndependant (l_interaction_found, d_count):
    
    
    for interact in l_interaction_found : 
        d_count[interact] = d_count[interact] + 1

def CountDependant (l_interaction_found, d_count, type_subs):
    
    print l_interaction_found, type_subs
    
    # for group COO
    if type_subs == "COO" : 
        if "N" in l_interaction_found : 
            d_count["N"] = d_count["N"] + 1
        elif "NH" in l_interaction_found : 
            d_count["NH"] = d_count["NH"] + 1
        elif "H2O" in l_interaction_found : 
            d_count["H2O"] = d_count["H2O"] + 1
        elif "OH" in l_interaction_found : 
            d_count["OH"] = d_count["OH"] + 1
        elif "COO" in l_interaction_found : 
            d_count["COO"] = d_count["COO"] + 1
        else : 
            d_count["Other"] = d_count["Other"] + 1
    else : 
        if "COO" in l_interaction_found : 
            d_count["COO"] = d_count["COO"] + 1
        elif "H2O" in l_interaction_found : 
            d_count["H2O"] = d_count["H2O"] + 1    
        elif "OH" in l_interaction_found : 
            d_count["OH"] = d_count["OH"] + 1        
        elif "NH" in l_interaction_found : 
            d_count["NH"] = d_count["NH"] + 1
        elif "N" in l_interaction_found : 
            d_count["N"] = d_count["N"] + 1
        else : 
            d_count["Other"] = d_count["Other"] + 1       
        
        
    
      
        
def GetInteractions (l_atoms, subs, restrained = 1, debug = 1) : 
     
    if restrained == 1 : 
        st_angle = structure.criteraAngle(subs)
    else : 
        st_angle = structure.criteraAngle(subs, loose = 1)
    
    l_interact = []
    
     
    for atom in l_atoms : 
        flag_temp_dist = 0
        flag_temp_angle = 0
        type_atom = structure.classificationATOM(atom)
        #print atom.keys ()
        # print atom    
        l_angle = atom["angleSubs"]
        dist = atom["distance"]
        # chech distance
        if not dist >= st_angle["distance"][0] or not dist <= st_angle["distance"][1] :
            flag_temp_dist = 1
        # check angle
        else : 
            for angle in l_angle : 
                if not angle >= st_angle["angle"][0] or not angle <= st_angle["angle"][1] :
                    flag_temp_angle = 1
                    
            if l_angle == [] : 
                flag_temp_angle = 1 
        
        # in the criteria    
        if flag_temp_dist == 0 and flag_temp_angle == 0  : # angle and distance OK
            if type_atom == "Carg" or type_atom == "Nhis":
                if not "N" in l_interact : 
                    l_interact.append ("N")
            elif type_atom == "COxAcid" : 
                if not "COO" in l_interact : 
                    l_interact.append ("COO")
            elif type_atom == "H2O" : 
                if not "H2O" in l_interact : 
                    l_interact.append ("H2O")
            elif type_atom == "ODonAcc" : 
                if not "OH" in l_interact : 
                    l_interact.append ("OH")
            elif type_atom == "Ndonor" : 
                if not "NH" in l_interact : 
                    l_interact.append ("NH")
            else : 
                l_interact.append ("Other")
    
    return l_interact
     
     
def planarityImidazole (atom_interest_close, p_dir_result) : 
     
    l_imidazole_atom_central = atom_interest_close["IMD"]
     
     
    p_dir_result = pathManage.coplorIMD (p_dir_result)
    p_filout = p_dir_result + "coplarRing.txt"
    filout = open (p_filout, "w")
 
    nb_imd = len (l_imidazole_atom_central)
    i = 0
    while i < nb_imd : 
        PDB_ID = l_imidazole_atom_central[i]["PDB"]
        name_ligand =  l_imidazole_atom_central[i]["resName"]
         
        # load ligand
        l_at_lig = loadFile.ExtractInfoPDBID(PDB_ID)[name_ligand][0] # change not tested
         
        # load structure
        l_at_subs = retrieveAtom.substructure ("IMD", l_imidazole_atom_central[i], l_at_lig)
         
        # coplar
        try : d_coplar = calcul.coplanarPoint(l_at_subs[3], [l_at_subs[0],l_at_subs[1], l_at_subs[2]])
        except : 
            i = i + 1
            continue
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
     
    l_guanidium_atom_central = atom_interest_close["GAI"]
     
     
    p_dir_result = pathManage.coplorGUA (p_dir_result)
    p_filout = p_dir_result + "coplarRing.txt"
    filout = open (p_filout, "w")
 
    nb_gua = len (l_guanidium_atom_central)
    i = 0
    while i < nb_gua : 
        PDB_ID = l_guanidium_atom_central[i]["PDB"]
        name_ligand =  l_guanidium_atom_central[i]["resName"]
         
        # load ligand
        l_at_lig = loadFile.ExtractInfoPDBID(PDB_ID)[name_ligand][0] # change not tested
         
        # load structure
        l_at_subs = retrieveAtom.substructure ("GAI", l_guanidium_atom_central[i], l_at_lig)
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
     
     

def ListNeighborsType (l_neighbors, subs):
    
    l_out = []
    d_nb_neighbor = structure.nbNeighbor ()
    nb_neighbor = d_nb_neighbor[subs]
    l_temp = deepcopy(l_neighbors)

    for i_neighbors in range (1, nb_neighbor + 1) : 
        if l_temp == [] : 
            l_out.append ("empty")
        else : 
            type_close, d_atom_close = searchMoreClose (l_temp)
            l_out.append (type_close)
            
    return l_out


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
        name_ligand =  at_central["resName"] 
         
        # all atom ligand
        l_at_lig = loadFile.ExtractInfoPDBID(PDB_ID)[name_ligand][0] # change not tested
        l_at_subs = retrieveAtom.substructure (substruct, at_central, l_at_lig)
         
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
            


def EnvironmentInteraction (st_atom, pr_result, file_log, restrained = 1) : 
    
    d_neighbor_considered = structure.nbNeighbor()
    
    # first step fix if there are a salt bridges i.e. counter ion close to the substructure
    for type_subs in st_atom.keys () : 
        # create directory
        pr_count = pathManage.CreatePathDir(pr_result + str(type_subs) + "/Count/")
        pr_dist = pathManage.CreatePathDir(pr_result + str(type_subs) + "/Dist/")
        pr_angle = pathManage.CreatePathDir(pr_result + str(type_subs) + "/Angle/")
    
        # structure of count
        d_count_SB = {}
        d_count_noSB = {}
        d_dist = {}
        d_angle = {}
        
        for atom_central in st_atom[type_subs] : 
            l_interact = GetInteractions (atom_central["neighbors"], type_subs, restrained)
            
            # case -> stabilization by salt bridges
            nb_neighbor = d_neighbor_considered[type_subs]
            if type_subs == "COO" and "N" in l_interact or type_subs != "global" and "COO" in l_interact : 
                l_neighbors = deepcopy(atom_central["neighbors"])
                i = 1
                while i <= nb_neighbor :  
                    type_neighbor, atom_neighbor = searchMoreClose(l_neighbors, option_lcopy = 0)
                    
                    if not i in d_count_SB.keys () : 
                        d_count_SB[i] = structure.countClassificationAtoms()
                        d_count_SB[i]["out"] = 0 # case where there is not enough neighbor
                    if not i in d_dist.keys () : 
                        d_dist[i] = structure.EmptyListClassifAtom()
                    
                    if type_neighbor == None : 
                        d_count_SB[i]["out"] = d_count_SB[i]["out"] + 1
                    else : 
                        d_count_SB[i][type_neighbor] = d_count_SB[i][type_neighbor] + 1
                        d_dist[i][type_neighbor].append (atom_neighbor["distance"])
                        
                        # for angle
                        if i == 1 : 
                            d_angle_temp = {}
                            d_angle_temp[i] = atom_neighbor
                        else : 
                            d_angle_temp[i] = atom_neighbor
                            
                            # init angle dico
                            k_in = str ("angle" + str (i-1) + "_" + str (i))
                            if not k_in in d_angle.keys () : 
                                d_angle[k_in] = {}
                                d_angle[k_in]["distance"] = []
                                d_angle[k_in]["angle"] = []
                                d_angle[k_in]["type"] = [] 
                                
                            d_angle[k_in]["distance"].append (calcul.distanceTwoatoms(d_angle_temp[i-1], d_angle_temp[i]))
                            d_angle[k_in]["angle"].append (calcul.angleVector(d_angle_temp[i-1], atom_central, d_angle_temp[i]))
                            d_angle[k_in]["type"].append (str (structure.classificationATOM (d_angle_temp[i-1])) + "_" + str (structure.classificationATOM (d_angle_temp[i])))
                    i = i +1
            
            else : 
                l_neighbors = deepcopy(atom_central["neighbors"])
                i = 1
                while i <= nb_neighbor :  
                    type_neighbor, atom_neighbor = searchMoreClose(l_neighbors, option_lcopy = 0)
                    
                    k_in = str (i) + "-no"
                    if not k_in in d_count_noSB.keys () : 
                        d_count_noSB[k_in] = structure.countClassificationAtoms()
                        d_count_noSB[k_in]["out"] = 0 # case where there is not enough neighbor
                    
                    if type_neighbor == None : 
                        d_count_noSB[k_in]["out"] = d_count_noSB[k_in]["out"] + 1
                    else : 
                        d_count_noSB[k_in][type_neighbor] = d_count_noSB[k_in][type_neighbor] + 1
                        
                    i = i +1
        
        
        
        # statistic of substructure with salt bridge
        if d_count_SB != {} : 
            p_file_count = writeFile.CountByNeighbors(d_count_SB, pr_count + str (type_subs))
            runScriptR.AFC(p_file_count, nb_neighbor)
    
        if d_dist != {} : 
            d_file_dist = writeFile.DistByType(d_dist, type_subs, pr_dist)
            for neighbor in d_file_dist.keys () : 
                for type_atom in d_file_dist[neighbor].keys () : 
                    p_filin = d_file_dist[neighbor][type_atom].name
                    if type_atom == "density" : 
                        if path.getsize(p_filin) > 30 : 
                            runScriptR.plotDistanceDensity(p_filin, file_log, debug = 1)
                    else : 
                        if path.getsize(p_filin) > 10 : 
                            runScriptR.plotDistance(p_filin, file_log)
            
            
        if d_angle != {} : 
            l_p_file_angle = writeFile.AngleByType(d_angle, pr_angle)
            
            for p_file_angle in l_p_file_angle : 
                runScriptR.DistVSAngleNeighbor (p_file_angle)
        
        
        # AFC with substructure no stabilized by SB
        if d_count_noSB != {} : 
            p_file_count = writeFile.CountByNeighbors(d_count_noSB, pr_count + "no" + str (type_subs))
            runScriptR.AFC(p_file_count, nb_neighbor)
        
        # merge d_count with and without SB
        d_count_SBnoSB = MergeDicCount (d_count_SB, d_count_noSB)
        p_file_count = writeFile.CountByNeighbors(d_count_SBnoSB, pr_count + "SBandnoSB" + str (type_subs))
        runScriptR.AFC(p_file_count, nb_neighbor)
        

def MergeDicCount (d_count1, d_count2):
    
    d_out = deepcopy(d_count1)
    
    for k2 in d_count2.keys () :
        d_out[k2] = deepcopy(d_count2[k2])
    
    return d_out
    
