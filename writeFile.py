import structure
import repertory
from os import makedirs, path, mkdir
import tool
import numpy


def resultFilterLigandPDB(struct, dir_out):

    list_path_file_dataset = []
    for resolutionKey in struct.keys():
        # open filout dataset
        fileWriteXRay = open(dir_out + "dataset_" + resolutionKey, "w")
        list_path_file_dataset.append (dir_out + "dataset_" + resolutionKey)
        if resolutionKey != "NMR" and resolutionKey != "OUT":
            fileWriteXRayRMN = open(dir_out + "dataset_" + resolutionKey + "_RMN", "w")
            for ligandNMR in struct["NMR"].keys():
                if not ligandNMR in struct[resolutionKey]:
                    fileWriteXRayRMN.write(str(ligandNMR) + "\t")
                    for PDBNMR in struct["NMR"][ligandNMR]:
                        fileWriteXRayRMN.write(str(PDBNMR) + " ")
                    fileWriteXRayRMN.write("\n")
        for ligand in struct[resolutionKey].keys():
            fileWriteXRay.write(str(ligand) + "\t")
            if resolutionKey != "NMR" and resolutionKey != "OUT":
                fileWriteXRayRMN.write(str(ligand) + "\t")
            for PDB in struct[resolutionKey][ligand]:
                fileWriteXRay.write(str(PDB) + " ")
                if resolutionKey != "NMR" and resolutionKey != "OUT":
                    fileWriteXRayRMN.write(str(PDB) + " ")
            if resolutionKey != "NMR" and resolutionKey != "OUT":
                if ligand in struct["NMR"]:
                    for PDBMNR in struct["NMR"][ligand]:
                        fileWriteXRayRMN.write(str(PDBMNR) + " ")
            fileWriteXRay.write("\n")
            if resolutionKey != "NMR" and resolutionKey != "OUT":
                fileWriteXRayRMN.write("\n")
        fileWriteXRay.close()
        if resolutionKey != "NMR" and resolutionKey != "OUT":
            fileWriteXRayRMN.close()
        
    return list_path_file_dataset



def resultLigandInPDB(structResult, directory_out):

    filout = open(directory_out + "resultLigandInPDB", "w")

    for pdb in structResult:
        filout.write(pdb["name"] + "\t")
        if pdb["ligands"] != []:
            for ligand in pdb["ligands"]:
                filout.write(ligand + " ")
        filout.write("\n")
    filout.close()
    return directory_out + "resultLigandInPDB"



def listFloat (l_value, p_filin):
    
    filout = open (p_filin, "w")
    for v in l_value : 
        filout.write("%.3f\n" % v)
    filout.close ()
    return p_filin
    
    

def parsingDataSet(listCount, countAmine, numberPDB, path_file_dataset):

    filout = open(path_file_dataset  + ".stat", "w")
    maxRepetition = 0
    SumPDB = 0.0
    
    for count in listCount:
        filout.write(str(count["name"]) + "\t")
        filout.write(str(count["Number PDB"]))
        filout.write("\n")
        SumPDB = SumPDB + count["Number PDB"]
        if count["Number PDB"] > maxRepetition : 
            maxRepetition = count["Number PDB"] 
    
    nbCount = len(listCount) 
    filout.write("----------------------------------------\n")
    filout.write("Number ligand: " + str(nbCount) + "\n")
    filout.write("Number different PDB files: " + str(numberPDB) + "\n")
    filout.write("Repetition maximun ligand: " + str(maxRepetition) + "\n")
    filout.write("Average PDB by ligand: ")
    try : average = SumPDB / nbCount
    except : average = 0
    filout.write("%.3f\n" % average)
    filout.write("----------------------------------------\n")
    filout.write("Number Primary Amine: " + str(countAmine["Primary"]) + "\n")
    filout.write("Number Secondary Amine: " + str(countAmine["Secondary"]) + "\n")
    filout.write("Number Tertiary Amine: " + str(countAmine["Tertiary"]) + "\n")
    filout.write("Number Imidazole: " + str(countAmine["Imidazole"]) + "\n")
    filout.write("Number Guanidium: " + str(countAmine["Guanidium"]) + "\n")
    filout.write("Number Diamine: " + str(countAmine["Diamine"]) + "\n")
    filout.write("Number Pyridine: " + str(countAmine["Pyridine"]) + "\n")
    filout.write("Number Acid Carboxylic: " + str(countAmine["AcidCarboxylic"]) + "\n")
    
    filout.close()



def resultDistance(count, pr_result):

    dir_result = repertory.resultDistance(pr_result)
    filout = open(dir_result + "resultDistanceOx", "w")

    for element in count.keys():
        for distance in count[element]:
            filout.write(str("%.2f\t" % distance) + element + "\n")

    filout.close()
    
    
    return dir_result + "resultDistanceOx" 



def neighborStruct(struct_neighbor, struct_global_neighbor, files):

    # global -> append structure
    struct_neighbor["global"] = struct_global_neighbor
    for type_search in struct_neighbor.keys():
        if struct_neighbor[type_search] == [] : 
            continue
        for atom_central in struct_neighbor[type_search]:
            
            print atom_central.keys (), type_search, "check"
                
            lineWrite = str(atom_central["PDB"]) + "\t" + str(atom_central["serial"]) + "/" + str(atom_central["resName"]) + "/" + str(atom_central["x"]) + "/" + str(atom_central["y"]) +  "/" + str(atom_central["z"]) + "\t"
            for neighbor in atom_central["neighbors"]:
                lineWrite = lineWrite + str(neighbor["serial"]) + " " + str(neighbor["resSeq"]) + " " + str(neighbor["element"]) + " " + str(neighbor["name"]) + " " + str(neighbor["resName"]) + " " + str("%.2f" % neighbor["distance"]) + " " +str("%.3f" % neighbor["x"]) + " " +str("%.3f" % neighbor["y"]) + " " +str("%.3f" % neighbor["z"])  
                for angle in neighbor["angle"]:
                    lineWrite = lineWrite + " " + str("%.2f" % angle)
                lineWrite = lineWrite + "//"
            lineWrite = lineWrite + "\n"
            files[type_search].write(lineWrite)
    del struct_neighbor["global"] 


def openFileSummary(pr_out):

    # control if file exsixt
    if not path.isdir(pr_out):
        mkdir(pr_out)
    
    
    listS = structure.listStructure()
    listS.append("global")

    dictFile = {}

    for element in listS:
        filout = open(pr_out + "neighbor_" + element + ".sum", "w")
        dictFile[element] = filout

    return dictFile



def closeFileSummary(dictFile):

    for key in dictFile.keys():
        dictFile[key].close()


# # # # # # # def countGlobalCount(dist_max, stCount, dir_out):
# # # # # # # 
# # # # # # # #     resultDistanceOx(stCount[str(dist_max)]["distanceOx"], dir_out)
# # # # # # #     resultLigand(stCount[str(dist_max)]["ligand"], dir_out)
# # # # # # #     resultAtom(stCount[str(dist_max)]["atom"], dir_out)
# # # # # # #     resultByAA(stCount[str(dist_max)]["byAA"], dir_out)
# # # # # # #     resultAngle(stCount, dist_max, dir_out)
# # # # # # #     resultNeighbor (stCount[str(dist_max)]["threeAnalysis"], dir_out)
# # # # # # #     resultNumberNeighbor (stCount[str(dist_max)]["numberNeighbors"], dir_out)
# # # # # # # 
# # # # # # #     for distance in stCount.keys():
# # # # # # #         resultResidue(float(distance), stCount[str(distance)]["residue"], dir_out)
# # # # # # #         resultProportion(float(distance), stCount[str(distance)]["proportionAtom"], dir_out)
# # # # # # #         resultProportionType(float(distance), stCount[str(distance)]["proportionType"], dir_out)
# # # # # # # 
# # # # # # #         
# # # # # # #     resultProportionGlobalType(stCount, dist_max, dir_out)
# # # # # # #     resultGlobalResidue(stCount, dist_max, dir_out)
# # # # # # #     resultResidueDistance(stCount, dist_max, dir_out)
    

def resultCount(stCountLigand, type_count, pr_result):

    l_file_result = []
    for interestGroup in stCountLigand.keys ():
        p_filout = pr_result + str (interestGroup) + "_" + str(type_count)
        filout = open(p_filout, "w")
        
        for ligand in stCountLigand[interestGroup].keys():
            filout.write(str(ligand) + "\t" + str(stCountLigand[interestGroup][ligand]) + "\n")
        filout.close()
        l_file_result.append (p_filout)
    
    return l_file_result


def disributionNumberNeighbor (stCount, pr_result) : 
    
    l_p_filout = []
    
    for subs in stCount.keys () : 
        for distance in stCount[subs].keys () : 
            p_filout = pr_result + subs + distance
            l_p_filout.append (p_filout)
            filout = open (p_filout, "w")
            print stCount[subs][distance]
            filout.write ("\n".join ([str(i) for i in stCount[subs][distance]]))
        filout.close ()
    return l_p_filout
    
    

# def resultNumberNeighbor (count, directory_result) : 
#     
#     # define new fold
#     try : makedirs(directory_result + "NumberNeighbor/" , mode=0777)
#     except : pass
#     
#     filout_means = open(directory_result + "NumberNeighbor/means.result", "w")
#     for substruct in count.keys () : 
#         filout_means.write (str (substruct) + " " + str (numpy.mean(count[substruct])) + " " + str (numpy.std(count[substruct])) + "\n")
#         filout = open(directory_result + "NumberNeighbor/" + substruct, "w")
#         for nb_neighbor in count[substruct] : 
#             filout.write (str (nb_neighbor) + "\n")
#         filout.close ()
#     filout_means.close ()



# # # # # # # # def resultAtom(count, directory_result):
# # # # # # # # 
# # # # # # # #     listStructure = structure.listStructure()
# # # # # # # # 
# # # # # # # #     for element in listStructure:
# # # # # # # #         directory_in = repertory.typeSubStructure(directory_result, element)
# # # # # # # #         filout = open(directory_in + "statAtoms" + element, "w")
# # # # # # # #         for atom in count[element].keys():
# # # # # # # #             filout.write(str(atom) + "\t" + str(count[element][atom]) + "\n")
# # # # # # # # 
# # # # # # # #         filout.close()


def resultCountAA(stCountAA, pr_result):

    l_file_result = []
    for interestGroup in stCountAA.keys ():
        p_filout = pr_result + str (interestGroup) + "resProx"
        filout = open(p_filout, "w")
        
        for res in stCountAA[interestGroup].keys():
            filout.write(str(res) + "\t" + str(stCountAA[interestGroup][res]["main"]) + "\t" + str(stCountAA[interestGroup][res]["side"]) + "\n")
        filout.close()
        l_file_result.append (p_filout)
    
    return l_file_result


def resultByAA(stCount, max_distance ,pr_result):

    l_p_filout = []
    
    for substruct in stCount.keys ():
        pr_sub = repertory.resultSub(substruct, pr_result)
        for res in stCount[substruct].keys():
            p_filout = pr_sub + str (res)
            l_p_filout.append (p_filout)
            filout = open (p_filout, "w")
            
            print stCount[substruct][res]["<3.5"], res
            print stCount[substruct][res][">3.5"], res
            
            for name_atom in stCount[substruct][res]["<3.5"].keys():
                lineWrite = str(name_atom) + "\t" + str(stCount[substruct][res]["<3.5"][name_atom]) + "\t" + str(stCount[substruct][res][">3.5"][name_atom]) + "\n"
                filout.write(lineWrite)
            filout.close()

    return l_p_filout



def proportionByPositionNeighbors (stCount, pr_result):

    l_out = []    
    l_typeatom = structure.classificationATOM("", out_list= 1)

    # directory result
    for substruct in stCount.keys() : 
        p_filout = pr_result + "proportion_" + substruct
        l_out.append (p_filout)
        
        filout = open(p_filout, "w")
        filout.write ("\t".join(l_typeatom) + "\n")
        for nb_neighbor in range(1,8) : 
            if not nb_neighbor in stCount[substruct].keys () : 
                continue
            else : 
                filout.write(str(nb_neighbor))
                print substruct, nb_neighbor
                print stCount[substruct][nb_neighbor]
                sum_neighbor = tool.sumDict(stCount[substruct][nb_neighbor])
                for class_atom in l_typeatom :
                    filout.write("\t" + str(stCount[substruct][nb_neighbor][class_atom] / sum_neighbor)) 
                filout.write("\n")
        filout.close ()
    return l_out
            
            
def countFirstNeighbor (stCount, pr_result):
    
    l_typeatom = structure.classificationATOM("", out_list= 1)   
    filout = open (pr_result + "countFirst", "w")
    filout.write ("\t".join(l_typeatom) + "\n")
    
    for sub_struct in stCount.keys() : 
        filout.write (sub_struct)
        for class_atom in l_typeatom : 
            filout.write("\t" + str(stCount[sub_struct][1][class_atom])) # first neighbors
        filout.write("\n")
    filout.close ()
    return [pr_result + "countFirst"]



def distanceCountStruct(stCount, pr_result) : 
    
    l_filout = []
    
    print 
    for sub_struct in stCount.keys() : 
        filout = open (pr_result + "DistanceFirst" + str (sub_struct) + ".txt", "w")
        l_filout.append (pr_result + "DistanceFirst" + str (sub_struct) + ".txt")
        nb_fisrt = len (stCount[sub_struct][1]["distance"])
        
        print stCount[sub_struct][1]["distance"]
        
        i = 0
        while i < nb_fisrt : 
            print i
            filout.write (str(stCount[sub_struct][1]["distance"][i]) + "\t" + str (stCount[sub_struct][1]["classe"][i]) + "\n")
            print stCount[sub_struct][1]["distance"][i], stCount[sub_struct][1]["classe"][i]
            i = i + 1
        filout.close ()
        
    return l_filout
    




def countNeighborsAll(stCount, pr_result):
    
    l_typeatom = structure.classificationATOM("", out_list= 1)   
    filout = open (pr_result + "countAll", "w")
    filout.write ("\t".join(l_typeatom) + "\n")
    
    for sub_struct in stCount.keys() : 
        filout.write (sub_struct)
        for class_atom in l_typeatom : 
            count = 0
            for i_neighbor in stCount[sub_struct].keys () :
                print i_neighbor, ">-----<" 
                if type(stCount[sub_struct][i_neighbor]) == dict and class_atom in stCount[sub_struct][i_neighbor].keys () : 
                    count = count + stCount[sub_struct][i_neighbor][class_atom]
            filout.write("\t" + str(count)) # first neighbors
        filout.write("\n")
    filout.close ()
    return [pr_result + "countAll"]



        
def resultNeighbor (countStruct, pr_result) : 
    """
    Three neighbors analysis -> write files
    """
    # distance list
    l_typeatom = structure.classificationATOM("", out_list= 1)
    
    # directory result
    for sub_struct in countStruct.keys() : 
        filout_neighbor = open (pr_result + "neighbor_" + sub_struct, "w")
        filout_neighbor_count = open (pr_result + "neighbor_count_" + sub_struct, "w")
        filout_distance = open (pr_result + "distance_" + sub_struct, "w")
        filout_angle = open (pr_result + "angle_neighbors" + sub_struct, "w")
        filout_neighbor.write ("\t".join(l_typeatom) + "\n")
        filout_neighbor_count.write ("\t".join(l_typeatom) + "\n")
        # barplot class of neighbors -> but not dynamic nb neighbor
        for nb_neighbor in range(1,8) : 
            if nb_neighbor == "angle1_2" or nb_neighbor == "angle2_3" or nb_neighbor == "angle1_3" : 
                continue
            filout_neighbor.write(str(nb_neighbor))
            filout_neighbor_count.write(str(nb_neighbor))
            filout_distance.write(str(nb_neighbor))
            sum_neigbor = tool.sumDict(countStruct[sub_struct][nb_neighbor])
            for class_atom in l_typeatom : 
                filout_neighbor.write("\t" + str(countStruct[sub_struct][nb_neighbor][class_atom] / sum_neigbor)) 
                filout_neighbor_count.write("\t" + str(countStruct[sub_struct][nb_neighbor][class_atom])) 
            filout_neighbor.write("\n")
            filout_neighbor_count.write("\n")
            filout_distance.write("\t" + "\t".join(countStruct[sub_struct][nb_neighbor]["distance"]) + "\n")
            filout_distance.write ("Classe\t" + "\t".join(countStruct[sub_struct][nb_neighbor]["classe"]) + "\n")
    
    
        # angles between neighbors
        nb_angle = len (countStruct[sub_struct]["angle1_2"])
        i = 0
        while i < nb_angle : 
            filout_angle.write (str (countStruct[sub_struct]["angle1_2"][i]) + "\t" +str (countStruct[sub_struct]["angle1_3"][i]) + "\t" +str (countStruct[sub_struct]["angle2_3"][i])  + "\n" )
            i = i + 1
    
    filout_distance.close ()
    filout_neighbor.close ()
    filout_neighbor_count.close ()
    filout_angle.close ()
    
    # write barplot file
    barplotThreeAtomBarplot (countStruct, pr_result)   


def barplotThreeAtomBarplot (countStruct, dir_out):
    """
    Barplot for distance function type atoms
    """
    l_typeatom = structure.classificationATOM("", out_list= 1)
    
    for substruct in countStruct.keys () : 
        for nb_neighbor in countStruct[substruct].keys() :
            if type (nb_neighbor) != type(int()) or countStruct[substruct][nb_neighbor]["distance"] == []:
                continue 
            filout = open (dir_out + "barplot_" + substruct + "_" + str(nb_neighbor), "w")
            
            # header
            filout.write ("\t".join(l_typeatom) + "\n")
            
            #count
            d_cout = {}
#             min_distance = min(countStruct[substruct][nb_neighbor]["distance"])
            max_distance = float(max(countStruct[substruct][nb_neighbor]["distance"]))
            
            d_temp = 2
            l_dist = []
            
            while d_temp <= max_distance + 0.4 : 
                l_dist.append (d_temp)
                
                d_cout[d_temp] = {}
                for classe_atom in l_typeatom : 
                    d_cout[d_temp][classe_atom] = 0
                
                # implement count struct
                i = 0
                len_neighbor = len (countStruct[substruct][nb_neighbor]["classe"])
                print d_temp
                while i < len_neighbor : 
                    
                    if float(countStruct[substruct][nb_neighbor]["distance"][i]) <= d_temp  and float(countStruct[substruct][nb_neighbor]["distance"][i]) > d_temp - 0.2 : 
                        d_cout[d_temp][countStruct[substruct][nb_neighbor]["classe"][i]] = d_cout[d_temp][countStruct[substruct][nb_neighbor]["classe"][i]] + 1
                    else : 
                        pass
                    
                    i = i + 1
                    
                d_temp = d_temp + 0.2

            
            for dist in l_dist : 
                filout.write (str (dist))
                for class_atom in l_typeatom : 
                    filout.write ("\t" + str(d_cout[dist][class_atom]))
                filout.write ("\n")
        filout.close ()
    


# def resultProportion (distance, count, directory_result):
#     """Write file proportion number of neighbors"""
# 
# 
#     for type_substruct in count.keys():
#         dir_in = repertory.globalProportionAtom(directory_result)
#         filout = open (dir_in + "proportionAtom" + type_substruct + str("%.2f" % distance), "w")
#         for nbNeighbor in count[type_substruct].keys():
#             filout.write(str(nbNeighbor) + "\t" + str(count[type_substruct][nbNeighbor]["C"]) + "\t" + str(count[type_substruct][nbNeighbor]["O"]) + "\t" + str(count[type_substruct][nbNeighbor]["N"]) + "\t" + str(count[type_substruct][nbNeighbor]["S"]) + "\t" + str(count[type_substruct][nbNeighbor]["others"]) + "\n")
# 
# 
#         filout.close()
        
        
# def resultProportionType (distance, count, directory_result):
# 
#     listClasse = structure.classificationATOM("", out_list= 1)
#     for type in count.keys():
#         dir_in = repertory.globalProportionType(directory_result)
#         filout = open (dir_in + "proportionType" + type + str("%.2f" % distance), "w")
#         filout.write ("\t".join(listClasse) + "\n")
#         for nbNeighbor in count[type].keys():
#             if not nbNeighbor == "allNumberNeighbors" : 
#                 filout.write(str(nbNeighbor))
#                 
#                 for classe in listClasse : 
#                     filout.write("\t" + str(count[type][nbNeighbor][classe]))
#                 filout.write("\n")
#                 
#         if count[type].keys() == [] : 
#             filout.write(str(0) + "\t" + str(0) + "\t" + str(0) + "\t" + str(0) + "\t" + str(0) + "\t" + str(0) + "\n")
#         
#         filout.close()

# def resultProportionGlobalType(count, distanceGlobal,directory_result) : 
#     
#     listDistance = structure.listDistance(distanceGlobal)
#     listClasse = structure.classificationATOM("", out_list=1)
#     listStruct = structure.listStructure()
#     listStruct.append("Global") 
#     
#     for type_substructure in listStruct :
#         dir_in = repertory.globalProportionType(directory_result) 
#         filout = open (dir_in + "proportionType" + type_substructure , "w")
#         filout.write("\t".join(listClasse) + "\n")
#         for distance in listDistance :
#             filout.write(str(distance))
#             for classe in listClasse : 
#                 filout.write("\t" + str(count[distance]["proportionType"][type_substructure]["allNumberNeighbors"][classe]))
#             filout.write("\n")
#             
#         filout.close()
#     
    

    

# def resultResidueDistance(countGlobalAmine, distanceMax, directory_out):
# 
#     listDistance = structure.listDistance(distanceMax)
#     listAminoAcid = ["HOH", "ASP", "GLU", "THR", "SER", "ASN", "GLN", "TYR", "HIS", "LYS", "ARG", "PHE", "TRP", "ALA", "ILE", "LEU", "MET", "VAL", "CYS", "GLY", "PRO"]
# 
# 
#     for type_substructure in  countGlobalAmine["2.5"]["residue"].keys():
#         filout_side_chain = open(directory_out + str(type_substructure) + "/GlobalResidueSide" + str(type_substructure), "w")
#         filout_global = open(directory_out + str(type_substructure) + "/GlobalResidue" + str(type_substructure), "w")
#         for aminoAcid in listAminoAcid:
#             line_side = str(aminoAcid)
#             line_global = str(aminoAcid)
#             for distance in listDistance:
#                 line_side = line_side + "\t"
#                 line_global = line_global + "\t"
#                 if aminoAcid == "HOH":
#                     line_side = line_side + str(countGlobalAmine[distance]["residue"][type_substructure][aminoAcid]["main"])
#                     line_global = line_global + str(countGlobalAmine[distance]["residue"][type_substructure][aminoAcid]["main"])
#                 else:
#                     line_side = line_side + str(countGlobalAmine[distance]["residue"][type_substructure][aminoAcid]["side"])
#                     line_global = line_global + str(countGlobalAmine[distance]["residue"][type_substructure][aminoAcid]["side"] + countGlobalAmine[distance]["residue"][type_substructure][aminoAcid]["main"])
# 
#             line_side = line_side + "\n"
#             line_global = line_global + "\n"
#             filout_side_chain.write(line_side)
#             filout_global.write(line_global)
#         filout_side_chain.close()
#         filout_global.close ()



# def resultAtLeastOneGlobal(countGlobal, distanceMax, directory_out):
# 
#     listDistance = structure.listDistance(distanceMax)
#     listStudyAtLeastOne = countGlobal[str(distanceMax)]["atLeastOne"].keys()
#     
#     for StudyAtleastOne in listStudyAtLeastOne : 
#         filout_side_chain = open(directory_out + "atLeastOneGlobal_" + StudyAtleastOne, "w")
#         line_side = ""
#         for distance in listDistance :  
#             line_side = line_side + str(countGlobal[distance]["atLeastOne"][StudyAtleastOne][StudyAtleastOne]) + "\t" + str(countGlobal[distance]["atLeastOne"][StudyAtleastOne]["others"]) + "\t"
#         line_side = line_side + "\n" 
#         
#         filout_side_chain.write(line_side)
#         filout_side_chain.close()


def resultResProx(stCount, distance_max, pr_result):

    l_p_filout = []
    l_distance = structure.listDistance(distance_max)
    l_aa = ["HOH", "ASP", "GLU", "THR", "SER", "ASN", "GLN", "TYR", "HIS", "LYS", "ARG", "PHE", "TRP", "ALA", "ILE", "LEU", "MET", "VAL", "CYS", "GLY", "PRO"]
    
    for substruct in stCount.keys () : 
        print substruct, "-----"
        print stCount[substruct]
        p_filout = pr_result + substruct + "resCount"
        l_p_filout.append (p_filout)
        filout = open (p_filout, "w")
        for aa in l_aa : 
            line_w = str (aa)
            for distance in l_distance : 
                line_w = line_w + "\t" + str(stCount[substruct][distance][aa])
            line_w = line_w + "\n"
            filout.write (line_w)
        filout.close ()
    
    return l_p_filout
    

# def resultAtLeastOne(count_global, count_substructure, max_distance, directory_out):
# 
#     # distance list
#     listDistance = structure.listDistance(max_distance)
#     l_study =  structure.listStructure()
#     l_study.append("global")
# #     print l_study
#     
#     # directory result
#     dir_out = directory_out + "AtLeastOne/"
#     try : 
#         makedirs(dir_out, mode=0777)
#     except : 
#         pass
#     
#     # type at least one
# #     print count_substructure[listDistance[0]]["atLeastOne"].keys()
#     for type_atleastone in count_substructure[listDistance[0]]["atLeastOne"].keys() : 
#         filout = open (dir_out + type_atleastone + ".dat", "w")
#         
#         # header
#         for distance in listDistance : 
#             filout.write(str(distance) + "\t")
#         
#         # rownames
#         filout.write("\n")
#         
#         for type_study in l_study : 
#             filout.write (str(type_study) + "\t")
#         
#             # data
#             for distance in listDistance :
#                 # print type_study
#                 if type_study == "global" : 
#                     filout.write (str(count_global[distance]["atLeastOne"][type_atleastone][type_atleastone] / (count_global[distance]["atLeastOne"][type_atleastone][type_atleastone] +count_global[distance]["atLeastOne"][type_atleastone]["other"]) ) + "\t")
#                 else : 
#                     filout.write (str(count_substructure[distance]["atLeastOne"][type_atleastone][type_study][type_atleastone] / (count_substructure[distance]["atLeastOne"][type_atleastone][type_study][type_atleastone] +count_substructure[distance]["atLeastOne"][type_atleastone][type_study]["other"]) ) + "\t")
#             filout.write("\n")
#         filout.close()



def resultAngle(d_count, pr_out):
    
    l_p_file = []
    for type_substruct in d_count.keys():
        pr_final = repertory.resultAngle(pr_out, type_substruct)
        p_filout = pr_final + "angle_" + str(type_substruct)
        l_p_file.append (p_filout)
        filoutGlobal = open(p_filout, "w")
        for classe in  d_count[type_substruct].keys():
            nbDistance = len(d_count[type_substruct][classe]["distance"])
            for i in range(0, nbDistance) : 
                distanceAt = d_count[type_substruct][classe]["distance"][i]
                filoutGlobal.write("%.2f" % distanceAt)
                for angle in d_count[type_substruct][classe]["angles"][i] : 
                    filoutGlobal.write("\t%.2f" % angle)
                filoutGlobal.write("\t" + classe + "\n")
        filoutGlobal.close()
    return l_p_file
 

def dAngleType (count, directory_in):
    
    listClasse = structure.classificationATOM("", out_list = 1)
    for type_substruct in count.keys() : 
        for distance in count[type_substruct].keys() : 
            pr_angle_type = repertory.resultAngle(directory_in, type_substruct)
            filout = open(pr_angle_type + "angle_" + type_substruct + "_" + distance, "w")
            
            ### HEADER ###
            filout.write (str(listClasse[0]))
            for classe in listClasse[1:] : 
                filout.write ("\t" + classe)
            filout.write("\n")
            ##############
            
            for angle in count[type_substruct][distance].keys() : 
                filout.write(str(angle))
                for classe in listClasse : 
                    try : filout.write("\t" + str(count[type_substruct][distance][angle][classe]))
                    except : filout.write("\t0")
                filout.write("\n")
        filout.close()    


def openFilesWithoutSummary(distanceMax, directory_in):
    
    fileClass = {}
    
    listDistance = structure.listDistance(distanceMax)
    listStructureStudy = structure.listStructure()
    
    for distance in listDistance : 
        fileClass[distance] = {}
        for struct in listStructureStudy :
            dir_struct = repertory.withoutAtLeastOneSummary(directory_in)
            fileClass[distance][struct] = open(dir_struct + struct + "_<" + distance, "w")
    
    return fileClass
    
    
    
def withoutAtLeastOneSummary(stAtom, filesWithout, distance):
    
    listStudy = stAtom.keys()
    
    for study in listStudy : 
        for nitrogen in stAtom[study] : 
            flag = 0 
            for neighbors in nitrogen["neighbors"] : 
                if neighbors["classificationAtLeastOne"] != "others" : 
                    flag = 1
            if flag == 0 : 
                filesWithout[str(distance)][study].write(str(nitrogen["PDB"]) + "\t" + str(nitrogen["serial"]) + "\n")


def closeFilesWithoutSummary (filesWithoutSummary):
    
    for distance in filesWithoutSummary.keys() : 
        for study in filesWithoutSummary[distance].keys() : 
            filesWithoutSummary[distance][study].close()

           

                

def lenBondType (l_distance, l_first, path_filout) : 
    
    filout = open (path_filout, "w")
    
    nb_element = len (l_distance)
    # control same size of list
    if nb_element != len (l_first) : 
        filout.close()
        print "ERROR number neighbors -writeFile l.547"
        return
    else :
        i = 0
        while i < nb_element : 
            filout.write (str (l_distance[i]) + "\t" + str (l_first[i]) + "\n")
            i = i + 1
        
    filout.close ()
        
        

def coordinates3D (l_atom, p_filout, type_substruct) : 
    
    filout = open (p_filout, "w") 
    
    substruct =  structure.substructureCoord(type_substruct)
    
    for atom in substruct : 
        filout.write (str(atom["x"]) + "\t" + str(atom["y"]) + "\t" + str(atom["z"]) + "\t" + "REF" + "\n")
    
    for atom in l_atom : 
        if not "occupancy" in atom.keys () : 
            filout.write (str(atom["x"]) + "\t" + str(atom["y"]) + "\t" + str(atom["z"]) + "\t" + structure.classificationATOM (atom) + "\n")
    filout.close () 
    
    return  p_filout    
        
        
        
        
        
        
        


