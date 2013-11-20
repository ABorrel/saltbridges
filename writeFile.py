import structure
import repertory
from os import makedirs
import tool



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



def resultCoplanar(struct, name, directory_out):
    
    filout = open(directory_out + name, "w")

    for element in struct:
        filout.write("%.2f\n" % element)
    filout.close()



def resultLengthCNBond(structResult, file, directory_out):

    filout = open(directory_out + file, "w")

    for element in structResult:
        filout.write("%.2f\n" % element)

    filout.close()



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



def resultDistanceOx(count, directory_in):

    dir_result = repertory.resultDistance(directory_in)
    filout = open(dir_result + "resultDistanceOx", "w")

    for element in count.keys():
        for distance in count[element]:
            filout.write(str("%.2f\t" % distance) + element + "\n")

    filout.close()



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


def openFileSummary(directory_out):

    listS = structure.listStructure()
    listS.append("global")

    dictFile = {}

    for element in listS:
        filout = open(directory_out + "neighbor_" + element + ".sum", "w")
        dictFile[element] = filout

    return dictFile



def closeFileAmine(dictFile):

    for key in dictFile.keys():
        dictFile[key].close()


def countGlobalCount(distanceGlobal, countGlobalAmine, dir_out):

    resultDistanceOx(countGlobalAmine[str(distanceGlobal)]["distanceOx"], dir_out)
    resultLigand(countGlobalAmine[str(distanceGlobal)]["ligand"], dir_out)
    resultAtom(countGlobalAmine[str(distanceGlobal)]["atom"], dir_out)
    resultByAA(countGlobalAmine[str(distanceGlobal)]["byAA"], dir_out)
    resultAngle(countGlobalAmine, distanceGlobal, dir_out)
    resultNeighbor (countGlobalAmine[str(distanceGlobal)]["threeAnalysis"], dir_out)

    for distance in countGlobalAmine.keys():
        resultResidue(float(distance), countGlobalAmine[str(distance)]["residue"], dir_out)
        resultProportion(float(distance), countGlobalAmine[str(distance)]["proportionAtom"], dir_out)
        resultProportionType(float(distance), countGlobalAmine[str(distance)]["proportionType"], dir_out)

        
    resultProportionGlobalType(countGlobalAmine, distanceGlobal, dir_out)
    resultGlobalResidue(countGlobalAmine, distanceGlobal, dir_out)
    resultResidueDistance(countGlobalAmine, distanceGlobal, dir_out)
    

def resultLigand(count, directory_result):

    listStructure = structure.listStructure()

    for element in listStructure:
        directory_in = repertory.typeSubStructure(directory_result, element)
        filout = open(directory_in + "statLigands" + element, "w")
        
        for ligand in count[element].keys():
            filout.write(str(ligand) + "\t" + str(count[element][ligand]) + "\n")
        filout.close()


def resultAtom(count, directory_result):

    listStructure = structure.listStructure()

    for element in listStructure:
        directory_in = repertory.typeSubStructure(directory_result, element)
        filout = open(directory_in + "statAtoms" + element, "w")
        for atom in count[element].keys():
            filout.write(str(atom) + "\t" + str(count[element][atom]) + "\n")

        filout.close()


def resultResidue(distance, count, directory_result):

    listStructure = structure.listStructure()
    
    for element in listStructure:
        dir_in = repertory.typeSubStructure(directory_result, element)
        filout = open(dir_in + "statResidues" + element + str("%.2f" % distance), "w")
        for residue in count[element].keys():
            filout.write(str(residue) + "\t" + str(count[element][residue]["main"]) + "\t" + str(count[element][residue]["side"]) + "\n")

        filout.close()


def resultByAA(count, directory_result):

    listStructure = structure.listStructure()

    for element in listStructure:
        for aminoAcid in count[element].keys():
            dir_in = repertory.typeSubStructure(directory_result, element + "/aminoAcid/")
            filout = open(dir_in + element + aminoAcid, "w")
            for atom in count[element][aminoAcid].keys():
                lineWrite = str(atom) + "\t" + str(count[element][aminoAcid][atom]["3.5"]) + "\t" + str(count[element][aminoAcid][atom]["4.5"]) + "\n"
                filout.write(lineWrite)
            filout.close()

    for aminoAcid in count["global"].keys():
        dir_in = repertory.typeSubStructure(directory_result, "AminoAcidGlobal")
        filout = open(dir_in + "Global" + aminoAcid, "w")
        for atom in count["global"][aminoAcid].keys():
            lineWrite = str(atom) + "\t" + str(count["global"][aminoAcid][atom]["3.5"]) + "\t" + str(count["global"][aminoAcid][atom]["4.5"]) + "\n"
            filout.write(lineWrite)
        filout.close()


def resultProportion (distance, count, directory_result):
    """Write file proportion number of neighbors"""


    for type_substruct in count.keys():
        dir_in = repertory.globalProportionAtom(directory_result)
        filout = open (dir_in + "proportionAtom" + type_substruct + str("%.2f" % distance), "w")
        for nbNeighbor in count[type_substruct].keys():
            filout.write(str(nbNeighbor) + "\t" + str(count[type_substruct][nbNeighbor]["C"]) + "\t" + str(count[type_substruct][nbNeighbor]["O"]) + "\t" + str(count[type_substruct][nbNeighbor]["N"]) + "\t" + str(count[type_substruct][nbNeighbor]["S"]) + "\t" + str(count[type_substruct][nbNeighbor]["others"]) + "\n")


        filout.close()
        
        
def resultProportionType (distance, count, directory_result):

    listClasse = structure.classificationATOM("", out_list= 1)
    for type in count.keys():
        dir_in = repertory.globalProportionType(directory_result)
        filout = open (dir_in + "proportionType" + type + str("%.2f" % distance), "w")
        filout.write ("\t".join(listClasse) + "\n")
        for nbNeighbor in count[type].keys():
            if not nbNeighbor == "allNumberNeighbors" : 
                filout.write(str(nbNeighbor))
                
                for classe in listClasse : 
                    filout.write("\t" + str(count[type][nbNeighbor][classe]))
                filout.write("\n")
                
        if count[type].keys() == [] : 
            filout.write(str(0) + "\t" + str(0) + "\t" + str(0) + "\t" + str(0) + "\t" + str(0) + "\t" + str(0) + "\n")
        
        filout.close()

def resultProportionGlobalType(count, distanceGlobal,directory_result) : 
    
    listDistance = structure.listDistance(distanceGlobal)
    listClasse = structure.classificationATOM("", out_list=1)
    listStruct = structure.listStructure()
    listStruct.append("Global") 
    
    for type_substructure in listStruct :
        dir_in = repertory.globalProportionType(directory_result) 
        filout = open (dir_in + "proportionType" + type_substructure , "w")
        filout.write("\t".join(listClasse) + "\n")
        for distance in listDistance :
            filout.write(str(distance))
            for classe in listClasse : 
                filout.write("\t" + str(count[distance]["proportionType"][type_substructure]["allNumberNeighbors"][classe]))
            filout.write("\n")
            
        filout.close()
    
    
    

def resultResidueDistance(countGlobalAmine, distanceMax, directory_out):

    listDistance = structure.listDistance(distanceMax)
    listAminoAcid = ["HOH", "ASP", "GLU", "THR", "SER", "ASN", "GLN", "TYR", "HIS", "LYS", "ARG", "PHE", "TRP", "ALA", "ILE", "LEU", "MET", "VAL", "CYS", "GLY", "PRO"]


    for type_substructure in  countGlobalAmine["2.5"]["residue"].keys():
        filout_side_chain = open(directory_out + str(type_substructure) + "/GlobalResidueSide" + str(type_substructure), "w")
        filout_global = open(directory_out + str(type_substructure) + "/GlobalResidue" + str(type_substructure), "w")
        for aminoAcid in listAminoAcid:
            line_side = str(aminoAcid)
            line_global = str(aminoAcid)
            for distance in listDistance:
                line_side = line_side + "\t"
                line_global = line_global + "\t"
                if aminoAcid == "HOH":
                    line_side = line_side + str(countGlobalAmine[distance]["residue"][type_substructure][aminoAcid]["main"])
                    line_global = line_global + str(countGlobalAmine[distance]["residue"][type_substructure][aminoAcid]["main"])
                else:
                    line_side = line_side + str(countGlobalAmine[distance]["residue"][type_substructure][aminoAcid]["side"])
                    line_global = line_global + str(countGlobalAmine[distance]["residue"][type_substructure][aminoAcid]["side"] + countGlobalAmine[distance]["residue"][type_substructure][aminoAcid]["main"])

            line_side = line_side + "\n"
            line_global = line_global + "\n"
            filout_side_chain.write(line_side)
            filout_global.write(line_global)
        filout_side_chain.close()
        filout_global.close ()



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


def resultGlobalResidue(count, distanceMax, directory_out):

    listDistance = structure.listDistance(distanceMax)
    listAminoAcid = ["HOH", "ASP", "GLU", "THR", "SER", "ASN", "GLN", "TYR", "HIS", "LYS", "ARG", "PHE", "TRP", "ALA", "ILE", "LEU", "MET", "VAL", "CYS", "GLY", "PRO"]
    filout_side = open(directory_out + "globalResidueAllAtomsSide", "w")
    filout_all = open(directory_out + "globalResidueAllAtoms", "w")
    for aminoAcid in listAminoAcid:
        line_side = str(aminoAcid)
        line_global = str(aminoAcid)
        for distance in listDistance:
            if aminoAcid == "HOH":
                line_side = line_side + "\t" + str(count[distance]["ResidueAllAtom"][aminoAcid]["main"])
                line_global = line_global + "\t" + str(count[distance]["ResidueAllAtom"][aminoAcid]["main"])
            else:
                line_side = line_side + "\t" + str(count[distance]["ResidueAllAtom"][aminoAcid]["side"])
                line_global = line_global + "\t" + str(count[distance]["ResidueAllAtom"][aminoAcid]["side"] + count[distance]["ResidueAllAtom"][aminoAcid]["main"])

        line_side = line_side + "\n"
        line_global = line_global + "\n"
        filout_side.write(line_side)
        filout_all.write(line_global)


def resultAtLeastOne(count_global, count_substructure, max_distance, directory_out):

    # distance list
    listDistance = structure.listDistance(max_distance)
    l_study =  structure.listStructure()
    l_study.append("global")
#     print l_study
    
    # directory result
    dir_out = directory_out + "AtLeastOne/"
    try : 
        makedirs(dir_out, mode=0777)
    except : 
        pass
    
    # type at least one
#     print count_substructure[listDistance[0]]["atLeastOne"].keys()
    for type_atleastone in count_substructure[listDistance[0]]["atLeastOne"].keys() : 
        filout = open (dir_out + type_atleastone + ".dat", "w")
        
        # header
        for distance in listDistance : 
            filout.write(str(distance) + "\t")
        
        # rownames
        filout.write("\n")
        
        for type_study in l_study : 
            filout.write (str(type_study) + "\t")
        
            # data
            for distance in listDistance :
                # print type_study
                if type_study == "global" : 
                    filout.write (str(count_global[distance]["atLeastOne"][type_atleastone][type_atleastone] / (count_global[distance]["atLeastOne"][type_atleastone][type_atleastone] +count_global[distance]["atLeastOne"][type_atleastone]["other"]) ) + "\t")
                else : 
                    filout.write (str(count_substructure[distance]["atLeastOne"][type_atleastone][type_study][type_atleastone] / (count_substructure[distance]["atLeastOne"][type_atleastone][type_study][type_atleastone] +count_substructure[distance]["atLeastOne"][type_atleastone][type_study]["other"]) ) + "\t")
            filout.write("\n")
        filout.close()



def resultAngle(count, distanceMax, dir_in):
    
    for type_substruct in count[str(distanceMax)]["angle"].keys():
        dir_angle = repertory.resultAngle(type_substruct, dir_in)
        filoutGlobal = open(dir_angle + "angle_" + str(type_substruct), "w")
        
        for classe in  count[str(distanceMax)]["angle"][type_substruct].keys():
            nbDistance = len(count[str(distanceMax)]["angle"][type_substruct][classe]["distance"])
            for i in range(0, nbDistance) : 
                distanceAt = count[str(distanceMax)]["angle"][type_substruct][classe]["distance"][i]
                filoutGlobal.write("%.2f" % distanceAt)
                
                for angle in count[str(distanceMax)]["angle"][type_substruct][classe]["angles"][i] : 
                    filoutGlobal.write("\t%.2f" % angle)
                    
                filoutGlobal.write("\t" + classe + "\n")
        filoutGlobal.close()
    countType = structure.countAngleType(count)  
    countAngle(countType, dir_in) 
    

def countAngle (count, directory_in):
    
    listClasse = structure.classificationATOM("", out_list = 1)
    for type_substruct in count.keys() : 
        for distance in count[type_substruct].keys() : 
            dir_in = repertory.resultAngle(type_substruct, directory_in)
            filout = open(dir_in + "angle_" + type_substruct + "_" + distance, "w")
            
            ### HEADER ###
            filout.write (str(listClasse[0]))
            for classe in listClasse[1:] : 
                filout.write ("\t" + classe)
            filout.write("\n")
            ##############
            
            for angle in count[type_substruct][distance].keys() : 
                filout.write(str(angle))
                for classe in listClasse : 
                    filout.write("\t" + str(count[type_substruct][distance][angle][classe]))
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
    
    
    
def withoutAtLeastOneSummary(amine, filesWithout, distance):
    
    listStudy = amine.keys()
    
    for study in listStudy : 
        for nitrogen in amine[study] : 
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

           

def resultNeighbor (countStruct, dir_out) : 
    """
    Three neighbors analysis -> write files
    """
    # distance list
    l_typeatom = structure.classificationATOM("", out_list= 1)
    
    # directory result
    dir_out = dir_out + "neigbhor/"
    try : 
        makedirs(dir_out, mode=0777)
    except : 
        pass

    for sub_struct in countStruct.keys() : 
        filout_neighbor = open (dir_out + "neighbor_" + sub_struct, "w")
        filout_distance = open (dir_out + "distance_" + sub_struct, "w")
        filout_angle = open (dir_out + "angle_neighbors" + sub_struct, "w")
        filout_neighbor.write ("\t".join(l_typeatom) + "\n")
        # barplot class of neighbors
        for nb_neighbor in range(1,4) : 
            if nb_neighbor == "angle1_2" or nb_neighbor == "angle2_3" or nb_neighbor == "angle1_3" : 
                continue
            filout_neighbor.write(str(nb_neighbor))
            filout_distance.write(str(nb_neighbor))
            sum_neigbor = tool.sumDict(countStruct[sub_struct][nb_neighbor])
            for class_atom in l_typeatom : 
                filout_neighbor.write("\t" + str(countStruct[sub_struct][nb_neighbor][class_atom] / sum_neigbor)) 
            filout_neighbor.write("\n")
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
    filout_angle.close ()
    
    # write barplot file
    barplotThreeAtomBarplot (countStruct, dir_out)   


def barplotThreeAtomBarplot (countStruct, dir_out):
    """
    Barplot for distance function type atoms
    """
    l_typeatom = structure.classificationATOM("", out_list= 1)
    
    for substruct in countStruct.keys () : 
        for nb_neighbor in countStruct[substruct].keys() :
            if type (nb_neighbor) != type(int()) :
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
        
        
        
        
        
        
        


