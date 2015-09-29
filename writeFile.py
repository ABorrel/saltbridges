import structure
import pathManage
import writePDBfile
from os import makedirs, path, mkdir
import tool
import numpy


def resultFilterLigandPDB(d_dataset, dir_out):

    l_p_dataset = []
    for RX in d_dataset.keys():
        # open filout dataset
        filout_dataset = open(dir_out + "dataset_" + RX + ".txt", "w")
        l_p_dataset.append (dir_out + "dataset_" + RX + ".txt")
        for ligand in d_dataset[RX].keys():
            filout_dataset.write(str(ligand) + "\t")
            filout_dataset.write (" ".join(d_dataset[RX][ligand]) + "\n")
        filout_dataset.close()
        
    return l_p_dataset



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
    
    

def AnalysisDataSet(l_d_lig, d_count_sub, numberPDB, path_file_dataset):


    filout = open(path_file_dataset[0:-4]  + ".stat", "w")
    maxRepetition = 0
    SumPDB = 0.0
    
    for count in l_d_lig:
        filout.write(str(count["name"]) + "\t")
        filout.write(str(count["Number PDB"]))
        filout.write("\n")
        SumPDB = SumPDB + count["Number PDB"]
        if count["Number PDB"] > maxRepetition : 
            maxRepetition = count["Number PDB"] 
    
    nb_lig = len(l_d_lig) 
    filout.write("----------------------------------------\n")
    filout.write("Number ligand: " + str(nb_lig) + "\n")
    filout.write("Number different PDB files: " + str(numberPDB) + "\n")
    filout.write("Repetition maximun ligand: " + str(maxRepetition) + "\n")
    filout.write("Average PDB by ligand: ")
    try : average = SumPDB / nb_lig
    except : average = 0
    filout.write("%.3f\n" % average)
    filout.write("----------------------------------------\n")
    for sub in d_count_sub.keys () : 
        filout.write("Number of " + sub + ": " + str(d_count_sub[sub][sub]) + "\n")
        filout.write("Number of PDB for " + sub + ": " + str(d_count_sub[sub]["PDB"]) + "\n")
        filout.write("Number of ligand for " + sub + ": " + str(d_count_sub[sub]["ligand"]) + "\n")
    filout.close()



def resultDistance(count, pr_result):

    dir_result = pathManage.resultDistance(pr_result)
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
            lineWrite = str(atom_central["PDB"]) + "\t" + str(atom_central["serial"]) + "/" + str(atom_central["resName"]) + "/" + str(atom_central["x"]) + "/" + str(atom_central["y"]) +  "/" + str(atom_central["z"]) + "\t"
            for neighbor in atom_central["neighbors"]:
                lineWrite = lineWrite + str(neighbor["serial"]) + " " + str(neighbor["resSeq"]) + " " + str(neighbor["element"]) + " " + str(neighbor["name"]) + " " + str(neighbor["resName"]) + " " + str("%.2f" % neighbor["distance"]) + " " +str("%.3f" % neighbor["x"]) + " " +str("%.3f" % neighbor["y"]) + " " +str("%.3f" % neighbor["z"])  
                for angleSubs in neighbor["angleSubs"]:
                    try : lineWrite = lineWrite + " " + str("%.2f" % angleSubs)
                    except : lineWrite = lineWrite + " NA"
                lineWrite = lineWrite + "//"
            lineWrite = lineWrite + "\n"
            files[type_search].write(lineWrite)
    del struct_neighbor["global"] 


def openFileSummary(pr_out):

    # control if file exsixt
    if not path.isdir(pr_out):
        mkdir(pr_out)
    
    
    l_sub = structure.ListSub()
    l_sub.append("global")

    d_filin = {}

    for element in l_sub:
        filout = open(pr_out + "neighbor_" + element + ".sum", "w")
        d_filin[element] = filout

    return d_filin



def closeFileSummary(dictFile):

    for key in dictFile.keys():
        dictFile[key].close()



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
            filout.write ("\n".join ([str(i) for i in stCount[subs][distance]]))
        filout.close ()
    return l_p_filout
    
    


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
        pr_sub = pathManage.resultSub(substruct, pr_result)
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
            #print "***", stCount[sub_struct]
            try : filout.write("\t" + str(stCount[sub_struct][1][class_atom])) # first neighbors
            except : filout.write("\tNA") # first neighbors
        filout.write("\n")
    filout.close ()
    return [pr_result + "countFirst"]



def distanceCountStruct(stCount, pr_result) : 
    
    l_filout = []
    
    for sub_struct in stCount.keys() : 
        try : nb_fisrt = len (stCount[sub_struct][1]["distance"])
        except : continue

        filout = open (pr_result + "DistanceFirst" + str (sub_struct) + ".txt", "w")
        l_filout.append (pr_result + "DistanceFirst" + str (sub_struct) + ".txt")
        
        # print stCount[sub_struct][1]["distance"]
        
        i = 0
        while i < nb_fisrt : 
            # print i
            filout.write (str(stCount[sub_struct][1]["distance"][i]) + "\t" + str (stCount[sub_struct][1]["classe"][i]) + "\n")
            # print stCount[sub_struct][1]["distance"][i], stCount[sub_struct][1]["classe"][i]
            i = i + 1
        filout.close ()
        
    return l_filout


def CountByNeighbors (d_count, p_result):
    
    filout = open (p_result, "w")
    
    l_neighbor = d_count.keys ()
    l_type_atom = d_count[l_neighbor[0]].keys ()
    
    print l_type_atom
    print l_neighbor
    
    filout.write ("\t".join(l_type_atom) + "\n")
    
    for neighbor in l_neighbor : 
        filout.write (str (neighbor))
        for type_atom in l_type_atom : 
            filout.write ("\t" + str (d_count[neighbor][type_atom]))
        filout.write ("\n")
    
    filout.close ()
    
    return p_result
    
    
def CountByDist (d_count, p_result):
    
    
    filout = open (p_result, "w")
    
    l_neighbor = d_count.keys ()
    l_type_atom = d_count[l_neighbor[0]].keys ()
    
    print l_type_atom
    print l_neighbor
    
    filout.write ("\t".join(l_type_atom) + "\n")
    
    for neighbor in l_neighbor : 
        filout.write (str (neighbor))
        for type_atom in l_type_atom : 
            filout.write ("\t" + str (d_count[neighbor][type_atom]))
        filout.write ("\n")
    
    filout.close ()
    
    return p_result
    



def DistByType (d_dist, subs, pr_result):
    
    
    d_file = {}
    
    l_neighbor = d_dist.keys ()
    l_type_atom = d_dist[l_neighbor[0]].keys ()
    
    print l_neighbor
    print l_type_atom
    
    
    
    for neighbor in l_neighbor :
        if not neighbor in d_file.keys () : 
            d_file[neighbor] = {}
            d_file[neighbor]["density"] = open ( pr_result + str (neighbor) + "_density", "w")
            
            
            for type_atom in l_type_atom : 
                if not type_atom in d_file[neighbor].keys () : 
                    d_file[neighbor][type_atom] = open(pr_result + str (neighbor) + "_" + str (type_atom), "w")
                    if d_dist[neighbor][type_atom] == [] : 
                        continue
                    else : 
                        d_file[neighbor][type_atom].write ("\n".join([str(dist) + "\t" + str (type_atom) + "\t" + str(subs) for dist in d_dist[neighbor][type_atom]]) + "\n")
                        d_file[neighbor]["density"].write ("\n".join([str(dist) + "\t" + str (type_atom) + "\t" + str(subs) for dist in d_dist[neighbor][type_atom]]) + "\n")
    
    
    # close file
    for neighbor in l_neighbor :
        d_file[neighbor]["density"].close ()
        for type_atom in l_type_atom : 
            d_file[neighbor][type_atom].close ()

    return d_file
    
    


def AngleByType(d_angle, pr_angle) : 
    
    l_out = []
    
    for type_angle in d_angle.keys () : 
        print type_angle
        
        p_filout = pr_angle + type_angle
        l_out.append (p_filout)
        filout = open (p_filout, "w")
        
        filout.write ("\n".join([str(d_angle[type_angle]["distance"][i]) + "\t" + str(d_angle[type_angle]["angle"][i]) + "\t" +  str(d_angle[type_angle]["type"][i]) for i in range (0, len(d_angle[type_angle]["type"]))]) + "\n")
        
        filout.close ()
    
    return l_out
    
    
    
    


def countNeighborsAll(stCount, pr_result):
    """
    -> need optimized to pass in argument the folder
    """
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



def CountNeighborRes (st_count, pr_result):
    
    l_res = structure.l_res
    l_sub = structure.ListSub()
    l_sub.append ("global")
    
    d_file_out = {}
    d_file_out["res"] = []
    d_file_out["count"] = []
    
    p_filout_byres = pr_result + "res_count" 
    filout_byres = open (p_filout_byres, "w")
    
    filout_byres.write ("\t".join (l_res) + "\n")
    
    d_file_out["res"].append (p_filout_byres)
    
    for sub in l_sub : 
        p_filout_count = pr_result + "global_count_" + str (sub) 
        
        filout_count = open (p_filout_count, "w")
        
        # global
        filout_count.write ("\n".join ([str (i) for i in st_count[sub]["count"]]) + "\n")
        filout_count.close ()
        # by res
        l_count = []
        for res in l_res : 
            l_count.append (str (st_count[sub]["res"][res]))
        filout_byres.write (str (sub) + "\t" + "\t".join (l_count) + "\n")
        
        d_file_out["count"].append (p_filout_count)
    
    filout_byres.close ()
    
    return d_file_out
        
        
        
        
        
        
        
        
        
        
        
        
         
         
        
        
        
    
    
    pr_result
    
    
    return 
    
    
    


        
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
    


#         filout_side_chain.close()


def resultResProx(stCount, distance_max, pr_result):

    l_p_filout = []
    l_distance = structure.listDistance(distance_max)
    l_aa = ["HOH", "ASP", "GLU", "THR", "SER", "ASN", "GLN", "TYR", "HIS", "LYS", "ARG", "PHE", "TRP", "ALA", "ILE", "LEU", "MET", "VAL", "CYS", "GLY", "PRO"]
    
    for substruct in stCount.keys () : 
        print substruct, "-----"
        print stCount[substruct]
        if stCount[substruct] == {} : 
            continue
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
    




def resultAngle(d_count, pr_out):
    
    l_p_file = []
    for subs in d_count.keys():
        pr_final = pathManage.resultAngle(pr_out, subs)
        p_filout = pr_final + "angle_" + str(subs)
        l_p_file.append (p_filout)
        filout = open(p_filout, "w")
        for type_atom in  d_count[subs].keys():
            nb_dist = len(d_count[subs][type_atom]["distance"])
            for i in range(0, nb_dist) : 
                dist = d_count[subs][type_atom]["distance"][i]
                #filout.write("%.2f" % dist)
                
                for angle in d_count[subs][type_atom]["angles"][i] : 
                    if angle == "NA" :
                        continue
                    else : 
                        angle = "%.2f" % angle
                        filout.write (str (dist) + "\t" + str (angle) + "\t" + str (type_atom) + "\n")
                        
        filout.close()
        
    return l_p_file
 

def dAngleType (count, directory_in):
    
    listClasse = structure.classificationATOM("", out_list = 1)
    for type_substruct in count.keys() : 
        for distance in count[type_substruct].keys() : 
            pr_angle_type = pathManage.resultAngle(directory_in, type_substruct)
            filout = open(pr_angle_type + "angle_" + type_substruct + "_" + distance, "w")
            
            ### HEADER ###
            filout.write (str(listClasse[0]))
            for classe in listClasse[1:] : 
                filout.write ("\t" + classe)
            filout.write("\n")
            ##############
            
            for angleSubs in count[type_substruct][distance].keys() : 
                filout.write(str(angleSubs))
                for classe in listClasse : 
                    try : filout.write("\t" + str(count[type_substruct][distance][angleSubs][classe]))
                    except : filout.write("\t0")
                filout.write("\n")
        filout.close()    


def openFilesWithoutSummary(distanceMax, directory_in):
    
    fileClass = {}
    
    listDistance = structure.listDistance(distanceMax)
    listStructureStudy = structure.ListSub()
    
    for distance in listDistance : 
        fileClass[distance] = {}
        for struct in listStructureStudy :
            dir_struct = pathManage.withoutAtLeastOneSummary(directory_in)
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
        
        
        
def coordinates3DPDBbyNeighborType (l_atom_in, ll_atom_subs, subs, pr_init) : 
    
    d_file = {}
    l_type_neighbors = structure.classificationATOM (out_list = 1)
    l_atom_sub_ref =  structure.substructureCoord(subs)
    for type_neighbors in l_type_neighbors : 
        d_file[type_neighbors] = open (pr_init + subs + "_" + type_neighbors + ".pdb", "w")
        writePDBfile.coordinateSection(d_file [type_neighbors], l_atom_sub_ref , "HETATM") 
        
        #writePDBfile.coordinateSection(d_file [type_neighbors], ll_atom_subs[0] , "HETATM")
        
        #l_write = []
        #for l_atom_subs in ll_atom_subs : 
        #    writePDBfile.coordinateSection(d_file [type_neighbors], l_atom_subs , "HETATM")
#             l_write = l_write + l_atom_subs
#         writePDBfile.coordinateSection(d_file [type_neighbors], l_write , "HETATM")
    
    
    for atom in l_atom_in : 
        type_neighbor = structure.classificationATOM (atom)
        writePDBfile.coordinateSection(d_file [type_neighbor], [atom] , "ATOM", header = 0)
        
    for type_neighbors in l_type_neighbors : 
        d_file[type_neighbors].close ()



def coordinates3DPDB (l_atom_in, l_atom_subs = [], subs = "", p_filout = "") : 
    
    if p_filout == "" : 
        return
    
    filout = open (p_filout, "w")
    
    l_atom_ref =  structure.substructureCoord(subs)
    writePDBfile.coordinateSection(filout, l_atom_ref , "HETATM")
    writePDBfile.coordinateSection(filout, l_atom_subs , "HETATM")
    writePDBfile.coordinateSection(filout, l_atom_in, "ATOM")
    
    filout.close ()



def RelationAngleDistNeighbors (d_relation_neighbors, pr_result) : 
    
    l_filout = []
    for subs in d_relation_neighbors.keys () : 
        if subs == "global" : continue
        
        pr_sub = pr_result + subs + "/"
        pathManage.CreatePathDir(pr_sub)
        
        filout_1_2 = open (pr_sub + subs + "_angles1_2VSDist", "w")
        filout_1_3 = open (pr_sub + subs + "_angles1_3VSDist", "w")
        filout_2_3 = open (pr_sub + subs + "_angles2_3VSDist", "w")
        
        l_filout.append (pr_sub + subs + "_angles1_2VSDist")
        l_filout.append (pr_sub + subs + "_angles1_3VSDist")
        l_filout.append (pr_sub + subs + "_angles2_3VSDist")
        
        i = 0
        try : nb_neighbor = min ([len(d_relation_neighbors[subs][2]["classe"]), len(d_relation_neighbors[subs][1]["classe"]), len(d_relation_neighbors[subs][3]["classe"])])
        except : nb_neighbor = 0
        while i < nb_neighbor :   
            
            print "distance", len (d_relation_neighbors[subs]["distance1_2"]), len (d_relation_neighbors[subs]["distance1_3"]), len (d_relation_neighbors[subs]["distance2_3"])
            print "angleSubs", len (d_relation_neighbors[subs]["angle1_2"]), len (d_relation_neighbors[subs]["angle1_3"]), len (d_relation_neighbors[subs]["angle2_3"])
            print "type", len(d_relation_neighbors[subs][2]["classe"]), len(d_relation_neighbors[subs][1]["classe"]), len(d_relation_neighbors[subs][3]["classe"])
            
            filout_1_2.write (str(d_relation_neighbors[subs]["distance1_2"][i]) + "\t" + str(d_relation_neighbors[subs]["angle1_2"][i]) + "\t" + str(d_relation_neighbors[subs][1]["classe"][i]) + "_" + str(d_relation_neighbors[subs][2]["classe"][i]) + "\n")
            filout_2_3.write (str(d_relation_neighbors[subs]["distance2_3"][i]) + "\t" + str(d_relation_neighbors[subs]["angle2_3"][i]) + "\t" + str(d_relation_neighbors[subs][2]["classe"][i]) + "_" + str(d_relation_neighbors[subs][3]["classe"][i]) + "\n")
            filout_1_3.write (str(d_relation_neighbors[subs]["distance1_3"][i]) + "\t" + str(d_relation_neighbors[subs]["angle1_3"][i]) + "\t" + str(d_relation_neighbors[subs][1]["classe"][i]) + "_" + str(d_relation_neighbors[subs][3]["classe"][i]) + "\n")
    
            i = i + 1
        
        filout_1_2.close ()
        filout_1_3.close ()
        filout_2_3.close ()
    
    return l_filout
    


