'''
Created on 12 mai 2015

@author: borrel
'''

# personnal
import loadFile
import parsing
import searchPDB
import writeFile
import writePDBfile
import statistic

# 
from os import listdir, path
from copy import deepcopy
import pathManage



def GlobalAnalysisGPCR (pr_data, pr_result, dist_thresold = 5.0, chemical_search = 0, option_stat = 0, option_model = 0):

    # search chemical substructures in docked ligand
    if chemical_search == 1 : 
        SearchChemicalSubstruct(pr_data, pr_result, control = 1)
    # build structure
    struct_neighbor = SearchNeighbor (pr_data, pr_result, dist_thresold = dist_thresold, debug = 1)
    
    for drug in struct_neighbor.keys () : 
        pr_drug = pathManage.CreatePathDir(pr_result + str(drug) + "/")
        if option_stat == 1: 
            # statistic -> same statistic as global database 
            statistic.globalRunStatistic(struct_neighbor[drug], dist_thresold, pr_drug)
        
        if option_model == 1 :
            SabilizationBySB(struct_neighbor[drug], pr_drug)
            
    
    



def SabilizationBySB(d_neighbors, pr_result) : 

    for sub in d_neighbors.keys () :
        # global
        if sub == "global" : 
            continue
        pr_resultPS = pathManage.CreatePathDir(pr_result + "StabilizeSB/" + str (sub) + "/")
        print pr_resultPS
        p_filout = pr_resultPS + "summarize"
        filout = open (p_filout, "w")
        filout.write ("Name file\tAtom\tStabilize restrained\tStabilize no restrained\tEnvironment\n") 
        
        print "*******"
        print d_neighbors[sub][1].keys ()
        
        for model in d_neighbors[sub] : 
            
            interaction_restrained = statistic.retrieveInteraction (model["neighbors"], sub, restrained = 1)
            interaction_norestrained = statistic.retrieveInteraction (model["neighbors"], sub, restrained = 0)
            
            environment = statistic.ListNeighborsType(model["neighbors"], sub)


            filout.write (model["PDB"] + "\t" + str(model["resName"]) + "_" + str (model["serial"]) + "\t" + str (interaction_restrained) + "\t" + str (interaction_norestrained) + "\t" + "_".join(environment) + "\n")
        filout.close () 
            




def ControlHETATOMModelFile (pr_data, pr_result, l_lig_control = []):
    
    # Write file ligand  
    
    l_ligand = listdir(pr_data)
    for ligand in l_ligand : 
        # control if file exist and size
        p_file_lig = pathManage.CreatePathDir(pr_result + str (ligand) + "/" ) + "ligInModel"
        
        if path.exists(p_file_lig) : 
            return
        
        else : 
            file_lig = open (p_file_lig, "w")
            print "ligand", ligand
        
            l_file = listdir(pr_data + ligand + "/")
        
            d_temp = {}
            for model_file in l_file :
            
                l_lig_PDB = parsing.retrieveListLigand(pr_data + ligand + "/" + model_file)
                print "***", l_lig_PDB
                if len (l_lig_PDB) == 0 : 
                    continue
                else : 
                    for lig in l_lig_PDB :
                        if not l_lig_control == [] and lig in l_lig_control :  
                            if not lig in d_temp.keys () : 
                                d_temp[lig] = []
                            d_temp[lig].append (pr_data + ligand + "/" + model_file)
                        else : 
                            continue
                
            # write temp control
            for lig in d_temp.keys () : 
                file_lig.write (str (lig) + "\t" + " ".join(d_temp[lig]) + "\n")
    
            file_lig.close ()
        
    return 
        
    
 
 
def SearchNeighbor (pr_data, pr_result, dist_thresold = 5.0, debug = 1): 
    
    # write file with HETATM in PDB -> to used same function
    l_lig_model = ["ZMA", "UNK", "RES", "CAU", "P0G", "CY8", "1KS"] # goal of this list is to remove the modified amino acid
    ControlHETATOMModelFile (pr_data, pr_result, l_lig_control = l_lig_model)

    l_drug = listdir(pr_data)
    d_out = {}
    
    for drug in l_drug :
         
        # pr result
        pr_drug = pathManage.CreatePathDir(pr_result + str(drug) + "/")
        l_lig = loadFile.resultFilterPDBLigand(pr_drug + "ligInModel")
        
        if debug == 1 : 
            print "List ligand for drug->"
            print drug
            print len (l_lig)
            print l_lig
        
        
        if l_lig == [] : 
            print "ERROR: find HET in model -> l68 GPCRDockAnalysis"
            return {}
        else : 
            #file summary
            pr_summary = pr_drug + "Sum/"
            
            # load structure in summary ---> if use need place option one PDB by ligand
            d_neighbor = loadFile.loadCloseStruct (pr_summary,  control_empty_file = 0)
            
            if d_neighbor == None : 
                
                print "Write Sum File"
                #write summary file
                d_files_summary = writeFile.openFileSummary(pr_summary)
            
                # structure of stock
                
                d_neighbor = {}
                l_neighbor_global = []
            
            
                # search neighbor
                nb_lig = len (l_lig)
                i = 0
                
                while i < nb_lig : 
                    j = 0
                    nb_model = len (l_lig[i]["PDB"])
                    print nb_model
                    while j < nb_model : 
                        
                        l_atom_ligand = loadFile.ExtractInfoPDBID(l_lig[i]["PDB"][j])[l_lig[i]["name"]][0] # change not tested
                    
                        # search neighbor for every atom in ligand selected
                        searchPDB.globalNeighbors(dist_thresold, l_atom_ligand, l_lig[i]["PDB"][j], l_neighbor_global)
                        # search neighbor for interest 
                        searchPDB.interestGroup(dist_thresold, l_atom_ligand, l_lig[i]["PDB"][j], d_neighbor, more_flex = 1)
                    
                        j = j + 1
                    i = i + 1
                    
            
                writeFile.neighborStruct(d_neighbor, l_neighbor_global, d_files_summary)
                writeFile.closeFileSummary(d_files_summary)
            
                # case where load directly substructure => why do not load directly in dictionnary
                d_neighbor["global"] = l_neighbor_global
        
        
        #stock
        print "Stock neighbor l164:", d_neighbor.keys()
        if not drug in d_out.keys () :
            d_out[drug] = deepcopy(d_neighbor)
            
        
    return d_out
    
 
 
 
def ControlPDBFormat (p_PDB):
    """
    Rewrite model PDB file 
    """
    
    l_res = ["ILE", "LEU", "LYS", "PHE", "TYR", "VAL", "SER", "MET", "ARG", "TRP", "PRO", "GLY", "GLU", "ASN", "HIS", "ALA", "ASP", "GLN", "THR", "CYS"]
    l_atom_PDB = parsing.loadCoordSectionPDB(p_PDB)
    
    filout = open (p_PDB, "w")
    
    l_atom_het = []
    for atom_PDB in l_atom_PDB : 
        if atom_PDB["resName"] in l_res : 
            writePDBfile.coordinateSection(filout, [atom_PDB], "ATOM", header = 0)
        else :
            l_atom_het.append (deepcopy(atom_PDB)) 
    
    writePDBfile.coordinateSection(filout, l_atom_het, recorder = "HETATM", header = 0, connect_matrix = 1)
    
    filout.close ()
            
    

def SearchChemicalSubstruct (pr_data, pr_result, control = 0):
    
    # substructure search 
    
    p_filout = pr_result + "findSruct"
    
    filout = open (p_filout, "w")
    l_ligand = listdir(pr_data)
    
    print l_ligand
    
    for ligand in l_ligand : 
        print ligand
    
        l_file = listdir(pr_data + ligand + "/")
        
        for f in l_file : 
            if ligand == "ZM241385" or ligand == "source": 
                group = f
            else : 
                group = f.split ("_")[-2]
            p_file_PDB = pr_data + ligand + "/" + f
            print p_file_PDB
            
            if control == 1 : 
                ControlPDBFormat(p_file_PDB)
            
            l_atom = parsing.loadCoordSectionPDB(p_file_PDB, remove_H = 1)
            if ligand == "ZM241385" : 
                ll_atom_lig = parsing.retrieveLigand(l_atom, "ZMA")
            else : 
                ll_atom_lig = parsing.retrieveLigand(l_atom, "RES")
    
            for l_atom_lig in ll_atom_lig : 
                
                l_subs = searchPDB.interestStructure(l_atom_lig, more_flex = 1)
                
                filout.write (str (ligand) + "\t" + str (f) + "\t" + str (group) + "\t" + " ".join(l_subs) + "\n")
    
    filout.close ()
       

# best group ergomine -> 2128 (2013) && 6882 (2013)
# best group taladegid -> 7527 (2013) Vanderblit
# best group eticlopride -> 1285 (2010)
# ZM241385 => mod7msp
    