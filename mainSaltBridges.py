# personnel
import searchPDB
import datasetFinal
import tool
import hetCloseAnalysis
import test
import statistic
import runScriptR
import volumeFonction
import managePDB
import pathManage
import runOtherSoft
import waterAnalysis
import superimpose
import writeICMScript
import loadFile
import GPCRDockAnalysis
import bondLength


# global
import os
from time import sleep



def main (name_database, max_distance = 5.0, RX = 3.0, RFree = 0.25, option_superimpose = 0, option_on_complexes_by_ligand = 0, option_bond = 0, option_stat = 0, option_stat_dataset = 0,  verbose = 1):
    
    
    #format input
    max_distance = float (max_distance)
     
    # run one database
    pr_result = pathManage.result (name_database)
    
    # search ligand in PDB
    searchPDB.ligands(name_database, pr_result)
    
    # dataset with resolution
    l_p_dataset = datasetFinal.Builder(name_database, RX, RFree, option_on_complexes_by_ligand)
#
#     ########################
#     #   Parsing dataset   #
#     ########################
# # 
    if option_stat_dataset == 1 : 
        for p_dataset in l_p_dataset : 
            statistic.ParseDataSet(p_dataset)
        
     
#     ####################
#     # result directory #
#     ####################
#     
#
    # run for every dataset -> with diffrent resolution
    # short cut
#     l_p_dataset = ["/home/borrel/saltBridgesProject/result/PDB/3.0_0.25_uniquePDB/dataset_3.00.txt" ]
# #     
    for p_dataset in l_p_dataset : 
        
        pr_result = pathManage.CreatePathDir(p_dataset[:-4] + "/")
        pr_hetion = pathManage.CreatePathDir(p_dataset[:-4] + "/HET/")
        
        if verbose == 1 :  
            print "== control path Main =="
            print pr_result
            print pr_hetion
            print "======================="
        
        
#         # stat -> build structure, not filter is !!!
        d_sub_neighbor = searchPDB.globalSearch(max_distance, p_dataset, pr_result)
    
        # remove iron close -> statistic before 
        # Becarful because the dictionnary change
        print "control-1", len(d_sub_neighbor["Primary"])
        d_close_het = hetCloseAnalysis.removeNeighborIron (d_sub_neighbor, pr_hetion + "ionSummarySubstruct.txt")
        print "control-2", len(d_sub_neighbor["Primary"])
        
        if option_superimpose == 1 : 
            # superimpose neighbors -> refaire a Helsinki car MAJ de de la PDB
            superimpose.globalNeighbor (d_sub_neighbor, "Primary", pr_result)
            superimpose.globalNeighbor (d_sub_neighbor, "Secondary", pr_result)
            superimpose.globalNeighbor (d_sub_neighbor, "Tertiary", pr_result)
            superimpose.globalNeighbor (d_sub_neighbor, "Imidazole", pr_result)
            superimpose.globalNeighbor (d_sub_neighbor, "Guanidium", pr_result)
            superimpose.globalNeighbor (d_sub_neighbor, "AcidCarboxylic", pr_result)
            
            # superimpose neighbors -> with het first stabilization 
#             superimpose.globalNeighbor (d_close_het, "Primary", pr_hetion)
#             superimpose.globalNeighbor (d_close_het, "Secondary", pr_hetion)
#             superimpose.globalNeighbor (d_close_het, "Tertiary", pr_hetion)
#             superimpose.globalNeighbor (d_close_het, "Imidazole", pr_hetion)
#             superimpose.globalNeighbor (d_close_het, "Guanidium", pr_hetion)
        
        if option_bond == 1 : 
        
            # check planarity imidazole + guanidium
            statistic.planarityImidazole (d_sub_neighbor, pr_result)
            statistic.planarityGuanidium (d_sub_neighbor, pr_result)
            
            statistic.lenBondAnalysis(d_sub_neighbor, "Primary", pr_result)
            statistic.lenBondAnalysis(d_sub_neighbor, "Secondary", pr_result)
            statistic.lenBondAnalysis(d_sub_neighbor, "Tertiary", pr_result)

        if option_stat == 1: 
            # statistic
            statistic.globalRunStatistic(d_sub_neighbor, max_distance, pr_result)
#             statistic.globalRunStatistic(d_close_het, max_distance, pr_hetion)
    
 
 


def waterGlobal (name_database, limit_acc = 00.0):
    """
    Number of water molecules in PDB
    arg: -> Path folder database
         -> name folder result
         -> limit acc
    return: NONE
    """
    
    pr_result = pathManage.result (name_database + "/water")
    
    # retrieve list PDB file
    l_PDBID = managePDB.retriveListPDB(name_database)
    # calcul acc with NACESS
    if limit_acc != 0.0 : 
        for PDB_ID in l_PDBID :
            p_PDB = pathManage.pathDitrectoryPDB () + PDB_ID + ".pdb"
            runOtherSoft.runNACESS(p_PDB, pathManage.pathDitrectoryPDB (), multi_run = 0)
        
    p_filout = waterAnalysis.resolutionWater(l_PDBID, pr_result, limit_acc)

    runScriptR.waterPlotResolution (p_filout)
    


def waterFamily (name_database):

    # water with summary file
    path_files_water = waterAnalysis.resolutionByStructure (name_database)
    for path_file in path_files_water : 
        runScriptR.waterType (path_file)


##############################
#           MAIN             #
##############################


###########################
#   Manage PDB file       #
###########################

# extract file and uncompress
# managePDB.formatFilePDB("/home/borrel/PDB/")


#####################################
#   Dataset building PDB file       #
#####################################
#Parameters
max_distance = 5.0
# option_on_complexes_by_ligand = 0
# option_angle = 1

RX_thresold = 3.0
RFree_thresold = 0.25


#RUN all
# main ("PDB50", max_distance = max_distance, option_on_complexes_by_ligand = 1, RX = RX_thresold, RFree = RFree_thresold, option_superimpose = 0, option_bond = 0, option_stat = 0, option_stat_dataset = 1)

# # PDB
main ( "PDB", max_distance = max_distance, option_on_complexes_by_ligand = 1, RX = RX_thresold, RFree = RFree_thresold, option_superimpose = 0, option_bond = 0,  option_stat = 1, option_stat_dataset = 0)
# main ( "PDB", max_distance = max_distance, option_on_complexes_by_ligand = 0, RX = RX_thresold, RFree = RFree_thresold, option_superimpose = 0, option_bond = 0,  option_stat = 0, option_stat_dataset = 0)

# test
# main ("test", max_distance = max_distance, option_on_complexes_by_ligand = 1, RX = RX_thresold, RFree = RFree_thresold, option_superimpose = 1, option_bond = 0, option_stat = 0, option_stat_dataset = 0)



##############################
#       Volume function      # -> rewrite with new criterion OK
##############################
# pr_result = pathManage.result()
# volumeFonction.AreaPrimary(pr_result)
# volumeFonction.AreaSecondary(pr_result)
# volumeFonction.AeraTertiary(pr_result)
# volumeFonction.AreaImidazole(pr_result)
# volumeFonction.AreaGuanidium(pr_result)


###########################
#      Bond length        # check with new function in the main // write again
###########################

# bondLength.BondCNAndCO("PDB", 1.5)




############################
#     Water analysis       # -> +++ to look up the solvant exposition
############################

# waterGlobal ("PDB20", limit_acc = 20.0)
# waterGlobal ("PDB50", limit_acc = 20.0)
# waterGlobal ("PDB50", limit_acc = 0.0)
# waterGlobal ("PDB", limit_acc = 20.0)
# waterFamily("PDB")
# waterFamily("PDB50")


#################
#   GPCR dock   #
#################

# pr_GPCRDock2010 = "/home/borrel/saltBridgesProject/GPCRDock2010/PDB_conserved/"
# pr_data = "/home/borrel/saltBridgesProject/GPCRDock/data/"
# pr_result = pathManage.result("GPCRDock")

# convertion ICM 
# writeICMScript.ScriptConvertICBtoPDB(pr_GPCRDock2010, pr_result + "convertGPCRDock2010.txt")

#GPCRDockAnalysis.GlobalAnalysisGPCR (pr_data, pr_result, dist_thresold = 5.0, chemical_search = 1, option_stat = 1, option_model = 1)

