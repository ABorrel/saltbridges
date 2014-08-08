import searchPDB
import datasetFinal
import tool
import hetCloseAnalysis
import test
import statistic
import runScriptR
import volumeFonction
import managePDB
import repertory
import os
import runOtherSoft
from time import sleep
import waterAnalysis
import superimpose



import loadFile




def main (name_database, max_distance = 5.0, option_on_complexes_by_ligand = 0, option_angle = 0, RX = 3.0, RFree = 0.25, verbose = 1):
    
    
    #format input
#     max_distance = float (max_distance)
     
    # run one database
    path_dir_result_global = repertory.result (name_database)
    
    # search ligand in PDB
    searchPDB.ligands(name_database, path_dir_result_global)
    
    # dataset with resolution
    list_path_file_dataset = datasetFinal.construction(name_database,RX, RFree )
#     
#
#     ########################
#     #   Parsing dataset   #
#     ########################
#     
#     for path_dataSet in list_path_file_dataset : 
#         statistic.parseDataSet(path_dataSet, 1)
#         statistic.parseDataSet(path_dataSet, 0)
     
#     ####################
#     # result directory #
#     ####################
#     
#
    # run for every dataset -> with diffrent resolution
    # short cut
    list_path_file_dataset = ["/home/borrel/saltBridgesProject/result/PDBTest/dataset_3.00"]
    print list_path_file_dataset

    
# #     
    for path_file_dataset in list_path_file_dataset : 
        
        name_folder =  path_file_dataset.split("_")[-1]
        if option_angle == 1 : 
            name_folder = name_folder + "_angle"
        else : 
            name_folder = name_folder + "_noangle"
                
        if option_on_complexes_by_ligand == 1 : 
            name_folder = name_folder + "_onecomplexe"
        else : 
            name_folder = name_folder + "_morecomplexe"
            
        path_dir_result = repertory.result (name_database + "/" + name_folder)
        p_dir_result_het = repertory.result (name_database + "/" + name_folder + "het")
            
           
        print "########"
        print path_dir_result
        print p_dir_result_het
        print "#########"
        
        
#         # stat -> build structure, not filter is !!!
    atom_interest_close, global_atom_close = searchPDB.globalSearch(max_distance, path_file_dataset, path_dir_result_global)
        
        # remove iron close -> statistic before 
        # Becarful because the dictionnary change
#     atom_interest_het = hetCloseAnalysis.removeNeighborIron (atom_interest_close)
#     global_atom_het = hetCloseAnalysis.removeNeighborIron (global_atom_close)
        
#         superimpose neighbors -> refaire a Helsinki car MAJ de de la PDB
#         superimpose.globalNeighbor (atom_interest_close, "Primary", path_dir_result)
#         superimpose.globalNeighbor (atom_interest_close, "Secondary", path_dir_result)
#         superimpose.globalNeighbor (atom_interest_close, "Tertiary", path_dir_result)
#         superimpose.globalNeighbor (atom_interest_close, "Imidazole", path_dir_result)
#         superimpose.globalNeighbor (atom_interest_close, "Guanidium", path_dir_result)
#         superimpose.globalNeighbor (atom_interest_close, "AcidCarboxylic", path_dir_result)
        
# #         superimpose neighbors -> with het first stabilization 
#         superimpose.globalNeighbor (atom_interest_het, "Primary", p_dir_result_het)
#         superimpose.globalNeighbor (atom_interest_het, "Secondary", p_dir_result_het)
#         superimpose.globalNeighbor (atom_interest_het, "Tertiary", p_dir_result_het)
#         superimpose.globalNeighbor (atom_interest_het, "Imidazole", p_dir_result_het)
#         superimpose.globalNeighbor (atom_interest_het, "Guanidium", p_dir_result_het)
        
        # check planarity imidazole + guanidium
# # # #         statistic.planarityImidazole (atom_interest_close, path_dir_result)
# # # #         statistic.planarityGuanidium (atom_interest_close, path_dir_result)
             
        # analyse length bond not use correctly because limited by crystallo quality !!!!
# # #         statistic.lenBondAnalysis(atom_interest_close, "Primary",path_dir_result)
# # #         statistic.lenBondAnalysis(atom_interest_close, "Secondary",path_dir_result)
# # #         statistic.lenBondAnalysis(atom_interest_close, "Tertiary",path_dir_result)


        # statistic
    statistic.globalRunStatistic(atom_interest_close, global_atom_close, max_distance, option_angle, path_dir_result)
#     statistic.globalRunStatistic(atom_interest_het, global_atom_het, max_distance, option_angle, p_dir_result_het)
    
    # draw plot
#     runScriptR.globalStat(max_distance, path_dir_result)
#     runScriptR.globalStat(max_distance,  p_dir_result_het)
    
    
        


def waterGlobal (name_database, limit_acc = 20.0):
    """
    Number of water molecules in PDB
    arg: -> Path folder database
         -> name folder result
         -> limit acc
    return: NONE
    """
    
    path_dir_result_global = repertory.result (name_database)
    
    # retrieve list PDB file
    list_PDBID = managePDB.retriveListPDB(name_database)
    print len (list_PDBID)
    
    # calcul acc with NACESS
    for PDB_ID in list_PDBID :
        path_file_PDB = repertory.pathDitrectoryPDB () + PDB_ID + ".pdb"
        runOtherSoft.runNACESS(path_file_PDB)
    
    sleep(10)
    try :
        os.system("mv *.asa " +  repertory.pathDitrectoryPDB ())
        os.system("mv *.rsa " +  repertory.pathDitrectoryPDB ())
        os.system("rm *.log")
    except : 
        pass
    
    path_file_result = waterAnalysis.resolutionWater(list_PDBID, path_dir_result_global, limit_acc)

    runScriptR.waterPlotResolution (path_file_result)
    


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
# path_folder_database = "/home/borrel/saltBridgesProject/PDB20/"
# path_folder_database = "/home/borrel/saltBridgesProject/PDB50/"
# name_folder_result = "PDB50"
max_distance = 4.0
# option_on_complexes_by_ligand = 0
# option_angle = 1

RX_thresold = 3.0
RFree_thresold = 0.25


#RUN all
#PDB 50
# main ("PDB50", max_distance = max_distance, option_on_complexes_by_ligand = 0, option_angle = 0, RX = RX_thresold, RFree = RFree_thresold )
# main ("PDB50", max_distance = max_distance, option_on_complexes_by_ligand = 1, option_angle = 0)
# 
# # PDB
# main ( "PDB", max_distance = max_distance, option_on_complexes_by_ligand = 0, option_angle = 0)

# test
main ("PDBTest", max_distance = max_distance, option_on_complexes_by_ligand = 0, option_angle = 0, RX = RX_thresold, RFree = RFree_thresold )



##############################
#       Volume function      #
##############################

# volumeFonction.primaryAmine("/home/borrel/saltBridgesProject/result/Primary.pdb", 90, 150, "primary")
# volumeFonction.secondaryAmine("/home/borrel/saltBridgesProject/result/Secondary.pdb", 90, 150, "secondary")
# volumeFonction.tertiaryAmine("/home/borrel/saltBridgesProject/result/Tertiary.pdb", 90, 150, "tertiary")
# volumeFonction.imidazole("/home/borrel/saltBridgesProject/result/Imidazole.pdb", 1, 30, "imidazole")
# volumeFonction.guanidium("/home/borrel/saltBridgesProject/result/Guanidium.pdb", 90, 150,1,30, "guanidium")
# volumeFonction.pyridine("/home/borrel/saltBridgesProject/result/Pyridine.pdb", 1, 30, "pyridine")
# volumeFonction.diamine("/home/borrel/saltBridgesProject/result/Diamine.pdb", 90, 150, "diamine")


############################
#     Water analysis       #
############################

# waterGlobal ("PDB20", limit_acc = 20.0)
# waterGlobal ("PDB50", limit_acc = 20.0)
# waterGlobal ("PDB", limit_acc = 20.0)
# waterFamily("PDB")
# waterFamily("PDB50")








"""
##############################
#       Help fonctions       # -> found IMD
##############################

import loadFile
import repertory

ligandInPDB = loadFile.resultLigandPDB(repertory.result() + "resultLigandInPDB")

list_ligand = ligandInPDB.keys()

nb_ligand = len(list_ligand)

i = 0
while i < nb_ligand : 
    if list_ligand[i] == "IMD" : 
        
        print i
        i = nb_ligand
    else :
        i = i + 1


for dataset in listDataSet : 
    for onlyOnePDBByLigand in range(0,2) :
        for optionAngle in range(0,1) :  
            statistic.neighborsAmine(distanceMax, dataset,onlyOnePDBByLigand, optionAngle)
            runScriptR.globalStat(3.5, distanceMax)
            datasetRep = dataset
            if onlyOnePDBByLigand == 1 : 
                datasetRep = datasetRep + "_only1PDBbyligand_" + str(distanceMax)
            if optionAngle == 1 : 
                datasetRep = datasetRep + "_SelectAngles"
            
            tool.moveResult(datasetRep)

"""
