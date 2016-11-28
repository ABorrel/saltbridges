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
from re import search



def main (name_database, max_distance = 5.0, RX = 3.00, RFree = 0.25, option_superimpose = 0, option_on_complexes_by_ligand = 0, option_bond = 0, option_stat = 0, option_stat_dataset = 0, option_merge = 0, verbose = 1):


    #format input
    max_distance = float (max_distance)

    # run one database
    pr_result = pathManage.result(name_database)
    # search ligand in PDB
    searchPDB.ligands(name_database, pr_result)

    # dataset with resolution
    l_p_dataset = datasetFinal.Builder(name_database, RX, RFree, option_on_complexes_by_ligand)


########################
#   Parsing dataset   #
########################
#
    if option_stat_dataset == 1:
        for p_dataset in l_p_dataset:
#             statistic.ParseDataSet(p_dataset)

            # identity sequence
            name_dataset = p_dataset.split ("/")
            name_dataset = name_dataset[-1][0:-4]
            l_PDB = []
            l_dataset = loadFile.resultFilterPDBLigand(p_dataset)
            for lig in l_dataset:
                for PDB in lig["PDB"]:
                    if not PDB in l_PDB:
                        l_PDB.append (PDB)
            datasetFinal.RedondanceAnalysis(l_PDB, name_dataset)

####################
# result directory #
####################


    # run for every dataset -> with diffrent resolution
    # short cut
    l_p_dataset = ["/home/borrel/saltBridgesProject/result/PDB/3.0_0.25_uniquePDB/dataset_3.00.txt", "/home/borrel/saltBridgesProject/result/PDB/3.0_0.25_uniquePDB/dataset_1.50.txt"]
#    l_p_dataset = ["/home/borrel/saltBridgesProject/result/PDB/3.0_0.25_uniquePDB/dataset_3.00.txt"] #, "/home/borrel/saltBridgesProject/result/PDB/3.0_0.25_uniquePDB/dataset_1.50.txt"]

    for p_dataset in l_p_dataset:
        pr_result = pathManage.CreatePathDir(p_dataset[:-4] + "/")
        pr_hetion = pathManage.CreatePathDir(p_dataset[:-4] + "/HET/")

        if verbose == 1:
            print "== control path Main =="
            print pr_result
            print pr_hetion
            print "======================="


#       # stat -> build structure, not filter is !!!
        d_sub_neighbor = searchPDB.globalSearch(max_distance, p_dataset, pr_result)

        # remove iron close -> statistic before
        # Becarful because the dictionnary change
        print "control-1", len(d_sub_neighbor["I"])
        d_close_het = hetCloseAnalysis.removeNeighborIron (d_sub_neighbor, pr_hetion + "ionSummarySubstruct.txt")
        print "control-2", len(d_sub_neighbor["I"])

        if option_superimpose == 1:
            # superimpose neighbors -> refaire a Helsinki car MAJ de de la PDB
            superimpose.globalNeighbor (d_sub_neighbor, "I", pr_result)
            superimpose.globalNeighbor (d_sub_neighbor, "II", pr_result)
            superimpose.globalNeighbor (d_sub_neighbor, "III", pr_result)
            superimpose.globalNeighbor (d_sub_neighbor, "IMD", pr_result)
            superimpose.globalNeighbor (d_sub_neighbor, "GAI", pr_result)
            superimpose.globalNeighbor (d_sub_neighbor, "COO", pr_result)

#             superimpose neighbors -> with het first stabilization
#             superimpose.globalNeighbor (d_close_het, "I", pr_hetion)
#             superimpose.globalNeighbor (d_close_het, "II", pr_hetion)
#             superimpose.globalNeighbor (d_close_het, "III", pr_hetion)
#             superimpose.globalNeighbor (d_close_het, "IMD", pr_hetion)
#             superimpose.globalNeighbor (d_close_het, "GAI", pr_hetion)

        if option_bond == 1:
            # check planarity imidazole + guanidium
            statistic.planarityImidazole (d_sub_neighbor, pr_result)
            statistic.planarityGuanidium (d_sub_neighbor, pr_result)

            statistic.lenBondAnalysis(d_sub_neighbor, "I", pr_result)
            statistic.lenBondAnalysis(d_sub_neighbor, "II", pr_result)
            statistic.lenBondAnalysis(d_sub_neighbor, "III", pr_result)

        if option_stat == 1:
            # statistic
            statistic.globalRunStatistic(d_sub_neighbor, max_distance, pr_result)
#             statistic.globalRunStatistic(d_close_het, max_distance, pr_hetion)

    if option_merge == 1:
        if option_on_complexes_by_ligand == 1:
            statistic.MergeDataSet(pathManage.result(name_database + "/" + str (RX) + "_" + str (RFree) + "_uniquePDB"), "dataset_1.50", "dataset_3.00")
            statistic.MergeDataSet(pathManage.result(name_database + "/" + str (RX) + "_" + str (RFree) + "_uniquePDB"), "dataset_1.50", "dataset_3.00", arom = 1)

        else:
            statistic.MergeDataSet(pathManage.result (name_database + "/" + str (RX) + "_" + str (RFree)), "dataset_1.50", "dataset_3.00")




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



def ProteinStat (name_database, option_on_complexes_by_ligand, RX, RFree, max_distance, nb_cases = 150000, nb_run = 1):


    l_p_dataset = datasetFinal.Builder(name_database, RX, RFree, option_on_complexes_by_ligand)
    print l_p_dataset

    for p_dataset in l_p_dataset:
        if search (str (RX), p_dataset.split ("/")[-1]) :
            l_PDB = []
            l_dataset = loadFile.resultFilterPDBLigand(p_dataset)
            for lig in l_dataset:
                for PDB in lig["PDB"]:
                    if not PDB in l_PDB:
                        l_PDB.append (PDB)
            print len (l_PDB), "NB PDB l-178 mainSaltbridge"
            pr_out = pathManage.result("ProtStat" + str (name_database) + str (RX) + "-" + str (RFree))

    # statistic -> identic to protein-ligand
    for i in range (nb_run):
        d_protein = searchPDB.SearchEnvironmentSaltBridgeProt(pr_out, l_PDB, max_distance, nb_cases, 1)
        statistic.globalRunStatistic(d_protein, max_distance, pathManage.CreatePathDir(pr_out + str (nb_cases) + "_" + str (i + 1) + "/"))

    # control different based on countAllpercent in neighbor area
    #l_pcompare3 = []
    #l_pcompare4 = []
    #for i in range (nb_run):
    #    pr_run = pathManage.CreatePathDir(pr_out + str (nb_cases) + "_" + str (i + 1) + "/")
    #    p_file_csv3 = pathManage.twoArea(pr_run, "neighborArea1_3.0") + "countAllpercent.csv"
    #    l_pcompare3.append (p_file_csv3)
    #    p_file_csv4 = pathManage.twoArea(pr_run, "neighborArea1_4.0") + "countAllpercent.csv"
    #    l_pcompare4.append (p_file_csv4)

    #statistic.CompareMultiRun (l_pcompare3, pr_out + "compare_d3.txt")
    #statistic.CompareMultiRun (l_pcompare4, pr_out + "compare_d4.txt")

#     l_p_file_byaa = statistic.SaltBridgeProt (l_PDB, pr_out, thresold_interact, thresold_max)
#     for p_file_byaa in l_p_file_byaa:
#         print p_file_byaa
#         runScriptR.ProtAnalysis (p_file_byaa)



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
max_distance = 6.0
# option_on_complexes_by_ligand = 0
# option_angle = 1

RX_thresold = 3.0
RFree_thresold = 0.25


#RUN all
#PDB50
# main ("PDB50", max_distance = max_distance, option_on_complexes_by_ligand = 1, RX = RX_thresold, RFree = RFree_thresold, option_superimpose = 0, option_bond = 0, option_stat = 1, option_stat_dataset = 1)
# main ("PDB50", max_distance = max_distance, option_on_complexes_by_ligand = 1, RX = RX_thresold, RFree = RFree_thresold, option_superimpose = 0, option_bond = 0, option_stat = 1, option_stat_dataset = 1)

# # PDB
main ("PDB", max_distance = max_distance, option_on_complexes_by_ligand = 1, RX = RX_thresold, RFree = RFree_thresold, option_superimpose = 0, option_bond = 0,  option_stat = 1, option_stat_dataset = 0, option_merge = 1)
# main ( "PDB", max_distance = max_distance, option_on_complexes_by_ligand = 0, RX = RX_thresold, RFree = RFree_thresold, option_superimpose = 0, option_bond = 0,  option_stat = 0, option_stat_dataset = 0, option_merge = 0)

# test
# main ("test", max_distance = max_distance, option_on_complexes_by_ligand = 1, RX = RX_thresold, RFree = RFree_thresold, option_superimpose = 1, option_bond = 0, option_stat = 0, option_stat_dataset = 0)


###################################
#    Stat salt bridge protein     #
###################################


#ProteinStat("PDB", option_on_complexes_by_ligand = 1, RX = 1.5, RFree = RFree_thresold, max_distance = max_distance, nb_cases = 20000, nb_run= 1)
#ProteinStat("PDB", option_on_complexes_by_ligand = 1, RX = 3.0, RFree = RFree_thresold, max_distance = max_distance, nb_cases = 20000, nb_run= 1)

# ProteinStat("PDB50", option_on_complexes_by_ligand = 1, RX = 1.5, RFree = RFree_thresold, max_distance = max_distance, nb_cases = 20000, nb_run= 5)
# ProteinStat("PDB50", option_on_complexes_by_ligand = 1, RX = 3.0, RFree = RFree_thresold, max_distance = max_distance, nb_cases = 20000, nb_run= 5)


##############################
#       Volume function      # -> rewrite with new criterion OK
##############################
# pr_result = pathManage.result()
# volumeFonction.AreaPrimary(pr_result)
# volumeFonction.AreaSecondary(pr_result)
# volumeFonction.AeraTertiary(pr_result)
# volumeFonction.AreaImidazole(pr_result)
# volumeFonction.AreaGuanidium(pr_result)
# volumeFonction.AreaCOO(pr_result)

#######################################
#      Bond length -> criteria        #
#######################################

# bondLength.GlobalBondLength("PDB", 1.5)
# bondLength.GlobalBondLength("PDB", 3.0)


############################
#     Water analysis       # -> +++ to look up the solvant exposition
############################

# waterGlobal ("PDB20", limit_acc = 20.0)
# waterGlobal ("PDB50", limit_acc = 20.0)
# waterGlobal ("PDB50", limit_acc = 0.0)
# waterGlobal ("PDB", limit_acc = 0.0)
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

