import searchPDB
import datasetFinal
import tool
import test
import statistic
import runScriptR
import volumeFonction
import managePDB
import repertory






def constructDataSet (path_folder_PDB, name_folder_result):
    """
    Search ligand in dataset and select with different resolution filters
    arg: - path folder database
         - name folder result
         - parsing dataset, simple analysis
    return: NONE
    """
    
    path_dir_result = repertory.result (name_folder_result)
    searchPDB.ligands(path_folder_PDB, path_dir_result)
    list_path_file_dataset = datasetFinal.construction(name_folder_result)
    print list_path_file_dataset, "check"
    

    
    ########################
    #   Parsing dataset   #
    ########################
    
    for path_dataSet in list_path_file_dataset : 
        statistic.parseDataSet(path_dataSet)
    
    return list_path_file_dataset
        

def statisticWithGroup (path_file_dataset, max_distance = 5.0, option_on_complexes_by_ligand = 0, option_angle = 0):
    
    # format 
    max_distance = float (max_distance)
    
    # stat
    statistic.neighborsAmine(max_distance, path_file_dataset, option_on_complexes_by_ligand, option_angle)



        


##############################
#     statistic neighbors    #
##############################

#########################################################################################
# dataset = listDataSet[0]
# onlyOnePDBByLigand = 0
#
#print tool.searchLigandInDataSetFile(dataset, "IMD")
#
# optionAngle = 1
# statistic.neighborsAmine(5.0, dataset, onlyOnePDBByLigand, optionAngle)
# runScriptR.globalStat(3.5, 5.0)

 ##if onlyOnePDBByLigand == 1 : 
# ##dataset = dataset + "_only1PDBbyligand"
# tool.moveResult("test")


##########################################################################################
"""distanceMax = 5.0
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





##############################
#       Volume fonction      #
##############################

# volumeFonction.primaryAmine("Primary.pdb", 90, 150, "primary")
# volumeFonction.secondaryAmine("Secondary.pdb", 90, 150, "secondary")
# volumeFonction.tertiaryAmine("Tertiary.pdb", 90, 150, "tertiary")
# volumeFonction.imidazole("IMD.pdb", 1, 30, "imidazole")
# volumeFonction.guanidium("Guanidium.pdb", 90, 150,1,30, "guanidium")
# volumeFonction.pyridine("Pyridine.pdb", 1, 30, "pyridine")
# volumeFonction.diamine("Diamine.pdb", 90, 150, "diamine")










##############################
#           MAIN             #
##############################
#Parameters
path_folder_PDB = "/home/borrel/saltBridgesProject/PDBTest/"
name_folder_result = "PDBTest"


###########################
#   Manage PDB file       #
###########################

# extract file and uncompress
#managePDB.formatFilePDB()


#####################################
#   Dataset building PDB file       #
#####################################
list_path_dataset = constructDataSet (path_folder_PDB, name_folder_result)

print list_path_dataset

# statisticWithGroup (path_file_dataset, max_distance = 5.0, option_on_complexes_by_ligand = 0, option_angle = 0)







"""

##############################
#       Help fonctions       # -> found IMZ
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

"""




