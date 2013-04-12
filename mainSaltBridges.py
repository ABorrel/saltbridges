import searchPDB
import datasetFinal
import tool
import test
import statistic
import runScriptR
import volumeFonction
import managePDB
import repertory
import os




def main (path_folder_database, name_folder_result, max_distance = 5.0, option_on_complexes_by_ligand = 0, option_angle = 0, distanceAtoms=3.5,distanceResidues= 5.0):
    
    
    #format input
    max_distance = float (max_distance)
    name_database = path_folder_database.split("/")[-2]
    
    # run one database
    path_dir_result_global = repertory.result (name_folder_result)
    searchPDB.ligands(path_folder_database, path_dir_result_global)
    list_path_file_dataset = datasetFinal.construction(name_folder_result)
    

    
    ########################
    #   Parsing dataset   #
    ########################
    
    for path_dataSet in list_path_file_dataset : 
        statistic.parseDataSet(path_dataSet)
    
    ####################
    # result directory #
    ####################
    
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
        
        
        print "########"
        print path_dir_result
        print "#########"
        # stat -> build structure
        statistic.neighborsAmine(max_distance, path_file_dataset, option_on_complexes_by_ligand, option_angle, path_dir_result)
        
        # draw graph
        runScriptR.globalStat(distanceAtoms, distanceResidues,path_dir_result)





##############################
#           MAIN             #
##############################


###########################
#   Manage PDB file       #
###########################

# extract file and uncompress
#managePDB.formatFilePDB()

# managePDB.retrievePDB( "/home/borrel/saltBridgesProject/PDB/", "/home/borrel/saltBridgesProject/PDB20.dat")
# managePDB.retrievePDB( "/home/borrel/saltBridgesProject/PDB/", "/home/borrel/saltBridgesProject/PDB50.dat")



#####################################
#   Dataset building PDB file       #
#####################################
#Parameters
# path_folder_database = "/home/borrel/saltBridgesProject/PDB20/"
path_folder_database = "/home/borrel/saltBridgesProject/PDB50/"
name_folder_result = "PDB50"
max_distance = 5.0
option_on_complexes_by_ligand = 1
option_angle = 1
distanceAtoms= 3.5
distanceResidues= 5.0



main (path_folder_database, name_folder_result, max_distance = max_distance, option_on_complexes_by_ligand = option_on_complexes_by_ligand, option_angle = option_angle, distanceAtoms=distanceAtoms,distanceResidues= distanceResidues)



##############################
#       Volume function      #
##############################

# volumeFonction.primaryAmine("Primary.pdb", 90, 150, "primary")
# volumeFonction.secondaryAmine("Secondary.pdb", 90, 150, "secondary")
# volumeFonction.tertiaryAmine("Tertiary.pdb", 90, 150, "tertiary")
# volumeFonction.imidazole("IMD.pdb", 1, 30, "imidazole")
# volumeFonction.guanidium("Guanidium.pdb", 90, 150,1,30, "guanidium")
# volumeFonction.pyridine("Pyridine.pdb", 1, 30, "pyridine")
# volumeFonction.diamine("Diamine.pdb", 90, 150, "diamine")


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




