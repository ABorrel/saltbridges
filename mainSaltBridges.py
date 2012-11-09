#! /usr/bin/env python

import searchPDB
import datasetFinal
import tool
import test
import statistic
import runScriptR
import volumeFonction
import managePDB

"""Typage code : 
-> programme procedural, separation des modules en fonction de leur utilite, calcul pour faire les calculs ....
-> Je ne rappel pas dans le nom des fonctions le modules par exemple je ne mets pas calculDistancetoAtom mais distanceTowAtom
-> Pour le nom de fonction je n utilise pas de caracteres speciaux type _ mais je met des MAJ pour separer les noms
-> J essais d avoir des noms de varibles explicites et de preciser le typage a la sortie des fonctions par exemple listAtom
-> Je commente tres peu :) mais je mets toujours au debut de fonction le in et out et descriptif rapide

NB : 1) Si tu veux lancer le script ailleur modifie bien les repertoires et garde la meme architecture de dossier et si tu veux faire tourner 
le script sur une autre database tu dois le pouvoir en changant juste le repertoire de la PDB (je pense nottament sur le projet de Leslie)
     2) Eclipse plante un peu chez Youe, je n ai pas eu le temps de le debuger mais pour faire tourner mon script sur eclipse je te conseil de faire
un nouveau projet et de ensuite copie coller mes scripts dans le dossier et ca devrait passer sans trop de pb. Et pour specifier linterpreteur,
jai explique a Youe mais pas sur quil est compris, il installe python 2.7 dans /usr/local tu as juste a preciser ce repertoire
 
 Bonne continuation
"""
    


listDataSet = ["dataset_3.00", "dataset_2.00"]

###########################
#   Manage PDB file       #
###########################

#managePDB.formatFilePDB()


############################
#    dataset construction  #
############################


#searchPDB.ligands()
datasetFinal.construction()

########################
#   Parsing dataset   #
########################

#for dataSet in listDataSet : 
#    statistic.parseDataSet(dataSet)

##############################
#     statistic neighbors    #
##############################

#########################################################################################
#dataset = listDataSet[1]
#onlyOnePDBByLigand = 0
#
#print tool.searchLigandInDataSetFile(dataset, "IMD")
#
#optionAngle = 1
#statistic.neighborsAmine(5.0, dataset,onlyOnePDBByLigand, optionAngle)
#runScriptR.globalStat(3.5, 5.0)
###if onlyOnePDBByLigand == 1 : 
###dataset = dataset + "_only1PDBbyligand"
#tool.moveResult("test")


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

#volumeFonction.primaryAmine("Primary.pdb", 90, 150, "primary")
#volumeFonction.secondaryAmine("Secondary.pdb", 90, 150, "secondary")
#volumeFonction.tertiaryAmine("Tertiary.pdb", 90, 150, "tertiary")
#volumeFonction.imidazole("IMD.pdb", 1, 30, "imidazole")
#volumeFonction.guanidium("Guanidium.pdb", 90, 150,1,30, "guanidium")
#volumeFonction.pyridine("Pyridine.pdb", 1, 30, "pyridine")
#volumeFonction.diamine("Diamine.pdb", 90, 150, "diamine")








##############################
#       Help fonctions       #
##############################

import loadFile
import repertory

ligandInPDB = loadFile.resultLigandPDB(repertory.result() + "resultLigandInPDB")

list_ligand = ligandInPDB.keys()

nb_ligand = len(list_ligand)

i = 0
while i< nb_ligand : 
    if list_ligand[i] == "IMD" : 
        
        print i
        i = nb_ligand
    else :
        i = i+1







