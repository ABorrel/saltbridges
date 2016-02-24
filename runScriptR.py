from os import system, listdir, path
from re import search

import pathManage
import structure
import log
import tool

def barplotQuantity(distance, type_study, p_filin, logFile, debug = 1):
    """Plot count statistic
    in: distance for legend plot, type histogram (coplar or length), p_filin with data, type of study, logFile
    out: CMD in terminal -> plot """
    
    if tool.checkFileEmpty(p_filin) == 1 : 
        return
    rep = pathManage.scriptR()
    cmd = rep + "barplotQuantity.R " + p_filin + " " + str("%.2f" % distance) + " " + type_study
    if debug : print cmd
    logFile.write(cmd + "\n")
    system(cmd)


def barplotQuantityByAA(distance, aminoAcid, path_file, logFile):
    """to do"""

    if tool.checkFileEmpty(path_file) == 1 : 
        return
    rep = pathManage.scriptR()
    cmd = rep + "barplotQuantityByAA.R " + path_file + " " + str("%.2f" % distance) + " " + aminoAcid + " " + str("%.2f" % distance)
    print cmd
    logFile.write(cmd + "\n")
    system(cmd)


def plotDistance(p_filin, logFile, debug = 1):
    """Plot list distance
    in: distance for legend plot, type histogram (coplar or length), file with data, type of study, logFile
    out: CMD in terminal -> plot """

    repScript = pathManage.scriptR()
    if tool.checkFileEmpty(p_filin) == 1 : 
        return
    cmd = repScript + "plotDistance.R " + p_filin
    if debug : print cmd
    logFile.write(cmd + "\n")
    system(cmd)


def plotDistanceDensity(p_filin, logFile, debug = 1):
    """Plot density"""

    repScript = pathManage.scriptR()
    if tool.checkFileEmpty(p_filin) == 1 : 
        return
    cmd = repScript + "densityDistance.R " + p_filin
    if debug : print cmd
    logFile.write(cmd + "\n")
    system(cmd)





def plotNbNeighbor(p_filin, logFile, debug = 1):
    """Plot number neighbor 
    in: distance for legend plot, type histogram (coplar or length), file with data, type of study, logFile
    out: CMD in terminal -> plot """

    repScript = pathManage.scriptR()
    if tool.checkFileEmpty(p_filin) == 1 : 
        return
    cmd = repScript + "plotNumberNeighbor.R " + p_filin
    if debug : print cmd
    logFile.write(cmd + "\n")
    system(cmd)


def barplotCombination (p_filin, logFile, debug = 1):
    
    repScript = pathManage.scriptR()
    if tool.checkFileEmpty(p_filin) == 1 : 
        return
    cmd = repScript + "barplotCombination.R " + p_filin
    if debug : print cmd
    logFile.write(cmd + "\n")
    system(cmd)




def histProportionType(distanceMax, dir_out, logFile):
    """Draw proportion type graphe
    in: Distance Max study, log file
    out: excecute CMD -> draw plot"""
    
    
    repScript = pathManage.scriptR()
    
    cmd = repScript + "barplotPercentClasse.R " + dir_out
    print cmd
    logFile.write(cmd + "\n")
    system (cmd)

    
def histProportionTypeNeighbors(distance, dir_out, logFile):
    """Draw proportion type neighbor
    in: Distance Max study, log file_data
    out: excecute CMD -> draw plot"""
    
    listStruct = structure.ListSub()
    listStruct.append("GlobalAmine")
    listStruct.append("Global")
    for type in listStruct : 
        file_data = dir_out + "globalProportionType/" + "proportionType" + type + str("%.2f" % distance)
        if tool.checkFileEmpty(file_data) == 1 : 
            continue
        cmd = pathManage.scriptR() + "barplotTypeNumberOfneighbors.R " + file_data + " " + type + " " + str("%.2f" % distance)
        
        print cmd
        logFile.write(cmd + "\n")
        system(cmd)
    
    
def histProportion (distance, dir_in, logFile):
    """Draw proportion type neighbor
    in: Distance Max study, log path_file
    out: excecute CMD -> draw plot"""
    
    repScript = pathManage.scriptR()
    listType = structure.ListSub()
    listType.append("Global")
    listType.append("GlobalAmine")


    for type in listType:
        path_file = dir_in + "globalProportionAtom/"+ "proportionAtom" + type + str("%.2f" % distance)
        if tool.checkFileEmpty(path_file) == 1 : 
            continue
        cmd = repScript + "barplotQuantityGlobalAnalysis.R " + type + " " + path_file + " " + str(distance)
        print cmd
        system(cmd)

        logFile.write(cmd + "\n")


def histGlobalProportion(dir_out, logFile):
    """Proportion for each atom
    in: log file
    out: execute -> CMD
    """
    
    repScript = pathManage.scriptR()
    file = dir_out + "GlobalproportionCounterIonGlobal"
    if tool.checkFileEmpty(file) == 1 : 
        return
    cmd = repScript + "barplotProportionGlobalRef.R " + file
    logFile.write(cmd + "\n")
    print cmd
    system (cmd)


def barplotResDist (p_filin, logFile):
    """Draw barplot proportion residue
    in: log file_side
    out: execute CMD"""
    
    cmd = pathManage.scriptR() + "barplotResidueDistance.R " + p_filin + " all"
    logFile.write (cmd + "\n")
    system(cmd)
    print cmd
    

def histDistance(p_filin, type_distance):
    """Draw CN length histogram and coplar histogram
    in: nameFile data, option for execute R script
    out: execute R script -> draw histogram distance"""

    if tool.checkFileEmpty(p_filin) == 1 : 
        return
    cmd = pathManage.scriptR() + "distance.R " + str(p_filin) + " " + str(type_distance) 
    print cmd
    system(cmd)


def plotAngle(l_p_filin, debug = 1):
    """Excecute commande for draw angle plot
    in: distance MAX and log file
    out: Execute CMD -> draw plot
    """
    for p_filin in l_p_filin : 
# # # #         substruct = p_filin.split ("_")[-1]
        
        # do not need with the new criterion
# #         if substruct == "Imidazole" or substruct == "Pyridine" or substruct == "Secondary" : 
# #             cmd_3d = pathManage.scriptR() + "angle3D_Secondary.R " + p_filin
# #         elif substruct == "Tertiary" : 
# #             cmd_3d = pathManage.scriptR() + "angle3D_Tertiary.R " + p_filin    

        # every structure  
        cmd_density = pathManage.scriptR() + "angle_density.R " +  p_filin + "&" 
        cmd_distribution = pathManage.scriptR() + "angle_distribution.R " + p_filin + "&"
        cmdBarplot = pathManage.scriptR() + "angle_barplot.R " + p_filin + " 5"
        cmd_densityandhist = pathManage.scriptR() + "angleVSDistance.R " + p_filin
        
        if debug == 1: 
            print cmd_density
            print cmdBarplot
            print cmd_distribution
            print cmd_densityandhist
        
        system(cmdBarplot)
        system(cmd_density)
        system(cmd_distribution)
        system (cmd_densityandhist)
        
        

        
# # # # # # # # # # # #         if "cmd_3d" in locals() : 
# # # # # # # # # # # #             print cmd_3d
# # # # # # # # # # # #             system(cmd_3d)
# # # # # # # # # # # #             if debug ==1 : print cmd_3d

        
def waterPlotResolution (path_filin, verbose = 1) : 
    
    cmd = pathManage.scriptR() + "plotWater.R " + path_filin
    if verbose == 1 : 
        print cmd
    system (cmd)
    
        
        
def waterType (path_file, verbose = 1) : 
    
    cmd = pathManage.scriptR() + "plotWaterQuantity.R " + path_file
    if verbose == 1 : 
        print cmd
    system (cmd)    
    

def proportionAtomClassNeighbor (p_filin, logFile, verbose = 1):
    
    cmd = pathManage.scriptR() + "AnalysisNeighbor.R " + p_filin + " " + p_filin.split ("/")[-1].split ("_")[-1]
    system (cmd)
    if verbose == 1 : print cmd
    logFile.write (cmd + "\n")
    




def AFCBarplot (p_filin, logFile, verbose = 1):
    
    cmd = pathManage.scriptR() + "AFCNeighbor.R " + p_filin
    system (cmd)
    if verbose == 1 : print cmd
    logFile.write (cmd)

    

    
def histNeigbor (dir_in, logFile) : 
    
    l_study = structure.ListSub()
    l_study.append("global")
    
    for substruct in l_study : 
        cmd = pathManage.scriptR() + "AnalysisNeighbor.R " + dir_in + "neigbhor/" + "neighbor_" + substruct + " " + dir_in + "neigbhor/" + "distance_" + substruct + " " + substruct
        cmd_hist =  pathManage.scriptR() + "histAngle.R " + dir_in + "neigbhor/" + "angle_neighbor_" + substruct
        
        
        if substruct == "Primary" : 
            nb = 5
        if substruct == "Secondary" or substruct == "Imidazole" : 
            nb = 4
        if substruct == "Tertiary" : 
            nb = 3
        else : 
            nb = 7
        
        cmd_barplot = pathManage.scriptR() + "barplotNeighbor.R " + dir_in + "neigbhor/" + "barplot_" + substruct + " " + str (nb)
        
        print cmd 
        print cmd_hist
        print cmd_barplot
        
        system (cmd_barplot)
        system (cmd)
        system (cmd_hist)
        
    for i in range (1,8) : 
        AFC (dir_in + "neigbhor/", str (i))
    

def barplotLenBond (path_filin) : 
    
    cmd = pathManage.scriptR() + "boxplotBond.R " + path_filin
    print cmd 
    system (cmd)
    
    
def plot3D (p_file_coord, option = "local") : 
    
    cmd = pathManage.scriptR() + "scatter3D.R " + p_file_coord + " " + option
    print cmd 
    system (cmd)
    
    # remove png temp 
    system ("rm " + path.dirname(p_file_coord) + "/*.png")
    

def plotAngleVs (path_filin):    
    
    cmd = pathManage.scriptR() + "plotAngleVS.R " + path_filin
    print cmd
    system (cmd)
    
    
    
def AFC (pr_neighbors, number_neighbor):
    
    cmd = pathManage.scriptR() + "AFCNeighbor.R " + pr_neighbors + " " + str (number_neighbor)
    print cmd
    system (cmd)


def multiHist (pr_filin):
    
    cmd = pathManage.scriptR() + "multihist.R " + pr_filin 
    print cmd
    system (cmd)
    
    

def InteractionProportion (p_filin):
    
    
    cmd_by_sub =  pathManage.scriptR() + "proportionBySub.R " + p_filin
    print cmd_by_sub
    system (cmd_by_sub)
    
    cmd_all =  pathManage.scriptR() + "proportionAllSubs.R " + p_filin
    print cmd_all
    system (cmd_all)


def WaterMediated (p_filin):
    
    cmd =  pathManage.scriptR() + "barplotWaterMediated.R " + p_filin
    print cmd
    system (cmd)


def MergeProportionInteractAtLeasNotAtLeast (p_atleast, p_notatleast, pr_out):
    
    cmd_R = pathManage.scriptR() + "mergeBarplot.R " + str (p_atleast) + " " + str (p_notatleast) + " " + pr_out
    print cmd_R
    system (cmd_R)
    
    
def MergeProportionAndDataset (p_atleast1, p_notatleast1, p_atleast2, p_notatleast2, pr_out):
    
    
    cmd_R = pathManage.scriptR() + "mergeBarplotAndDataset.R " + str (p_atleast1) + " " + str (p_notatleast1) + " " + str (p_atleast2) + " " + str (p_notatleast2) + " " + str (pr_out)
    print cmd_R
    system (cmd_R)

def DistVSAngleNeighbor (p_filin):
    
    cmd_area = pathManage.scriptR() + "plotAngleVSDistCombi.R " + p_filin 
    print cmd_area
    
    system (cmd_area)


def HistAngle (p_filin):
    
    cmd_hist = pathManage.scriptR() + "plotAngleCriteria.R " + p_filin
    print cmd_hist
    
    system (cmd_hist)
    
     
def MeansNumberNeighbors (p_filout) : 
    
    cmd_run = pathManage.scriptR() + "barplotMeans.R " + str (p_filout)
    print cmd_run
    system(cmd_run)
    
    
        
def CorInteractionVSNbNeighbours (p_file_interact, p_file_nb_neighbor, pr_result) : 
    
    cmd_run = pathManage.scriptR() + "CorInteractVSneighbour.R " + p_file_interact + " " + p_file_nb_neighbor + " " + pr_result
    print cmd_run
    system (cmd_run)
    
def CorpKaVSCI (p_file_pka) : 
    
    cmd_run = pathManage.scriptR() + "corCIVSPKA.R " + p_file_pka
    
    print cmd_run
    system (cmd_run)
    
def CorpKaVSNb (p_file_pka):
    
    cmd_run = pathManage.scriptR() + "corNbVSPKA.R " + p_file_pka
    
    print cmd_run
    system (cmd_run)        



def ProtAnalysis (p_filin) : 
    
    p_filout = p_filin + "_sum"
    cmd_run = pathManage.scriptR() + "protSaltBridge.R " + p_filin + " > " + str (p_filout) 
    print cmd_run
    system (cmd_run)




