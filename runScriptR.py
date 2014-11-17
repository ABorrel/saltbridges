from os import system, listdir, path
from re import search

import repertory
import structure
import log
import tool

def barplotQuantity(distance, type_study, p_filin, logFile, debug = 1):
    """Plot count statistic
    in: distance for legend plot, type histogram (coplar or length), p_filin with data, type of study, logFile
    out: CMD in terminal -> plot """
    
    if tool.checkFileEmpty(p_filin) == 1 : 
        return
    rep = repertory.scriptR()
    cmd = rep + "barplotQuantity.R " + p_filin + " " + str("%.2f" % distance) + " " + type_study
    if debug : print cmd
    logFile.write(cmd + "\n")
    system(cmd)


def barplotQuantityByAA(distance, aminoAcid, path_file, logFile):
    """to do"""

    if tool.checkFileEmpty(path_file) == 1 : 
        return
    rep = repertory.scriptR()
    cmd = rep + "barplotQuantityByAA.R " + path_file + " " + str("%.2f" % distance) + " " + aminoAcid + " " + str("%.2f" % distance)
    print cmd
    logFile.write(cmd + "\n")
    system(cmd)


def plotDistance(p_filin, logFile, debug = 1):
    """Plot list distance
    in: distance for legend plot, type histogram (coplar or length), file with data, type of study, logFile
    out: CMD in terminal -> plot """

    repScript = repertory.scriptR()
    if tool.checkFileEmpty(p_filin) == 1 : 
        return
    cmd = repScript + "plotDistance.R " + p_filin
    if debug : print cmd
    logFile.write(cmd + "\n")
    system(cmd)


def plotDistanceDensity(p_filin, logFile, debug = 1):
    """Plot density"""

    repScript = repertory.scriptR()
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

    repScript = repertory.scriptR()
    if tool.checkFileEmpty(p_filin) == 1 : 
        return
    cmd = repScript + "plotNumberNeighbor.R " + p_filin
    if debug : print cmd
    logFile.write(cmd + "\n")
    system(cmd)


def barplotCombination (p_filin, logFile, debug = 1):
    
    repScript = repertory.scriptR()
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
    
    
    repScript = repertory.scriptR()
    
    cmd = repScript + "barplotPercentClasse.R " + dir_out
    print cmd
    logFile.write(cmd + "\n")
    system (cmd)

    
def histProportionTypeNeighbors(distance, dir_out, logFile):
    """Draw proportion type neighbor
    in: Distance Max study, log file_data
    out: excecute CMD -> draw plot"""
    
    listStruct = structure.listStructure()
    listStruct.append("GlobalAmine")
    listStruct.append("Global")
    for type in listStruct : 
        file_data = dir_out + "globalProportionType/" + "proportionType" + type + str("%.2f" % distance)
        if tool.checkFileEmpty(file_data) == 1 : 
            continue
        cmd = repertory.scriptR() + "barplotTypeNumberOfneighbors.R " + file_data + " " + type + " " + str("%.2f" % distance)
        
        print cmd
        logFile.write(cmd + "\n")
        system(cmd)
    
    
def histProportion (distance, dir_in, logFile):
    """Draw proportion type neighbor
    in: Distance Max study, log path_file
    out: excecute CMD -> draw plot"""
    
    repScript = repertory.scriptR()
    listType = structure.listStructure()
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
    
    repScript = repertory.scriptR()
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
    
    cmd = repertory.scriptR() + "barplotResidueDistance.R " + p_filin + " all"
    logFile.write (cmd + "\n")
    system(cmd)
    print cmd
    

def histDistance(p_filin, type_distance):
    """Draw CN length histogram and coplar histogram
    in: nameFile data, option for execute R script
    out: execute R script -> draw histogram distance"""

    if tool.checkFileEmpty(p_filin) == 1 : 
        return
    cmd = repertory.scriptR() + "distance.R " + str(p_filin) + " " + str(type_distance) 
    print cmd
    system(cmd)


def histAtleastOne(dir_out, logFile):
    """Draw at least one plot
    in: distance max, log file
    out: Excecute CMD -> draw all at least one plot
    """
    repResult = dir_out + "AtLeastOne/"
    l_file = listdir(repResult)
    
    for name_file in l_file:
        if search(".dat", name_file) :
            type_atleastne = name_file.split(".")[0]
            cmd = repertory.scriptR() + "barplotAtLeastOne.R " + repResult + name_file + " " + type_atleastne
            logFile.write(cmd + "\n")
            print cmd
            system(cmd)


def plotAngle(l_p_filin, logFile, debug = 1):
    """Excecute commande for draw angle plot
    in: distance MAX and log file
    out: Execute CMD -> draw plot
    """
    for p_filin in l_p_filin : 
        substruct = p_filin.split ("_")[-1]
        
        if substruct == "Imidazole" or substruct == "Pyridine" or substruct == "Secondary" : 
            cmd_3d = repertory.scriptR() + "angle3D_Secondary.R " + p_filin
        elif substruct == "Tertiary" : 
            cmd_3d = repertory.scriptR() + "angle3D_Tertiary.R " + p_filin    

        # every structure  
        cmd_density = repertory.scriptR() + "angle_density.R " +  p_filin + "&" 
        cmd_distribution = repertory.scriptR() + "angle_distribution.R " + p_filin + "&"
        cmdBarplot = repertory.scriptR() + "angle_barplot.R " + p_filin + " 5"
        
        logFile.write(cmdBarplot + "\n")
        logFile.write(cmd_density + "\n")
        logFile.write(cmd_distribution + "\n")
        
        system(cmdBarplot)
        system(cmd_density)
        system(cmd_distribution)
        
        if debug == 1: 
            print cmd_density
            print cmdBarplot
            print cmd_distribution
        
        if "cmd_3d" in locals() : 
            logFile.write(cmd_3d + "\n")
            system(cmd_3d)
            if debug ==1 : print cmd_3d

        
def waterPlotResolution (path_filin, verbose = 1) : 
    
    cmd = repertory.scriptR() + "plotWater.R " + path_filin
    if verbose == 1 : 
        print cmd
    system (cmd)
    
        
        
def waterType (path_file, verbose = 1) : 
    
    cmd =    repertory.scriptR() + "plotWaterQuantity.R " + path_file
    if verbose == 1 : 
        print cmd
    system (cmd)    
    

def proportionAtomClassNeighbor (p_filin, logFile, verbose = 1):
    
    cmd = repertory.scriptR() + "AnalysisNeighbor.R " + p_filin + " " + p_filin.split ("/")[-1].split ("_")[-1]
    system (cmd)
    if verbose == 1 : print cmd
    logFile.write (cmd + "\n")
    

def AFCPieFirstNeighbor (p_filin, logFile, verbose = 1):
    
    cmd = repertory.scriptR() + "AFCNeighbor.R " + p_filin
    system (cmd)
    if verbose == 1 : print cmd
    logFile.write (cmd)

    

    
def histNeigbor (dir_in, logFile) : 
    
    l_study = structure.listStructure()
    l_study.append("global")
    
    for substruct in l_study : 
        cmd = repertory.scriptR() + "AnalysisNeighbor.R " + dir_in + "neigbhor/" + "neighbor_" + substruct + " " + dir_in + "neigbhor/" + "distance_" + substruct + " " + substruct
        cmd_hist =  repertory.scriptR() + "histAngle.R " + dir_in + "neigbhor/" + "angle_neighbor_" + substruct
        
        
        if substruct == "Primary" : 
            nb = 5
        if substruct == "Secondary" or substruct == "Imidazole" : 
            nb = 4
        if substruct == "Tertiary" : 
            nb = 3
        else : 
            nb = 7
        
        cmd_barplot = repertory.scriptR() + "barplotNeighbor.R " + dir_in + "neigbhor/" + "barplot_" + substruct + " " + str (nb)
        
        print cmd 
        print cmd_hist
        print cmd_barplot
        
        system (cmd_barplot)
        system (cmd)
        system (cmd_hist)
        
    for i in range (1,8) : 
        AFC (dir_in + "neigbhor/", str (i))
    

def barplotLenBond (path_filin) : 
    
    cmd = repertory.scriptR() + "boxplotBond.R " + path_filin
    print cmd 
    system (cmd)
    
    
def plot3D (p_file_coord) : 
    
    cmd = repertory.scriptR() + "scatter3D.R " + p_file_coord
    print cmd 
    system (cmd)
    
    # remove png temp 
    system ("rm " + path.dirname(p_file_coord) + "/*.png")
    

def plotAngleVs (path_filin):    
    
    cmd = repertory.scriptR() + "plotAngleVS.R " + path_filin
    print cmd
    system (cmd)
    
    
    
def AFC (pr_neighbors, number_neighbor):
    
    cmd = repertory.scriptR() + "AFCNeighbor.R " + pr_neighbors + " " + str (number_neighbor)
    print cmd
    system (cmd)


def multiHist (pr_filin):
    
    cmd = repertory.scriptR() + "multihist.R " + pr_filin 
    print cmd
    system (cmd)
    
    

def saltBridgesProportion (p_filin):
    
    cmd =  repertory.scriptR() + "pieProportion.R " + p_filin
    print cmd
    system (cmd)
    
    
      
    
