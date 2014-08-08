from os import system, listdir, path
from re import search

import repertory
import structure
import log
import tool

def histStat(distance, type_studie, file, type_struct, logFile):
    """Plot count statistic
    in: distance for legend plot, type histogram (coplar or length), file with data, type of study, logFile
    out: CMD in terminal -> plot """
    
    if tool.checkFileEmpty(file) == 1 : 
        return
    rep = repertory.scriptR()
    cmd = rep + "barplotQuantity.R " + file + " " + str("%.2f" % distance) + " " + type_studie + " " + type_struct
    print cmd
    logFile.write(cmd + "\n")
    system(cmd)


def histAA(distance, aminoAcid, path_file, logFile):
    """Plot file count amino acid proportion
    in: distance for legend plot, type histogram (coplar or length), file with data, type of study, logFile
    out: CMD in terminal -> plot """

    if tool.checkFileEmpty(path_file) == 1 : 
        return
    rep = repertory.scriptR()
    cmd = rep + "barplotQuantityAA.R " + path_file + " " + str("%.2f" % distance) + " " + aminoAcid + " " + str("%.2f" % distance)
    print cmd
    logFile.write(cmd + "\n")
    system(cmd)


def plotDistanceOx(p_filin, logFile, debug = 1):
    """Plot list distance
    in: distance for legend plot, type histogram (coplar or length), file with data, type of study, logFile
    out: CMD in terminal -> plot """

    repScript = repertory.scriptR()
    if tool.checkFileEmpty(p_filin) == 1 : 
        return
    cmd = repScript + "plotDistanceOx.R " + p_filin
    if debug : print cmd
    logFile.write(cmd + "\n")
    system(cmd)


def globalStat(d_max, dir_in):
    """MAIN draw plot
    in: distance Atom study, distance Residues study"""

    timeStart, logFile = log.initAction("Run R Scripts")
    listStruct = structure.listStructure()
    listType = ["Ligands", "Residues", "Atoms"]
    listAminoAcid = ["ILE", "LEU", "LYS", "PHE", "TYR", "VAL", "SER", "MET", "ARG", "TRP", "PRO", "GLY", "GLU", "ASN", "HIS", "ALA", "ASP", "GLN", "THR", "CYS"]
 
 
    for type_struct in listStruct:
        repStruct = repertory.typeSubStructure(dir_in, type_struct)
        for type_studies in listType:
            path_file = repStruct + "stat" + type_studies + type_struct
            if type_studies == "Atoms":
                histStat(d_max, type_studies, path_file, type_struct, logFile)
            elif type_studies == "Residues":
                distance = 2.00
                while distance <= d_max:
                    fileTrace = path_file + str("%.2f" % distance)
                    histStat(distance, type_studies, fileTrace, type_struct, logFile)
                    histProportion(distance, dir_in, logFile)
                    histProportionTypeNeighbors(distance, dir_in, logFile)
  
                    distance = distance + 0.5
 
            else:
                histStat(d_max, type_studies, path_file, type_struct, logFile)
            for aminoAcid in listAminoAcid:
                path_file_aa = dir_in + type_struct + "/aminoAcid/" + type_struct + aminoAcid
                histAA(d_max, aminoAcid, path_file_aa, logFile)
#  
        for aminoAcid in listAminoAcid:
            path_file = dir_in + "AminoAcidGlobal/Global" + aminoAcid
            histAA(d_max, aminoAcid, path_file, logFile)
 
    plotDistanceOx(repertory.resultDistance(dir_in), logFile)
    ################################################## histGlobalProportion(dir_in, logFile)
    histGlobalResidue(dir_in, logFile)
    histProportionType(d_max, dir_in + "globalProportionType/", logFile)
    histAtleastOne(dir_in, logFile)  
    histNeigbor (dir_in, logFile)
    plotAngle(d_max, dir_in, logFile)
    
    log.endAction("Run R Scripts", timeStart, logFile)


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


def histGlobalResidue(dir_out, logFile):
    """Draw barplot proportion residue
    in: log file_side
    out: execute CMD"""
    
    listStructure = structure.listStructure()
    for element in listStructure:
        file_side = repertory.resultStruct(element, dir_out) + "GlobalResidueSide" + element
        file_global = repertory.resultStruct(element, dir_out) + "GlobalResidue" + element
        if tool.checkFileEmpty (file_side) == 1 and tool.checkFileEmpty (file_global): 
            continue
        cmd = repertory.scriptR() + "barplotResidueDistance.R " + file_side + " " + element 
        cmd_global = repertory.scriptR() + "barplotResidueDistance.R " + file_global + " " + element 
        logFile.write(cmd + "\n")
        logFile.write(cmd_global + "\n")
        system(cmd)
        system (cmd_global)
    file_all_atom_side = dir_out + "globalResidueAllAtomsSide"
    file_all_atom = dir_out + "globalResidueAllAtoms"
    if tool.checkFileEmpty(file_all_atom_side) == 1 or tool.checkFileEmpty(file_all_atom) == 1: 
        return
    cmd_all_side = repertory.scriptR() + "barplotResidueDistance.R " + dir_out + "globalResidueAllAtomsSide" + " all"
    cmd_all = repertory.scriptR() + "barplotResidueDistance.R " + dir_out + "globalResidueAllAtoms" + " all"
    logFile.write(cmd_all_side + "\n")
    logFile.write(cmd_all + "\n")
    system(cmd_all_side)
    system(cmd_all)
    print cmd
    print cmd_global
    print cmd_all_side
    print cmd_all


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


    
