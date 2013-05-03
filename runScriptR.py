from os import system, listdir
from re import search

import repertory
import structure
import log
from tool import checkFileEmpty as empty


def histStat(distance, type_studie, file, type_struct, logFile):
    """Plot count statistic
    in: distance for legend plot, type histogram (coplar or length), file with data, type of study, logFile
    out: CMD in terminal -> plot """
    
    if empty(file) == 1 : 
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

    if empty(path_file) == 1 : 
        return
    rep = repertory.scriptR()
    cmd = rep + "barplotQuantityAA.R " + path_file + " " + str("%.2f" % distance) + " " + aminoAcid + " " + str("%.2f" % distance)
    #print cmd
    logFile.write(cmd + "\n")
    system(cmd)


def plotDistanceOx(rep_out, logFile):
    """Plot list distance
    in: distance for legend plot, type histogram (coplar or length), file with data, type of study, logFile
    out: CMD in terminal -> plot """

    repScript = repertory.scriptR()
    file = rep_out + "resultDistanceOx"
    if empty(file) == 1 : 
        return
    cmd = repScript + "plotDistanceOx.R " + file
    logFile.write(cmd + "\n")
    system(cmd)


def globalStat(distanceAtoms, distanceResidues, dir_in):
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
                histStat(distanceAtoms, type_studies, path_file, type_struct, logFile)
            elif type_studies == "Residues":
                distance = 2.00
                while distance <= distanceResidues:
                    fileTrace = path_file + str("%.2f" % distance)
                    histStat(distance, type_studies, fileTrace, type_struct, logFile)
                    histProportion(distance, dir_in, logFile)
                    histProportionTypeNeighbors(distance, dir_in, logFile)

                    distance = distance + 0.5

            else:
                histStat(distanceResidues, type_studies, path_file, type_struct, logFile)
            for aminoAcid in listAminoAcid:
                path_file_aa = dir_in + type_struct + "/aminoAcid/" + type_struct + aminoAcid
                histAA(distanceResidues, aminoAcid, path_file_aa, logFile)

        for aminoAcid in listAminoAcid:
            path_file = dir_in + "AminoAcidGlobal/Global" + aminoAcid
            histAA(distanceResidues, aminoAcid, path_file, logFile)

    plotDistanceOx(repertory.resultDistance(dir_in), logFile)
    ################################################## histGlobalProportion(dir_in, logFile)
    histGlobalResidue(dir_in, logFile)
    histProportionType(distanceResidues, dir_in + "globalProportionType/", logFile)
    histAtleastOne(distanceResidues, dir_in, logFile)  ###################################3
    plotAngle(distanceResidues, dir_in, logFile)
    
    log.endAction("Run R Scripts", timeStart, logFile)


def histProportionType(distanceMax, dir_out, logFile):
    """Draw proportion type graphe
    in: Distance Max study, log file
    out: excecute CMD -> draw plot"""
    
    
    repScript = repertory.scriptR()
    listDistance = structure.listDistance(distanceMax)
    
    
    cmd = repScript + "barplotPercentClasse.R " + str(len(listDistance)) + " " + dir_out
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
        file_data = dir_out + "/globalProportionType/" + "proportionType" + type + str("%.2f" % distance)
        if empty(file_data) == 1 : 
            continue
        cmd = repertory.scriptR() + "barplotTypeNumberOfneighbors.R " + file_data + " " + type + " " + str("%.2f" % distance)
        
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
        path_file = dir_in + "/globalProportionAtom/"+ "proportionAtom" + type + str("%.2f" % distance)
        if empty(path_file) == 1 : 
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
    if empty(file) == 1 : 
        return
    cmd = repScript + "barplotProportionGlobalRef.R " + file
    logFile.write(cmd + "\n")
    system (cmd)


def histGlobalResidue(dir_out, logFile):
    """Draw barplot proportion residue
    in: log file
    out: execute CMD"""
    
    listStructure = structure.listStructure()
    for element in listStructure:
        file = repertory.resultStruct(element, dir_out) + "GlobalResidue" + element
        if empty (file) == 1 : 
            continue
        cmd = repertory.scriptR() + "barplotResidueDistance.R " + file + " " + element 
        logFile.write(cmd + "\n")
        system(cmd)
    fileGlobal = dir_out + "globalResidueAllAtoms"
    if empty(fileGlobal) == 1 : 
        return
    cmd2 = repertory.scriptR() + "barplotResidueDistance.R " + dir_out + "globalResidueAllAtoms" + " all"
    logFile.write(cmd2 + "\n")
    system(cmd2)
    print cmd
    print cmd2


def histDistance(nameFile, type_distance, base, dir_out):
    """Draw CN length histogram and coplar histogram
    in: nameFile data, option for execute R script
    out: execute R script -> draw histogram distance"""


    file_distance = str(dir_out + nameFile)
    if empty(file_distance) == 1 : 
        return
    cmd = repertory.scriptR() + "distance.R " + str(dir_out + nameFile) + " " + str(type_distance) + " " + str(base)
    print cmd
    system(cmd)


def histAtleastOne(distanceMax, dir_out, logFile):
    """Draw at least one plot
    in: distance max, log file
    out: Excecute CMD -> draw all at least one plot
    """
    repResult = dir_out
    listFile = listdir(repResult)
    
    for file in listFile:
        if search("^atLeastOne_", file) :
            type = file.split("One_")[1]
            fileGobal = "atLeastOneGlobal" + file[10:]
            cmd = repertory.scriptR() + "barplotAtLeastOne.R " + dir_out + file + " " + dir_out + fileGobal + " " + str(distanceMax) + " " + type.upper()[0] + type[1:]
            logFile.write(cmd + "\n")
            system(cmd)


def plotAngle(distanceMax, dir_out, logFile):
    """Excecute commande for draw angle plot
    in: distance MAX and log file
    out: Execute CMD -> draw plot
    """
    
    listStruct = structure.listStructure()
    for struct in listStruct : 
        if struct == "Imidazole" or struct == "Pyridine": 
            cmd = repertory.scriptR() + "angle_Secondary" + ".R " + repertory.resultAngle(struct, dir_out) + "angle_" + str(struct)
            cmdBarplot = repertory.scriptR() + "angle_barplot.R " + repertory.resultAngle(struct, dir_out) + "angle_" + str(struct) + " " + str(distanceMax)
        elif struct == "Diamine" or struct == "AcidCarboxylic" or struct ==  "Guanidium": 
            cmd = repertory.scriptR() + "angle_Primary.R " + repertory.resultAngle(struct, dir_out) + "angle_" + str(struct)
            cmdBarplot = repertory.scriptR() + "angle_barplot.R " + repertory.resultAngle(struct, dir_out) + "angle_" + str(struct) + " " + str(distanceMax)
        else : 
            cmd = repertory.scriptR() + "angle_" + str(struct) + ".R " + repertory.resultAngle(struct, dir_out) + "angle_" + str(struct)
            cmdBarplot = repertory.scriptR() + "angle_barplot.R " + repertory.resultAngle(struct, dir_out) + "angle_" + str(struct) + " " + str(distanceMax)
            
        cmd_density = repertory.scriptR() + "angle_density.R " + repertory.resultAngle(struct, dir_out) + "angle_" + str(struct)
        
#         print cmd
#         print cmdBarplot
        logFile.write(cmd + "\n")
        logFile.write(cmdBarplot + "\n")
        system(cmd)
        system(cmdBarplot)
        system(cmd_density)
        

def waterPlotResolution (path_filin, verbose = 1) : 
    
    cmd = repertory.scriptR() + "plotWater.R " + path_filin
    
    if verbose == 1 : 
        print cmd
        
    system (cmd)
    
        
        
    
    

