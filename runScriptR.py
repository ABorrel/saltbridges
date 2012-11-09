from os import system
from os import listdir
from re import search

import repertory
import structure
import log
from tool import checkFileEmpty as empty


def histStat(distance, type, file, typeStudyStruct, logFile):
    """Plot count statistic
    in: distance for legend plot, type histogram (coplar or length), file with data, type of study, logFile
    out: CMD in terminal -> plot """
    
    if empty(file) == 1 : 
        return
    rep = repertory.scriptR()
    cmd = rep + "barplotQuantity.R " + file + " " + str("%.2f" % distance) + " " + type + " " + typeStudyStruct
    logFile.write(cmd + "\n")
    system(cmd)


def histAA(distance, aminoAcid, file, logFile):
    """Plot file count amino acid proportion
    in: distance for legend plot, type histogram (coplar or length), file with data, type of study, logFile
    out: CMD in terminal -> plot """

    if empty(file) == 1 : 
        return
    rep = repertory.scriptR()
    repResult = repertory.resultAminoAcidGlobal()
    fileGlobal = repResult + "Global" + aminoAcid
    cmd = rep + "barplotQuantityAA.R " + file + " " + str("%.2f" % distance) + " " + aminoAcid + " " + str("%.2f" % distance)
    cmdGlobal = rep + "barplotQuantityAA.R " + fileGlobal + " " + str("%.2f" % distance) + " " + aminoAcid + " " + str("%.2f" % distance)
    logFile.write(cmd + "\n")
    logFile.write(cmdGlobal + "\n")
    system(cmd)
    system(cmdGlobal)


def plotDistanceOx(logFile):
    """Plot list distance
    in: distance for legend plot, type histogram (coplar or length), file with data, type of study, logFile
    out: CMD in terminal -> plot """

    repScript = repertory.scriptR()
    repFile = repertory.resultDistance()
    file = repFile + "resultDistanceOx"
    if empty(file) == 1 : 
        return
    cmd = repScript + "plotDistanceOx.R " + file
    logFile.write(cmd + "\n")
    system(cmd)


def globalStat(distanceAtoms, distanceResidues):
    """MAIN draw plot
    in: distance Atom study, distance Residues study"""

    timeStart, logFile = log.initAction("Run R Scripts")
    listStruct = structure.listStructure()
    listType = ["Ligands", "Residues", "Atoms"]
    listAminoAcid = ["ILE", "LEU", "LYS", "PHE", "TYR", "VAL", "SER", "MET", "ARG", "TRP", "PRO", "GLY", "GLU", "ASN", "HIS", "ALA", "ASP", "GLN", "THR", "CYS"]
    rep = repertory.result()


    for element in listStruct:
        repStruct = rep + element + "/"
        for type in listType:
            file = repStruct + "stat" + type + element
            if type == "atom":
                histStat(distanceAtoms, type, file, logFile)
            elif type == "Residues":
                distance = 2.00
                while distance <= distanceResidues:
                    fileTrace = file + str("%.2f" % distance)
                    histStat(distance, type, fileTrace, element, logFile)
                    histProportion(distance, logFile)
                    histProportionTypeNeighbors(distance, logFile)

                    distance = distance + 0.5

            else:
                histStat(distanceResidues, type, file, element, logFile)

        for aminoAcid in listAminoAcid:
            repAA = repertory.aminoAcid(element)
            file = repAA + element + aminoAcid
            histAA(distanceResidues, aminoAcid, file, logFile)

    plotDistanceOx(logFile)
    histGlobalProportion(logFile)
    histGlobalResidue(logFile)
    histProportionType(distanceResidues, logFile)
    histAtleastOne(distanceResidues, logFile)###################################3
    plotAngle(distanceResidues, logFile)
    plotAngle(distanceResidues, logFile)
    
    log.endAction("Run R Scripts", timeStart, logFile)


def histProportionType(distanceMax, logFile):
    """Draw proportion type graphe
    in: Distance Max study, log file
    out: excecute CMD -> draw plot"""
    
    
    repScript = repertory.scriptR()
    repFile = repertory.globalProportionType()
    listDistance = structure.listDistance(distanceMax)
    
    
    cmd = repScript + "barplotPercentClasse.R " + str(len(listDistance)) + " " + repFile
    logFile.write(cmd + "\n")
    system (cmd)

    
def histProportionTypeNeighbors(distance, logFile):
    """Draw proportion type neighbor
    in: Distance Max study, log file
    out: excecute CMD -> draw plot"""
    
    listStruct = structure.listStructure()
    listStruct.append("GlobalAmine")
    listStruct.append("Global")
    for type in listStruct : 
        file = repertory.globalProportionType() + "proportionType" + type + str("%.2f" % distance)
        if empty(file) == 1 : 
            continue
        cmd = repertory.scriptR() + "barplotTypeNumberOfneighbors.R " + repertory.globalProportionType() + "proportionType" + type + str("%.2f" % distance) + " " + type + " " + str("%.2f" % distance)
        
        logFile.write(cmd + "\n")
        system(cmd)
    
    
def histProportion (distance, logFile):
    """Draw proportion type neighbor
    in: Distance Max study, log file
    out: excecute CMD -> draw plot"""
    
    repScript = repertory.scriptR()
    listType = structure.listStructure()
    listType.append("Global")
    listType.append("GlobalAmine")


    for type in listType:
        repFile = repertory.globalProportionAtom()
        file = repFile + "proportionAtom" + type + str("%.2f" % distance)
        if empty(file) == 1 : 
            continue
        cmd = repScript + "barplotQuantityGlobalAnalysis.R " + type + " " + file + " " + str(distance)
        system(cmd)

        logFile.write(cmd + "\n")


def histGlobalProportion(logFile):
    """Proportion for each atom
    in: log file
    out: execute -> CMD
    """
    
    repScript = repertory.scriptR()
    file = repertory.result() + "GlobalproportionCounterIonGlobal"
    if empty(file) == 1 : 
        return
    cmd = repScript + "barplotProportionGlobal.R " + file
    logFile.write(cmd + "\n")
    system (cmd)


def histGlobalResidue(logFile):
    """Draw barplot proportion residue
    in: log file
    out: execute CMD"""
    
    listStructure = structure.listStructure()
    for element in listStructure:
        file = repertory.resultStruct(element) + "GlobalResidue" + element
        if empty (file) == 1 : 
            continue
        cmd = repertory.scriptR() + "barplotResidueDistance.R " + file + " " + element 
        logFile.write(cmd + "\n")
        system(cmd)
    fileGlobal = repertory.result() + "globalResidueAllAtoms"
    if empty(fileGlobal) == 1 : 
        return
    cmd2 = repertory.scriptR() + "barplotResidueDistance.R " + repertory.result() + "globalResidueAllAtoms" + " all"
    logFile.write(cmd2 + "\n")
    system(cmd2)


def histDistance(nameFile, type, base):
    """Draw CN length histogram and coplar histogram
    in: nameFile data, option for execute R script
    out: execute R script -> draw histogram distance"""


    file = str(repertory.resultDistance() + nameFile)
    if empty(file) == 1 : 
        return
    cmd = repertory.scriptR() + "distance.R " + str(repertory.resultDistance() + nameFile) + " " + str(type) + " " + str(base)
    system(cmd)


def histAtleastOne(distanceMax, logFile):
    """Draw at least one plot
    in: distance max, log file
    out: Excecute CMD -> draw all at least one plot
    """
    repResult = repertory.result()
    listFile = listdir(repResult)
    
    for file in listFile:
        
        if search("^atLeastOne_", file) :
            type = file.split("One_")[1]
            fileGobal = "atLeastOneGlobal" + file[10:]
            cmd = repertory.scriptR() + "barplotAtLeastOne.R " + repertory.result() + file + " " + repertory.result() + fileGobal + " " + str(distanceMax) + " " + type.upper()[0] + type[1:]
            logFile.write(cmd + "\n")
            system(cmd)


def plotAngle(distanceMax, logFile):
    """Excecute commande for draw angle plot
    in: distance MAX and log file
    out: Execute CMD -> draw plot
    """
    
    listStruct = structure.listStructure()
    for struct in listStruct : 
        if struct == "Imidazole" or struct == "Pyridine": 
            cmd = repertory.scriptR() + "angle_Secondary" + ".R " + repertory.resultAngle(struct) + "angle_" + str(struct)
            cmdBarplot = repertory.scriptR() + "angle_barplot.R " + repertory.resultAngle(struct) + "angle_" + str(struct) + " " + str(distanceMax)
        elif struct == "Diamine" : 
            cmd = repertory.scriptR() + "angle_Primary" + ".R " + repertory.resultAngle(struct) + "angle_" + str(struct)
            cmdBarplot = repertory.scriptR() + "angle_barplot.R " + repertory.resultAngle(struct) + "angle_" + str(struct) + " " + str(distanceMax)
        else : 
            cmd = repertory.scriptR() + "angle_" + str(struct) + ".R " + repertory.resultAngle(struct) + "angle_" + str(struct)
            cmdBarplot = repertory.scriptR() + "angle_barplot.R " + repertory.resultAngle(struct) + "angle_" + str(struct) + " " + str(distanceMax)
        
        logFile.write(cmd + "\n")
        logFile.write(cmdBarplot + "\n")
        system(cmd)
        system(cmdBarplot)
        
