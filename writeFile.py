import structure



def resultFilterLigandPDB(struct, dir_out):

    list_path_file_dataset = []
    for resolutionKey in struct.keys():
        # open filout dataset
        fileWriteXRay = open(dir_out + "dataset_" + resolutionKey, "w")
        list_path_file_dataset.append (dir_out + "dataset_" + resolutionKey)
        if resolutionKey != "NMR" and resolutionKey != "OUT":
            fileWriteXRayRMN = open(dir_out + "dataset_" + resolutionKey + "_RMN", "w")
            for ligandNMR in struct["NMR"].keys():
                if not ligandNMR in struct[resolutionKey]:
                    fileWriteXRayRMN.write(str(ligandNMR) + "\t")
                    for PDBNMR in struct["NMR"][ligandNMR]:
                        fileWriteXRayRMN.write(str(PDBNMR) + " ")
                    fileWriteXRayRMN.write("\n")
        for ligand in struct[resolutionKey].keys():
            fileWriteXRay.write(str(ligand) + "\t")
            if resolutionKey != "NMR" and resolutionKey != "OUT":
                fileWriteXRayRMN.write(str(ligand) + "\t")
            for PDB in struct[resolutionKey][ligand]:
                fileWriteXRay.write(str(PDB) + " ")
                if resolutionKey != "NMR" and resolutionKey != "OUT":
                    fileWriteXRayRMN.write(str(PDB) + " ")
            if resolutionKey != "NMR" and resolutionKey != "OUT":
                if ligand in struct["NMR"]:
                    for PDBMNR in struct["NMR"][ligand]:
                        fileWriteXRayRMN.write(str(PDBMNR) + " ")
            fileWriteXRay.write("\n")
            if resolutionKey != "NMR" and resolutionKey != "OUT":
                fileWriteXRayRMN.write("\n")
        fileWriteXRay.close()
        if resolutionKey != "NMR" and resolutionKey != "OUT":
            fileWriteXRayRMN.close()
        
    return list_path_file_dataset



def resultLigandInPDB(structResult, directory_out):

    filout = open(directory_out + "resultLigandInPDB", "w")

    for pdb in structResult:
        filout.write(pdb["name"] + "\t")
        if pdb["ligands"] != []:
            for ligand in pdb["ligands"]:
                filout.write(ligand + " ")
        filout.write("\n")
    filout.close()
    return directory_out + "resultLigandInPDB"



def resultCoplanar(struct, name, directory_out):
    
    filout = open(directory_out + name, "w")

    for element in struct:
        filout.write("%.2f\n" % element)
    filout.close()



def resultLengthCNBond(structResult, file, directory_out):

    filout = open(directory_out + file, "w")

    for element in structResult:
        filout.write("%.2f\n" % element)

    filout.close()



def parsingDataSet(listCount, countAmine, numberPDB, path_file_dataset):

    filout = open(path_file_dataset  + ".stat", "w")
    maxRepetition = 0
    SumPDB = 0.0
    
    for count in listCount:
        filout.write(str(count["name"]) + "\t")
        filout.write(str(count["Number PDB"]))
        filout.write("\n")
        SumPDB = SumPDB + count["Number PDB"]
        if count["Number PDB"] > maxRepetition : 
            maxRepetition = count["Number PDB"] 
    
    nbCount = len(listCount) 
    filout.write("----------------------------------------\n")
    filout.write("Number ligand: " + str(nbCount) + "\n")
    filout.write("Number different PDB files: " + str(numberPDB) + "\n")
    filout.write("Repetition maximun ligand: " + str(maxRepetition) + "\n")
    filout.write("Average PDB by ligand: ")
    try : average = SumPDB / nbCount
    except : average = 0
    filout.write("%.3f\n" % average)
    filout.write("----------------------------------------\n")
    filout.write("Number Primary Amine: " + str(countAmine["Primary"]) + "\n")
    filout.write("Number Secondary Amine: " + str(countAmine["Secondary"]) + "\n")
    filout.write("Number Tertiary Amine: " + str(countAmine["Tertiary"]) + "\n")
    filout.write("Number Imidazole: " + str(countAmine["Imidazole"]) + "\n")
    filout.write("Number Guanidium: " + str(countAmine["Guanidium"]) + "\n")
    filout.write("Number Diamine: " + str(countAmine["Diamine"]) + "\n")
    filout.write("Number Pyridine: " + str(countAmine["Pyridine"]) + "\n")
    
    
    filout.close()



def resultDistanceOx(count, directory_out):

    filout = open(directory_out + "resultDistanceOx", "w")

    for element in count.keys():
        for distance in count[element]:
            filout.write(str("%.2f\t" % distance) + element + "\n")

    filout.close()



def amine(listAmine, files):

    for typeStructure in listAmine.keys():
        for nitrogen in listAmine[typeStructure]:
            lineWrite = str(nitrogen["PDB"]) + "\t" + str(nitrogen["serial"]) + "-" + str(nitrogen["resName"]) + "\t"
            for neighbor in nitrogen["neighbors"]:
                lineWrite = lineWrite + str(neighbor["serial"]) + " " + str(neighbor["name"]) + " " + str(neighbor["resName"]) + " " + str("%.2f" % neighbor["distance"])
                for angle in neighbor["angle"]:
                    lineWrite = lineWrite + " " + str("%.2f" % angle)
                lineWrite = lineWrite + "//"
            lineWrite = lineWrite + "\n"
            files[typeStructure].write(lineWrite)


def openFileAmine(directory_out):

    listS = structure.listStructure()

    dictFile = {}

    for element in listS:
        filout = open(directory_out + "summary" + element, "w")
        dictFile[element] = filout

    return dictFile



def closeFileAmine(dictFile):

    for key in dictFile.keys():
        dictFile[key].close()


def countGlobalAmine(distanceGlobal, countGlobalAmine, dir_out):

    resultDistanceOx(countGlobalAmine[str(distanceGlobal)]["distanceOx"], dir_out)
    resultLigand(countGlobalAmine[str(distanceGlobal)]["ligand"], dir_out)
    resultAtom(countGlobalAmine[str(distanceGlobal)]["atom"], dir_out)
    resultByAA(countGlobalAmine[str(distanceGlobal)]["byAA"], dir_out)
    resultAngle(countGlobalAmine, distanceGlobal)

    for distance in countGlobalAmine.keys():
        resultResidue(float(distance), countGlobalAmine[str(distance)]["residue"], dir_out)
        resultProportion(float(distance), countGlobalAmine[str(distance)]["proportionAtom"], dir_out)
        resultProportionType(float(distance), countGlobalAmine[str(distance)]["proportionType"], dir_out)

        
    resultAtLeastOne(countGlobalAmine, distanceGlobal, dir_out)
    resultProportionGlobalType(countGlobalAmine, distanceGlobal, dir_out)
    resultGlobalResidue(countGlobalAmine, distanceGlobal, dir_out)
    resultResidueDistance(countGlobalAmine, distanceGlobal, dir_out)
    

def resultLigand(count, directory_out):

    listStructure = structure.listStructure()

    for element in listStructure:
        filout = open(directory_out + "statLigands" + element, "w")
        for ligand in count[element].keys():
            filout.write(str(ligand) + "\t" + str(count[element][ligand]) + "\n")

        filout.close()


def resultAtom(count, directory_out):

    listStructure = structure.listStructure()

    for element in listStructure:
        filout = open(directory_out + "statAtoms" + element, "w")
        for atom in count[element].keys():
            filout.write(str(atom) + "\t" + str(count[element][atom]) + "\n")

        filout.close()


def resultResidue(distance, count, directory_out):

    listStructure = structure.listStructure()
    
    for element in listStructure:
        filout = open(directory_out + "statResidues" + element + str("%.2f" % distance), "w")
        for residue in count[element].keys():
            filout.write(str(residue) + "\t" + str(count[element][residue]["main"]) + "\t" + str(count[element][residue]["side"]) + "\n")

        filout.close()


def resultByAA(count, directory_out):

    listStructure = structure.listStructure()

    for element in listStructure:
        for aminoAcid in count[element].keys():
            filout = open(directory_out + element + aminoAcid, "w")
            for atom in count[element][aminoAcid].keys():
                lineWrite = str(atom) + "\t" + str(count[element][aminoAcid][atom]["3.5"]) + "\t" + str(count[element][aminoAcid][atom]["4.5"]) + "\n"
                filout.write(lineWrite)
            filout.close()

    for aminoAcid in count["global"].keys():
        filout = open(directory_out + "Global" + aminoAcid, "w")
        for atom in count["global"][aminoAcid].keys():
            lineWrite = str(atom) + "\t" + str(count["global"][aminoAcid][atom]["3.5"]) + "\t" + str(count["global"][aminoAcid][atom]["4.5"]) + "\n"
            filout.write(lineWrite)
        filout.close()


def resultProportion (distance, count, directory_out):
    """Write file proportion number of neighbors"""


    for type in count.keys():
        filout = open (directory_out + "proportionAtom" + type + str("%.2f" % distance), "w")
        for nbNeighbor in count[type].keys():
            filout.write(str(nbNeighbor) + "\t" + str(count[type][nbNeighbor]["C"]) + "\t" + str(count[type][nbNeighbor]["O"]) + "\t" + str(count[type][nbNeighbor]["N"]) + "\t" + str(count[type][nbNeighbor]["S"]) + "\t" + str(count[type][nbNeighbor]["others"]) + "\n")


        filout.close()
        
        
def resultProportionType (distance, count, directory_out):

    listClasse = ["OxAcid", "amphiprotic", "OxAccept", "Carom", "H2O", "Ndonnor", "Nbasic", "others"]
    for type in count.keys():
        filout = open (directory_out + "proportionType" + type + str("%.2f" % distance), "w")

        for nbNeighbor in count[type].keys():
            if not nbNeighbor == "allNumberNeighbors" : 
                filout.write(str(nbNeighbor))
                for classe in listClasse : 
                    filout.write("\t" + str(count[type][nbNeighbor][classe]))
                filout.write("\n")
                
        if count[type].keys() == [] : 
            filout.write(str(0) + "\t" + str(0) + "\t" + str(0) + "\t" + str(0) + "\t" + str(0) + "\t" + str(0) + "\n")

        filout.close()


def resultProportionGlobalType(count, distanceGlobal,directory_out) : 
    
    listDistance = structure.listDistance(distanceGlobal)
    listClasse = ["OxAcid", "amphiprotic", "H2O", "OxAccept", "Nbasic", "Ndonnor", "Carom", "others"]
    listStruct = structure.listStructure()
    listStruct.append("Global") 
    
    for type in listStruct : 
        filout = open (directory_out + "proportionType" + type , "w")
        
        for distance in listDistance :
            for classe in listClasse : 
                filout.write(str(count[distance]["proportionType"][type]["allNumberNeighbors"][classe]) + "\t")
            filout.write("\n")
            
        filout.close()
    
    
    

def resultResidueDistance(countGlobalAmine, distanceMax, directory_out):

    listDistance = structure.listDistance(distanceMax)
    listAminoAcid = ["HOH", "ASP", "GLU", "THR", "SER", "ASN", "GLN", "TYR", "HIS", "LYS", "ARG", "PHE", "TRP", "ALA", "ILE", "LEU", "MET", "VAL", "CYS", "GLY", "PRO"]

    for type in  countGlobalAmine["2.5"]["residue"].keys():
        filout = open(directory_out + str("Global") + "Residue" + str(type), "w")
        for aminoAcid in listAminoAcid:
            line = str(aminoAcid)
            for distance in listDistance:
                line = line + "\t"
                if aminoAcid == "HOH":
                    line = line + str(countGlobalAmine[distance]["residue"][type][aminoAcid]["main"])
                else:
                    line = line + str(countGlobalAmine[distance]["residue"][type][aminoAcid]["side"])

            line = line + "\n"
            filout.write(line)
        filout.close()



def resultAtLeastOneGlobal(countGlobal, distanceMax, directory_out):

    listDistance = structure.listDistance(distanceMax)
    listStudyAtLeastOne = countGlobal[str(distanceMax)]["atLeastOne"].keys()
    
    for StudyAtleastOne in listStudyAtLeastOne : 
        filout = open(directory_out + "atLeastOneGlobal_" + StudyAtleastOne, "w")
        line = ""
        for distance in listDistance :  
            line = line + str(countGlobal[distance]["atLeastOne"][StudyAtleastOne][StudyAtleastOne]) + "\t" + str(countGlobal[distance]["atLeastOne"][StudyAtleastOne]["others"]) + "\t"
        line = line + "\n" 
        
        filout.write(line)
        filout.close()


def resultGlobalResidue(count, distanceMax, directory_out):

    listDistance = structure.listDistance(distanceMax)
    listAminoAcid = ["HOH", "ASP", "GLU", "THR", "SER", "ASN", "GLN", "TYR", "HIS", "LYS", "ARG", "PHE", "TRP", "ALA", "ILE", "LEU", "MET", "VAL", "CYS", "GLY", "PRO"]
    filout = open(directory_out + str("global") + "ResidueAllAtoms", "w")
    for aminoAcid in listAminoAcid:
        line = str(aminoAcid)
        for distance in listDistance:
            if aminoAcid == "HOH":
                line = line + "\t" + str(count[distance]["ResidueAllAtom"][aminoAcid]["main"])
            else:
                line = line + "\t" + str(count[distance]["ResidueAllAtom"][aminoAcid]["side"])

        line = line + "\n"
        filout.write(line)


def resultAtLeastOne(count, distanceMax, directory_out):

    listDist = structure.listDistance(distanceMax)
    listStructureStudy = structure.listStructure()
    listStudyAtLeastOne = count[str(distanceMax)]["atLeastOne"].keys()
    
    for StudyAtleastOne in listStudyAtLeastOne : 
        filout = open(directory_out + "atLeastOne_" + StudyAtleastOne, "w")
        line = ""
        for structureStudy in listStructureStudy :
            for distance in listDist :  
                line = line + str(count[distance]["atLeastOne"][StudyAtleastOne][structureStudy][StudyAtleastOne]) + "\t" + str(count[distance]["atLeastOne"][StudyAtleastOne][structureStudy]["others"]) + "\t"
            line = line + "\n" 
        
        filout.write(line)
        filout.close()



def resultAngle(count, distanceMax, directory_out):
    
    for type in count[str(distanceMax)]["angle"].keys():
        filoutGlobal = open(directory_out + "angle_" + str(type), "w")
        
        for classe in  count[str(distanceMax)]["angle"][type].keys():
            nbDistance = len(count[str(distanceMax)]["angle"][type][classe]["distance"])
            for i in range(0, nbDistance) : 
                distanceAt = count[str(distanceMax)]["angle"][type][classe]["distance"][i]
                filoutGlobal.write("%.2f" % distanceAt)
                
                for angle in count[str(distanceMax)]["angle"][type][classe]["angles"][i] : 
                    filoutGlobal.write("\t%.2f" % angle)
                    
                filoutGlobal.write("\t" + classe + "\n")
        filoutGlobal.close()
    countType = structure.countAngleType(count)  
    countAngle(countType) 
    

def countAngle (count, directory_out):
    
    listClasse = ["OxAcid", "amphiprotic", "H2O", "OxAccept", "Nbasic", "Ndonnor", "Carom", "others"]
    for type in count.keys() : 
        for distance in count[type].keys() : 
            filout = open(directory_out + "angle_" + type + "_" + distance, "w")
            for angle in count[type][distance].keys() : 
                filout.write(str(angle))
                for classe in listClasse : 
                    filout.write("\t" + str(count[type][distance][angle][classe]))
                filout.write("\n")
                
        filout.close()    
           


def openFilesWithoutSummary(distanceMax, directory_out):
    
    fileClass = {}
    
    listDistance = structure.listDistance(distanceMax)
    listStructureStudy = structure.listStructure()
    
    for distance in listDistance : 
        fileClass[distance] = {}
        for struct in listStructureStudy :
            fileClass[distance][struct] = open(directory_out + struct + "_<" + distance, "w")
    
    
    return fileClass
    
    
    
def withoutAtLeastOneSummary(amine, filesWithout, distance):
    
    listStudy = amine.keys()
    
    for study in listStudy : 
        for nitrogen in amine[study] : 
            flag = 0 
            for neighbors in nitrogen["neighbors"] : 
                if neighbors["classificationAtLeastOne"] != "others" : 
                    flag = 1
            if flag == 0 : 
                filesWithout[str(distance)][study].write(str(nitrogen["PDB"]) + "\t" + str(nitrogen["serial"]) + "\n")


def closeFilesWithoutSummary (filesWithoutSummary):
    
    for distance in filesWithoutSummary.keys() : 
        for study in filesWithoutSummary[distance].keys() : 
            filesWithoutSummary[distance][study].close()

           
