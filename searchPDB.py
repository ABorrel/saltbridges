from os import listdir
from re import search
from copy import deepcopy
import calcul
import loadFile
import log
import parsing
import repertory
import retrieveAtom
import structure
import tool
import toolSubstructure
import writeFile



#############global search######################
def interestStructure (listAtomLigand):
    """For one nitrogen atom search substructure
    in : list atom ligand
    out : list of substructures found (list string)"""

    listSerialNitroge = listNitrogen (listAtomLigand)
    structureFound = []

    for nitrogen in listSerialNitroge:
        listAtomConnectNitrogen, connect = retrieveAtom.atomConnect(listAtomLigand, nitrogen)
        if imidazole(listAtomConnectNitrogen, listAtomLigand) == 1:
            structureFound.append("Imidazole")
        elif cn (listAtomConnectNitrogen, listAtomLigand) == 1:
            structureFound.append("Primary")
        elif cnc(listAtomConnectNitrogen, listAtomLigand) == 1:
            structureFound.append("Secondary")
        elif cncc(listAtomConnectNitrogen, listAtomLigand) == 1:
            structureFound.append("Tertiary")
        elif guanidium(listAtomConnectNitrogen, listAtomLigand) == 1 : 
            structureFound.append("Guanidium")
        elif pyridine(listAtomConnectNitrogen, listAtomLigand) == 1 : 
            structureFound.append("Pyridine")
        elif diAmine(listAtomConnectNitrogen, listAtomLigand) == 1 : 
            structureFound.append("Diamine") 


    return structureFound


def ligands():
    '''search ligands in PDB database
    out : list of ligands with PDB files associated'''

    start, fileLog = log.initAction("Search ligands in PDB")

    repIn = repertory.openPdbFile()
    listPDB = listdir(repIn)
    listPDBLigand = []

    for PDBFile in listPDB:
        namePDB = PDBFile.split(".")[0]
        fileLog.write(namePDB + "\n")
        ligandInPDB = structure.ligandPDB()
        linesPDB = loadFile.openPdbFile(namePDB)
        ligandInPDB["name"] = namePDB

        for linePDB in linesPDB:
            if(search ("^HETATM", linePDB)):
                atom = parsing.lineCoords(linePDB)
                if not atom["resName"] in ligandInPDB["ligands"]:
                    ligandInPDB["ligands"].append(atom["resName"])

        listPDBLigand.append(ligandInPDB)

    writeFile.resultLigandInPDB(listPDBLigand)
    log.endAction("Search ligands in PDB", start, fileLog)




################substructure search -> amine###############

def cn(listAtomConnectNitrogen, listAtomLigan):
    """search primary amine 
    in: list atom connected of nitrogen, list atom ligand
    out: 1 or 0"""

    amine = toolSubstructure.matrixElement(listAtomConnectNitrogen)
    if amine == ["N", "C"] or amine == ["N", "C"]:
        if toolSubstructure.checkSingleBond(listAtomConnectNitrogen[0], listAtomConnectNitrogen[1]) == 1:
            if toolSubstructure.checkConectOnlyC(listAtomConnectNitrogen[1], listAtomLigan) == 1 : 
                return 1
    return 0


def cnc(listAtomConnectNitrogen, listAtomLigand):
    """search secondary amine 
    in: list atom connected of nitrogen, list atom ligand
    out: 1 or 0"""

    amine = toolSubstructure.matrixElement(listAtomConnectNitrogen)

    if amine == ["N", "C", "C"]:
        if toolSubstructure.checkConectOnlyC(listAtomConnectNitrogen[1], listAtomLigand) == 1 and toolSubstructure.checkConectOnlyC(listAtomConnectNitrogen[2], listAtomLigand) == 1:
            if toolSubstructure.checkSingleBond(listAtomConnectNitrogen[0], listAtomConnectNitrogen[1]) == 1 and toolSubstructure.checkSingleBond(listAtomConnectNitrogen[0], listAtomConnectNitrogen[2]) == 1:
                return 1
    return 0


def cncc(listAtomConnectNitrogen, listAtomligand):
    """search tertiary amine 
    in: list atom connected of nitrogen, list atom ligand
    out: 1 or 0"""

    amine = toolSubstructure.matrixElement(listAtomConnectNitrogen)

    if amine == ["N", "C", "C", "C"]:
        if toolSubstructure.checkConectOnlyC(listAtomConnectNitrogen[1], listAtomligand) == 1 and toolSubstructure.checkConectOnlyC(listAtomConnectNitrogen[2], listAtomligand) == 1 and toolSubstructure.checkConectOnlyC(listAtomConnectNitrogen[3], listAtomligand) == 1:
            if toolSubstructure.checkCoplanar(listAtomConnectNitrogen[0], listAtomligand) == 1:
                if toolSubstructure.checkSingleBond(listAtomConnectNitrogen[0], listAtomConnectNitrogen[1]) == 1 and toolSubstructure.checkSingleBond(listAtomConnectNitrogen[0], listAtomConnectNitrogen[2]) == 1 and toolSubstructure.checkSingleBond(listAtomConnectNitrogen[0], listAtomConnectNitrogen[3]) == 1:
                    return 1
    return 0


def guanidium(listConnectNitrogen, listAtomLigand):
    """search guanidium 
    in: list atom connected of nitrogen, list atom ligand
    out: 1 or 0"""
       
    amine = toolSubstructure.matrixElement(listConnectNitrogen)

    if amine == ["N", "C"] :
        atomC = listConnectNitrogen[1]
        findNH = 0
        atomC, connectC = retrieveAtom.atomConnect(listAtomLigand , int (atomC["serial"]))
        if connectC == ["C", "N", "N", "N"] :
            
            groupAtomN1, conect_N1 = retrieveAtom.atomConnect(listAtomLigand , int (atomC[0]["connect"][1]))
            groupAtomN2, conect_N2 = retrieveAtom.atomConnect(listAtomLigand , int (atomC[0]["connect"][2]))
            groupAtomN3, conect_N3 = retrieveAtom.atomConnect(listAtomLigand , int (atomC[0]["connect"][3]))
            list_conect = [conect_N1, conect_N2, conect_N3]
            
            for conect in list_conect : 
                print conect
                if conect == ["N", "C"] : 
                    findNH = findNH + 1
                else :     
                    for elementConect in conect[1:] : 
                        if elementConect != "C"  : 
                            return 0
            if findNH >= 2  : 
                return 1
    return 0


def guanidiumNitrogenCentral(listConnectNitrogen, listAtomLigand):
    """Search guanidium 
    in: list Atom connected of nitrogen, list global atoms
    out: boolean"""
    
    findNH = 0
    for atomConnectNitrogen in listConnectNitrogen : 
        if atomConnectNitrogen ["element"] == "C" : 
            listAtomC, connectElementC = retrieveAtom.atomConnect(listAtomLigand, atomConnectNitrogen["serial"])
            if connectElementC == ["C", "N", "N", "N"] : 
                for atomConnectC in listAtomC[1:] : 
                    nitrogenConnect, connectNitrogen = retrieveAtom.atomConnect(listAtomLigand, atomConnectC["serial"])
                    if connectNitrogen == ["N", "C"] : 
                        findNH = findNH + 1
    
    if findNH == 2 : 
        return 1
    
    return 0
                
    

def pyridine(listAtomConnectNitrogen, listAtomLigand):
    """search pyridine 
    in: list atom connected of nitrogen, list atom ligand
    out: boolean"""
        
    amine = toolSubstructure.matrixElement(listAtomConnectNitrogen)
    
    if amine == ["N", "C", "C"] : 
        nitrogenInit = listAtomConnectNitrogen[0]
        if  cycleOnlyTestCarbon(nitrogenInit["serial"], nitrogenInit["serial"], nitrogenInit["serial"], listAtomLigand, 6, 0) == 1 : 
            return 1
        
         
def diAmine (listAtomConnectNitrogen, listAtomLigand):
    """search diamine 
    in: list atom connected of nitrogen, list atom ligand
    out: boolean"""
        
    amine = toolSubstructure.matrixElement(listAtomConnectNitrogen)
    if amine == ["N", "C"] :
        carbon1, connectC1 = retrieveAtom.atomConnect(listAtomLigand, listAtomConnectNitrogen[0]["connect"][1]) 
        if connectC1 == ["C", "C", "N", "N"] or connectC1 == ["C", "N", "C", "N"] or connectC1 == ["C", "N", "N", "C"] : 
            for atomConectC1 in carbon1[0]["connect"][1:] : 
                element, connectElement = retrieveAtom.atomConnect(listAtomLigand, int(atomConectC1))
                if connectElement == ["C", "C", "C", "C"]:
                    nbDoubleBond = 0
                    atomC1, connectC1 = retrieveAtom.atomConnect(listAtomLigand, element[0]["connect"][1])
                    atomC2, connectC2 = retrieveAtom.atomConnect(listAtomLigand, element[0]["connect"][2])
                    atomC3, connectC3 = retrieveAtom.atomConnect(listAtomLigand, element[0]["connect"][3])
                    
                    if toolSubstructure.checkSingleBond(element[0] , atomC1[0]) == 1  : 
                        nbDoubleBond = nbDoubleBond + 1 
                    if toolSubstructure.checkSingleBond(element[0] , atomC2[0]) == 1  : 
                        nbDoubleBond = nbDoubleBond + 1
                    if toolSubstructure.checkSingleBond(element[0] , atomC3[0]) == 1  : 
                        nbDoubleBond = nbDoubleBond + 1
                        
                    if nbDoubleBond >= 1 : 
                        return 1
                    
    return 0
                    
######################################ribose #########################################
                    
def riboseFromADPRibose (listAtomLigand):
    """search ribose global 
    in: list atom connected of nitrogen, list atom ligand
    out: boolean"""
            
    for atom in listAtomLigand : 
        if atom["element"] == "O" :
            listAtomConnectO, connectO = retrieveAtom.atomConnect(listAtomLigand, atom["serial"])
            if connectO == ["O", "C", "C"] : 
                serialInit = atom["serial"]
                listAtomConnectC1, connectC1 = retrieveAtom.atomConnect(listAtomLigand, atom["connect"][1])
                listAtomConnectC2, connectC2 = retrieveAtom.atomConnect(listAtomLigand, atom["connect"][2])
                
                if "C" in connectC1[1:] : 
                    if riboseFromADPRiboseAtom2(listAtomConnectC1, listAtomLigand, serialInit) == 1 : 
                        return 1
                if "C" in connectC2[1:] : 
                    if riboseFromADPRiboseAtom2(listAtomConnectC2, listAtomLigand, serialInit) == 1 : 
                        return 1
    
    return 0
            
    
def riboseFromADPRiboseAtom2 (listAtomConnectC, listAtomLigand, serialInit) :
    """search ribose atom2 
    in: list atom connected of nitrogen, list atom ligand
    out: boolean"""
        
    for atom in listAtomConnectC :
        if atom["element"] == "C" :
            listAtomC3, connectC = retrieveAtom.atomConnect(listAtomLigand, atom["serial"])
        
            if connectC == ["C", "O", "C", "C"] or connectC == ["C", "C", "O", "C"] or connectC == ["C", "C", "C", "O"] : 
                serialNext = listAtomC3[0]["serial"]
                if riboseFromADPRiboseAtom3 (serialNext, listAtomLigand, serialInit, atom["serial"]) == 1 : 
                    return 1
    return 0
    

def riboseFromADPRiboseAtom3 (serialAtom3, listAtomLigand, serialInit, serialPrevious) : 
    """search ribose atom3 
    in: list atom connected of nitrogen, list atom ligand
    out: boolean"""
         
    listAtomConnect3, connect3 = retrieveAtom.atomConnect(listAtomLigand , serialAtom3)
    for atom in listAtomConnect3 : 
        if atom["element"] == "O" : 
            listAtomConnectO, connectO = retrieveAtom.atomConnect(listAtomLigand, atom["serial"])
            if connectO != ["O", "C"] : 
                return 0
        
    for atom in listAtomConnect3 :
        if atom["element"] == "C" and atom["serial"] != serialPrevious and atom["serial"] != serialAtom3: 
            atomConnect4, connectAtom4 = retrieveAtom.atomConnect(listAtomLigand, atom["serial"])
            if connectAtom4 == ["C", "O", "C", "C"] or connectAtom4 == ["C", "C", "O", "C"] or connectAtom4 == ["C", "C", "C", "O"] : 
                if riboseFromADPRiboseAtom4 (atom["serial"], listAtomLigand, serialInit, serialAtom3) == 1 : 
                    return 1
    return 0
    
    
    
def riboseFromADPRiboseAtom4 (serial, listAtomLigand, serialInit, serialPrevious) : 
    """search ribose atom4 
    in: list atom connected of nitrogen, list atom ligand
    out: boolean"""    
    
    listAtomC, connectC = retrieveAtom.atomConnect(listAtomLigand, serial) 
    
    for atom in listAtomC[1:] : 
        if atom["element"] == "O" : 
            listAtomConnectO, connectO = retrieveAtom.atomConnect(listAtomLigand, atom["serial"])
            if connectO != ["O", "C"] : 
                return 0

    for atom in listAtomC :
        if atom["serial"] != serial and atom["serial"] != serialPrevious :  
            for serial in atom["connect"] : 
                if serial == serialInit : 
                    return 1
    
    return 0
    
    
    
##########Search Imidazole###############

def imidazole(listAtomConnectNitrogen, listAtomLigand):
    """search imidazole global
    in: list atom connected of nitrogen, list atom ligand
    out: boolean"""

    amine = toolSubstructure.matrixElement(listAtomConnectNitrogen)
    
    if amine == ["N", "C", "C"]:
        groupAtomC1, conect_C1 = retrieveAtom.atomConnect(listAtomLigand, int (listAtomConnectNitrogen[0]["connect"][1]))
        groupAtomC2, conect_C2 = retrieveAtom.atomConnect(listAtomLigand, int (listAtomConnectNitrogen[0]["connect"][2]))
        if conect_C1 == ["C", "C", "N", "N"] or conect_C1 == ["C", "N", "N", "C"] or conect_C1 == ["C", "N", "C", "N"] or conect_C1 == ["C", "N", "N"]:
            if imidazoleATOM3(groupAtomC1, listAtomConnectNitrogen[0]["serial"], listAtomLigand) == 1:
                return 1

        if conect_C2 == ["C", "C", "N", "N"] or conect_C2 == ["C", "N", "N", "C"] or conect_C2 == ["C", "N", "C", "N"] or conect_C2 == ["C", "N", "N"]:
            if imidazoleATOM3(groupAtomC2, listAtomConnectNitrogen[0]["serial"], listAtomLigand) == 1:
                return 1

    return 0


def imidazoleATOM3(listAtom, serialAtomInit, listAtomLigand):
    """Check the atom 3 in cercle of imidazole
    in : atoms ligand in list, serial of first nitrogen, name of ligand
    out : boolean"""

    for atom in listAtom:
        if atom["element"] == "N":
            if "conect_N1" in locals():
                groupAtomN2, conect_N2 = retrieveAtom.atomConnect(listAtomLigand, atom["serial"])
            else:
                groupAtomN1, conect_N1 = retrieveAtom.atomConnect(listAtomLigand, atom["serial"])

    if conect_N1 == ['N', 'C', 'C']:
        if imidazoleATOM4(groupAtomN1, serialAtomInit, listAtomLigand) == 1:
            return 1

    if conect_N2 == ['serialAtomInit', 'C', 'C']:
        if imidazoleATOM4(groupAtomN2, serialAtomInit, listAtomLigand) == 1:
            return 1

    return 0


def imidazoleATOM4(listAtomNitrogen, serialAtomInit, listAtomLigand):
    """Check the atom 4 in cercle of imidazole
    in : atoms ligand in list, serial of first nitrogen, name of ligand
    out : boolean"""
    
    for atom in listAtomNitrogen:
        if atom["element"] == "C":
            if "conect_C3" in locals():
                groupAtomC4, conect_C4 = retrieveAtom.atomConnect(listAtomLigand, atom["serial"])
            else:
                groupAtomC3, conect_C3 = retrieveAtom.atomConnect(listAtomLigand, atom["serial"])

    if  imidazoleATOM5( conect_C3, serialAtomInit, listAtomLigand) == 1:
        return 1
    if imidazoleATOM5( conect_C4, serialAtomInit, listAtomLigand) == 1:
        return 1

    return 0


def imidazoleATOM5(connectMatrixElement, serialAtomInit, listAtomLigand):
    """Check the atom 5 in cercle of imidazole
    in : atoms ligand in list, serial of first nitrogen, name of ligand
    out : boolean"""
    
    lengthGroup = len(listAtomLigand)

    for compound in connectMatrixElement:
        if compound != "N":
            if compound != "C":
                return 0

    for i in range(1, lengthGroup):
        if listAtomLigand[i]["element"] != "C":
            if listAtomLigand[i]["element"] != "N":
                return 0

        if listAtomLigand[i]["element"] == "C":
            if "conect_C5" in locals():
                groupAtomC6, conect_C6 = retrieveAtom.atomConnect(listAtomLigand, listAtomLigand[i]["serial"])
            else:
                groupAtomC5, conect_C5 = retrieveAtom.atomConnect(listAtomLigand, listAtomLigand[i]["serial"])

    if not "conect_C5" in locals():
        return 0
    else:
        for atomC5 in groupAtomC5:
            if atomC5["serial"] == serialAtomInit:
                return 1

        if "groupAtomC6" in locals():
            for atomC6 in groupAtomC6:
                if atomC6["serial"] == serialAtomInit:
                    return 1
    return 0

######################################################################################################################


def interestGroup (distanceMax, listAtomLigand, namePDB, angleOption):
    """Search different amine
    in : ligands in namePDB
    out : major nitrogen in the group of different structures"""
    
    amine = structure.amine()
    listSerialNitrogen = listNitrogen(listAtomLigand)
    atomImidazole = []
    atomGuanidium = []
    atomDiamine = []
    for serialNitrogen in listSerialNitrogen:
        listAtomConnectNitrogen, conect = retrieveAtom.atomConnect(listAtomLigand, serialNitrogen)
        if imidazole(listAtomLigand, listAtomConnectNitrogen[0]) == 1:
            atom = buildAtom(distanceMax, listAtomConnectNitrogen[0], namePDB, "Imidazole", listAtomLigand)
            if angleOption == 1 : 
                checkAngleInSearchNeighbor(atom, "Imidazole")
            if len(atom["neighbors"]) != 0 : 
                atomImidazole.append(atom)
        elif guanidium(listAtomConnectNitrogen, listAtomLigand) == 1:
            atom = buildAtom(distanceMax, listAtomConnectNitrogen[0], namePDB, "Guanidium", listAtomLigand)
            if angleOption == 1 : 
                checkAngleInSearchNeighbor(atom, "Guanidium")
            if len(atom["neighbors"]) != 0 : 
                atomGuanidium.append(atom)
        elif guanidiumNitrogenCentral(listAtomConnectNitrogen, listAtomLigand) == 1:
            atom = buildAtom(distanceMax, listAtomConnectNitrogen[0], namePDB, "Guanidium", listAtomLigand)
            if angleOption == 1 : 
                checkAngleInSearchNeighbor(atom, "Guanidium")
            if len(atom["neighbors"]) != 0 : 
                atomGuanidium.append(atom)
        elif diAmine(listAtomConnectNitrogen, listAtomLigand) == 1:
            atom = buildAtom(distanceMax, listAtomConnectNitrogen[0], namePDB, "Diamine", listAtomLigand)
            if angleOption == 1 : 
                checkAngleInSearchNeighbor(atom, "Diamine")
            if len(atom["neighbors"]) != 0 : 
                atomDiamine.append(atom)
        elif pyridine(listAtomConnectNitrogen, listAtomLigand) == 1:
            atom = buildAtom(distanceMax, listAtomConnectNitrogen[0], namePDB, "Pyridine", listAtomLigand)
            if angleOption == 1 :
                checkAngleInSearchNeighbor(atom, "Pyridine")
            if len(atom["neighbors"]) != 0 : 
                amine["Pyridine"].append(atom)
        elif cn(listAtomConnectNitrogen, listAtomLigand) == 1:
            atom = buildAtom(distanceMax, listAtomConnectNitrogen[0], namePDB, "Primary", listAtomLigand)
            if angleOption == 1 :
                checkAngleInSearchNeighbor(atom, "Primary")
            if len(atom["neighbors"]) != 0 : 
                amine["Primary"].append(atom)
        elif cnc(listAtomConnectNitrogen, listAtomLigand) == 1:
            atom = buildAtom(distanceMax, listAtomConnectNitrogen[0], namePDB, "Secondary", listAtomLigand)
            if angleOption == 1 :
                checkAngleInSearchNeighbor(atom, "Secondary")
            if len(atom["neighbors"]) != 0 : 
                amine["Secondary"].append(atom)
        elif cncc(listAtomConnectNitrogen, listAtomLigand) == 1:
            atom = buildAtom(distanceMax, listAtomConnectNitrogen[0], namePDB, "Tertiary", listAtomLigand)
            if angleOption == 1 :
                checkAngleInSearchNeighbor(atom, "Tertiary")
            if len(atom["neighbors"]) != 0 : 
                amine["Tertiary"].append(atom)
    
    
    if len(atomGuanidium) > 0 : 
        regroupAtomNeighborGuanidium(atomGuanidium, amine["Guanidium"], listAtomLigand)
    if len(atomImidazole) > 0 : 
        regroupAtomNeighbor(atomImidazole, amine["Imidazole"], listAtomLigand)
    if len(atomDiamine) > 0 : 
        regroupAtomNeighbor(atomDiamine, amine["Diamine"], listAtomLigand)
    
    return amine

#######regroup neighbors case of imidazole, guanidium and diamine###########

def regroupAtomNeighborGuanidium(listAtomRegroup, amine, listAtomLigand): ##A revoir tres tres lourds en temps 
    """Regroup neighbors for 3 nitrogen atoms guanidium
    in: list atom finds in search neighbors, count structure amine, list atom of ligand in pdb
    out: modification amine structure"""
    
    
    for atomRegroup in listAtomRegroup : 
        atomRegroupConnect, connectAtomRegroup = retrieveAtom.atomConnect(listAtomLigand, atomRegroup["serial"])
        if connectAtomRegroup == ["N", "C"] : 
            atom1 = atomRegroup 
            
    outAtom = {}
    outAtom["PDB"] = atom1["PDB"]
    outAtom["resName"] = atom1["resName"]
    outAtom["serial"] = atom1["serial"]
    outAtom["neighbors"] = deepcopy(atom1["neighbors"]) 
    
    listAtomConnectN, connectMatrixN = retrieveAtom.atomConnect(listAtomLigand, atom1["serial"])
    listAtomConnectC, connectMatrixC = retrieveAtom.atomConnect(listAtomLigand, listAtomConnectN[1]["serial"])
    
    for atomConnectC in listAtomConnectC : 
        for atom in listAtomRegroup :
            if atom["serial"] == atomConnectC["serial"] : 
                for neighbor in atom["neighbors"] : 
                    if not neighbor in outAtom["neighbors"] : 
                        outAtom["neighbors"].append(neighbor)
    
    amine.append(outAtom)
    



def regroupAtomNeighbor (atomRegroup, amine, listAtomLigand):
    """regroup equivalent atom
    in: atom no regroup, amine structure neighbor, list atom ligand
    out: amine modified"""
    
    if len(atomRegroup) == 2 : 
        regroup = regroupNeighbor(atomRegroup[0]["serial"], atomRegroup[1]["serial"], atomRegroup)
        amine.append(regroup)
    else : 
        listatom = []
        for atom in atomRegroup : 
            listatom.append(atom["serial"])
        nbSerial = len(listatom)
      
        i = 0
        while i < nbSerial :  
            flag = 0
            serialSearch = listatom[i]
            groupAtom, conect = retrieveAtom.atomConnect(listAtomLigand, listatom[i])
            for atomGroup in groupAtom :
                if flag == 1 : 
                    break
                if serialSearch != atomGroup["serial"] and serialSearch in atomGroup["connect"] : 
                    for serialConnect in atomGroup["connect"] : 
                        if serialConnect in listatom and serialConnect != serialSearch : 
                            regroup = regroupNeighbor(serialSearch, serialConnect, atomRegroup)
                            amine.append(regroup)
                            del listatom[i]
                            nbSerial = nbSerial - 1
                            flag = 1
                            break
            i = i + 1 
                    
                    
def regroupNeighbor(serial1, serial2, listAtom):
    """Regroup list neighbor that 2 atom
    in: atom serial 1, atom serial 2, list atoms
    out: list neighbors global"""
    
    for atom in listAtom : 
        if atom["serial"] == serial1 : 
            atom1 = atom
            
        if atom["serial"] == serial2 : 
            atom2 = atom

#    print len(atom1["neighbors"]), len(atom2["neighbors"])    
    outAtom = {}
    outAtom["PDB"] = atom1["PDB"]
    outAtom["resName"] = atom1["resName"]
    outAtom["serial"] = atom1["serial"]
    outAtom["neighbors"] = deepcopy(atom1["neighbors"])
    
    for neighbor in atom2["neighbors"] : 
        outAtom["neighbors"].append(neighbor)
    
#    print len(outAtom["neighbors"])
    return outAtom
    

######################################################################################


def buildAtom(rayon, serialAtomNitrogen, namePDB, typeStudy, listAtomLigand):
    """Build count structure for amine structure
    in: distance, serial atom nitrogen, type study, list atom ligand
    out: count strucutre"""
    
    atom = {}
    atom["PDB"] = namePDB
    atom["resName"] = serialAtomNitrogen["resName"]
    atom["serial"] = serialAtomNitrogen["serial"]
    atom["neighbors"] = neighbors(rayon, serialAtomNitrogen, namePDB, typeStudy, listAtomLigand)

    return atom


def buildAtomGlobal(rayon, atomFind, pdb):
    """Retrieve information atom find
    in : atom find, name of PDB
    out : structure of atoms
    """
    
    atom = {}
    atom["PDB"] = pdb
    atom["resName"] = atomFind["resName"]
    atom["serial"] = atomFind["serial"]
    atom["neighbors"] = neighborsGlobal(rayon, atomFind, pdb)

    return atom



def neighbors(rayon, atomN, pdb, typeStructure, ligandPDB):
    """Search neighbors for all ligand
    in : rayon where is atoms, central atom, pdb file
    out : list atoms found"""

    listAtom = []
    linesPDB = loadFile.openPdbFile(pdb)
    for line in linesPDB:
        if search("^ATOM", line) or search("^HETATM", line): 
            atom = parsing.lineCoords(line)
            if atom != {} and atom["element"] != "H":
                distance = calcul.distanceTwoatoms(atomN, atom)
                if distance <= rayon and distance != 0.0:
                    if atomN["resSeq"] != atom["resSeq"]:
                        if tool.atomInList(listAtom, atom) == 0:
                            atom["distance"] = distance
                            atom["angle"] = calcul.angle(atomN, atom, ligandPDB, typeStructure)
                            atom["classification"] = tool.classification(atom)
                            atom["classificationAtLeastOne"] = tool.atLeastOneClassification(atom)
                            listAtom.append(atom)

    return listAtom


def neighborsGlobal(rayon, atomFind, pdb):
    """Search neighbors for all ligand
    in : rayon where is atoms, central atom, pdb file
    out : list atoms found"""

    listAtom = []
    linesPDB = loadFile.openPdbFile(pdb)
    for line in linesPDB:
        if search("^ATOM", line) or search("^HETATM", line):
            atom = parsing.lineCoords(line)
            if atom != {} and atom["element"] != "H":
                distance = calcul.distanceTwoatoms(atomFind, atom)
                if distance <= rayon and distance != 0.0:
                    if atomFind["resSeq"] != atom["resSeq"]:
                        if tool.atomInList(listAtom, atom) == 0:
                            atom["distance"] = distance
                            atom["classification"] = tool.classification(atom)
                            atom["classificationAtLeastOne"] = tool.atLeastOneClassification(atom)
                            listAtom.append(atom)

    return listAtom


def repetitionInPBD(ligand, pdb):
    """Count number of same ligand in pdb file
    in : name of ligand, name of pdb file
    out : number of ligand in pdb"""

    groupAtom = loadFile.ligandInPDB(pdb, ligand)

    listLigand = []
    for atom in groupAtom:
        if not atom["resSeq"] in listLigand:
            listLigand.append(atom["resSeq"])

    return len(listLigand)


def globalNeighbors(distance, ligand, pdb):
    """Retrieve neighbor in sphere
    in : distance arroud  atom, atom in ligand, name of pdb
    out : list of atoms"""
    
    globalNeighbor = []

    for atomLigand in ligand:
        groupAtom, conect = retrieveAtom.atomConnect(ligand, atomLigand["serial"])
        atom = buildAtomGlobal(distance, groupAtom[0], pdb)
        globalNeighbor.append(atom)

    return globalNeighbor


def listNitrogen(groupAtom):
    """Search N in the atoms list
    return a list of number atom"""

    listN = []
    for i in groupAtom:
        if i["element"] == "N":
            serial = i["serial"]
            if not serial in listN:
                listN.append(i["serial"])
    return  listN



def checkAngleInSearchNeighbor(atomRetrieve, typeStructStudy):
    """Check for every neighbors of atom retrieve the angleVector
    in: atomRetrieve
    out: atomRetrieve modified"""
    
    conditionAngle = structure.selectionAngle()
    nbNeighbors = len(atomRetrieve["neighbors"])
        
    i = 0
    while i < nbNeighbors :
        
        for angleVector in atomRetrieve["neighbors"][i]["angleVector"] : 
            if angleVector < conditionAngle[typeStructStudy]["INF"] :
                del atomRetrieve["neighbors"][i]
                nbNeighbors = nbNeighbors - 1
                break
            elif angleVector > conditionAngle[typeStructStudy]["SUP"] : 
                del atomRetrieve["neighbors"][i]
                nbNeighbors = nbNeighbors - 1
                break
            
        i = i +1
                
    



##################################Serine Protease #############################################

#def cycle(atomInitial, firstPosition, atomLigand, nbTest, listRetrieve, serialTest):
#
#    #print nbTest
#    if atomInitial == 0 or len(atomInitial["connect"]) < 3 or atomInitial["serial"] == serialTest : 
##        print "OUT1"
#        del listRetrieve
#        return 0, []
#    else : 
#        appendAtomConnect(atomInitial, listRetrieve, atomLigand)
#    
#    if nbTest == 0 : 
#        if firstPosition in atomInitial["connect"] [1:] and firstPosition != serialTest:
#            appendAtomConnect(atomInitial, listRetrieve, atomLigand)
#            appendAtomConnect(retrieveAtom.serial(firstPosition, atomLigand), listRetrieve, atomLigand)
#            return 1 , listRetrieve
#        else : 
##            print "OUT2"
#            del listRetrieve
#            return 0 , []
#    else : 
#        #print "in"
#        #appendAtomConnect(atomInitial, listRetrieve, atomLigand)
#        connectMatrix = atomInitial["connect"]
#        for serialConnect in connectMatrix :
#            #print serialConnect
#            
#            if serialConnect != atomInitial["serial"] and serialConnect != serialTest:  
#                atomtest = retrieveAtom.serial(serialConnect, atomLigand)
#                result, listRetrieve = cycle(atomtest, firstPosition, atomLigand, nbTest - 1, listRetrieve, atomInitial["serial"])
#                if result == 1 : 
#                    return 1, listRetrieve
#    
#    
#        del listRetrieve            
#        return 0, []
#            
#
#
#
#def cycle2(atomInitial, firstPosition, atomLigand, nbTest, strSerialList, serialTest):
#
#    #print nbTest
#    if atomInitial == 0 or len(atomInitial["connect"]) < 3 or atomInitial["serial"] == serialTest : 
##        print "OUT1"
#        return 0, ""
#    else : 
#        strSerialList = strSerialList + "_" + str(atomInitial["serial"])
#    
#    if nbTest == 0 : 
#        if firstPosition in atomInitial["connect"] [1:] and firstPosition != serialTest:
#            return 1 , strSerialList
#        else : 
##            print "OUT2"
#            return 0 , ""
#    else : 
#        #print "in"
#        #appendAtomConnect(atomInitial, listRetrieve, atomLigand)
#        connectMatrix = atomInitial["connect"]
#        for serialConnect in connectMatrix :
#            #print serialConnect
#            if serialConnect != atomInitial["serial"] and serialConnect != serialTest:  
#                atomtest = retrieveAtom.serial(serialConnect, atomLigand)
#                result , strSerialList = cycle2(atomtest, firstPosition, atomLigand, nbTest - 1, strSerialList, atomInitial["serial"])
#                if result == 1 : 
#                    return 1, strSerialList
#    
#    
#        return 0, ""


        
def appendAtomConnect(atom, listRetrive, atomLigand):
    
    listSerial = []
    for atomRetrieve in listRetrive :
        listSerial.append(atomRetrieve["serial"])
    
    for serialAtom in atom["connect"] :
        if not serialAtom in listSerial : 
            atom = retrieveAtom.serial(serialAtom, atomLigand)
            if atom != 0 :
                listRetrive.append(atom)


#def checkCycle (atom, firstPosition, lastPosition, n, atomLigand):
#    
#    if firstPosition in atom["connect"][:1] : 
#        return n
#    
#    else : 
#        for atSerial in atom["connect"][:1] : 
#            if atSerial != lastPosition : 
#                return checkCycle(retrieveAtom.serial(atSerial, atomLigand), firstPosition, atSerial, n + 1 , atomLigand)
#    
#    return 0
#    
    
def cycleGlobal(atomLigand):
    
    for atom in atomLigand : 
        if cycle3(atom["serial"], atom["serial"], atom["serial"], atomLigand, 5) == 1 : 
            atom["cycle"] = 1 
        else : 
            atom["cycle"] = 0
    
    
def cycle3(serialFirst, serialTest, serialPrevious, atomLigand, test) :
    atomTest = retrieveAtom.serial(serialTest, atomLigand)
    if atomTest == 0 : 
        return 0
    
    if serialFirst in atomTest["connect"] and  serialFirst != serialTest and serialFirst != serialPrevious : 
        return 1
    
    if test == 0 : 
        return 0
    
    for serialConnect in atomTest["connect"][1:] : 
        if serialConnect != serialPrevious : 
            if cycle3(serialFirst, serialConnect, serialTest, atomLigand, test - 1) == 1 : 
                return 1
            
    return 0
    
    
    
def cycleOnlyTestCarbon(serialFirst, serialTest, serialPrevious, atomLigand, test, flagNitrogen) :
    atomTest = retrieveAtom.serial(serialTest, atomLigand)
    
    if serialFirst == atomTest["serial"] and test == 0: 
        return 1
    
    if atomTest == 0 : 
        return 0
    
    if atomTest ["element"] != "C" and test != 0 and serialFirst != serialTest:
        return 0
    

    for conect in atomTest["connect"] : 
        atomSecond = retrieveAtom.serial(conect , atomLigand)
        if atomSecond == 0 : 
            return 0
        
        elif  atomSecond["element"] == "N" : 
            flagNitrogen = flagNitrogen + 1
            if flagNitrogen > 3 :
                return 0
            else : 
                continue
            
        elif  atomSecond["element"] != "C" : 
            return 0    
         
    if test == 0 : 
        return 0
    
    for serialConnect in atomTest["connect"][1:] : 
        if serialConnect != serialPrevious : 
            if cycleOnlyTestCarbon(serialFirst, serialConnect, serialTest, atomLigand, test - 1, flagNitrogen) == 1 : 
                return 1

    return 0
