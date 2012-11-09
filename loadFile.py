from re import search
from re import sub

import calcul
import formatCharacter
import parsing
import repertory



def openPdbFile(namePDB):
    """Open PDB file, check if there are in file many model structure and chose first model
        in: name file
        out: lines interest in list"""

    rep = repertory.openPdbFile()
    filin = open(rep + namePDB + ".pdb")
    linesFile = filin.readlines()
    filin.close()

    nbLines = len(linesFile)
    for i in range(0, nbLines):
        if search("^ENDMDL", linesFile[i]):
            lineEnd = i
        if search("^CONECT", linesFile[i]):
            lineStartConect = i
            break

    if "lineStartConect" in locals() and "lineEnd" in locals():
        return linesFile[0: lineEnd] + linesFile[lineStartConect:nbLines]
    else:
        return linesFile



def resultFilterPDBLigand (nameFile):
    """load result file with ligand and PDB ID associated
    in: name file result
    out: list of ligand"""

    file_dataset = repertory.result() + nameFile
    file_open = open(file_dataset, "r")
    lines = file_open.readlines()
    file_open.close()
    len_file = len(lines)

    ligands = []
    for i in range(0, len_file):
        ligand = {}
        line = lines[i].split("\t")
        ligand["name"] = line[0]
        list_pdb = line[1].split("\n")[0]
        list_pdb = list_pdb.split(" ")
        ligand["PDB"] = []

        for pdb in list_pdb:
            if pdb != "":
                ligand["PDB"].append(pdb)

        ligands.append(ligand)

    return ligands




def ligandInPDB(namePDB, nameLigands):
    """load ligand in structure in PDB file on single ligand by pdb
    out : atom in ligands
    out: list of atoms that ligand"""

    linesPDB = openPdbFile(namePDB)
    ligand = []
    for line in linesPDB:
        if search ("^HETATM", line):
            atom = parsing.lineCoords(line)
            if atom["resName"] == nameLigands and atom["element"] != "H":
                ligand.append(atom)
    checkOnlyOneLigand(ligand)
    connectMatrix = connectMatrixInPDB(namePDB)

    if connectAtom(connectMatrix, ligand) == "False":
        print "Construt Matrix"
        calcul.buildConnectMatrix(ligand, namePDB)

    return ligand


def checkOnlyOneLigand(groupAtom):
    """Check if there is in the list of atoms only one ligand
    in: list of atom ligand
    out: NULL -> remove directly in list atom"""

    serialLigand = groupAtom[0]["resSeq"]
    chainID = groupAtom[0]["chainID"]

    nbAtom = len(groupAtom)

    i = 0
    while i < nbAtom:
        if serialLigand != groupAtom[i]["resSeq"] :
            del groupAtom[i]
            nbAtom = nbAtom - 1
        elif chainID != groupAtom[i]["chainID"] : 
            del groupAtom[i]
            nbAtom = nbAtom - 1
        else:
            i = i + 1


def ligandInPDBConnectMatrixLigand(pdbName, nameLigands):
    """load ligand in structure in PDB file
    out : list atom in ligands with connect matrix calculated"""

    linesPDB = openPdbFile(pdbName)
    ligand = []
    listSerial = []
    for line in linesPDB:
        if search ("^HETATM", line):
            atom = parsing.lineCoords(line)
            if atom["resName"] == nameLigands and atom["element"] != "H":
                if not atom["serial"] in listSerial : 
                    ligand.append(atom)
                    listSerial.append(atom["serial"])
    checkOnlyOneLigand(ligand)####retrieve only first ligand
    calcul.buildConnectMatrix(ligand, pdbName)

    return ligand



def connectMatrixInPDB(namePDB):
    """Retrieve connect matrix in PDB file
    in : name of PDB file
    out : Connect matrix
    """

    linesPDB = openPdbFile(namePDB)
    connectMatrix = []
    for line in linesPDB:
        if search("^CONECT", line):
            connectAtom = parsing.lineConnectMatrix(line)
            connectMatrix.append(connectAtom)

    return connectMatrix


def connectAtom(matrixConnect, listAtom):
    """Check if matrix connect is present
    in : - matrix connect
         - group atom
    out : false or true"""

    for atom in listAtom:
        for connect in matrixConnect:
            if connect[0] == atom["serial"]:
                atom["connect"] = connect

    for atom in listAtom:
        if atom["connect"] == []:
            return "False"

    return "True"




def resultLigandPDB(nameFile):
    """Load file result, PDB associated at ligand
    in: Name of file
    out: list ligand with PDB associated"""


    fileOpen = open(nameFile, "r")

    lineFile = fileOpen.readlines()
    outFile = {}

    for line in lineFile:
        line = line.split("\t")
        PDB = line[0]

        listLigand = line[1]
        listLigand = listLigand.split("\n")[0]
        listLigand = listLigand.split(" ")
        if listLigand == []:
            continue

        for ligand in listLigand:
            if ligand != "":
                try:
                    outFile[ligand].append(PDB)
                except:
                    outFile[ligand] = []
                    outFile[ligand].append(PDB)

    return outFile



def resolution (pdbfile):
    """Retrieve by PDB file the resolution if X-ray structure
    in : pdb file
    out : resolution -> format float, return 1000.00 if do not have resolution in file"""

    fileLines = openPdbFile(pdbfile)

    for line in fileLines:
        if search("^REMARK   2 RESOLUTION", line):
            line = sub('[ ]{2,}', ' ', line)

            try: resolution = formatCharacter.formatFloat(line.split(" ")[3])
            except: resolution = 1000.0

            return resolution

    return 1000.0


#def typeData(file):
#
#    fileLines = openPdbFile(file)
#    type = "NONE"
#
#    for line in fileLines:
#        if search("^EXPDTA", line):
#            line = sub('[ ]{2,}', ' ', line)
#            try:
#                type = line.split(" ")[1]
#            except:
#                pass
#
#    return type


def globalPDB(PDB, ligand):
    """Retrieve every lines coordinates in PDB file
    in: name PDB
    out: list with line PDB  """
    
    file = openPdbFile(PDB)
    out = []
    for line in file:
        if search("^ATOM", line) or search("^HETATM", line):
            atom = parsing.lineCoords(line)
            out.append(atom)

    connectMatrix = connectMatrixInPDB(PDB)
    print connectMatrix
    connectAtom(connectMatrix, out)
        #calcul.buildConnectMatrix(ligand, PDB)

    return out
