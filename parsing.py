import formatCharacter
import loadFile
from re import search
from re import sub


def lineCoords (line):
    """Parsing line of coordinate PDB File
    in: line
    out: dictionnary atom"""

    atom = {}
    atom["serial"] = formatCharacter.formatInt(line[6:11])
    atom["name"] = formatCharacter.suppSpace(line[12:16])
    atom["char"] = line[16]
    atom["resName"] = formatCharacter.suppSpace(line[17:20])
    atom["chainID"] = str(line[21])
    atom["resSeq"] = formatCharacter.formatInt(line[22:26])
    atom["iCode"] = str(line[26])
    atom["x"] = formatCharacter.formatFloat (line[30:38])
    atom["y"] = formatCharacter.formatFloat (line[38:46])
    atom["z"] = formatCharacter.formatFloat (line[46:54])
    atom["element"] = formatCharacter.suppSpace(line[76:78])
    atom["charge"] = formatCharacter.suppSpace(line[78:80])
    atom["occupancy"] = formatCharacter.suppSpace(line[54:60])
    atom["tempFactor"] = formatCharacter.suppSpace(line[60:66])
    
    atom["connect"] = []
    return atom

def lineConnectMatrix(line):
    """Retrieve connect matrix by atom
    in: line
    out: list of connect"""

    connectAtom = []
    try: connectAtom.append(formatCharacter.formatInt(line[6:11]))
    except: return connectAtom
    try: connectAtom.append(formatCharacter.formatInt(line[11:16]))
    except: return connectAtom
    try: connectAtom.append(formatCharacter.formatInt(line[16:21]))
    except: return connectAtom
    try: connectAtom.append(formatCharacter.formatInt(line[21:26]))
    except: return connectAtom
    try: connectAtom.append(formatCharacter.formatInt(line[26:31]))
    except: return connectAtom

    return connectAtom

def methods(namePDB):
    """Retrieve by PDB file the method of study / RMN or X-ray
    in : name PDB file
    out : methods -> format string lower"""


    fileLines = loadFile.openPdbFile(namePDB)
    type = fileLines[0][10:50]
    type = formatCharacter.suppSpace(type)

    return type.lower()


def resolution(namePDB):
    """Retrieve by PDB file the resolution if X-ray methods
    in : name of pdb file
    out : resolution -> format float, return 1000.00 if have not resolution in file"""

    fileLines = loadFile.openPdbFile(namePDB)

    for line in fileLines:
        if search("^REMARK   2 RESOLUTION", line):
            line = sub('[ ]{2,}', ' ', line)
            try:
                resolution = formatCharacter.formatFloat(line.split(" ")[3])
            except:
                resolution = 1000.0
            return resolution
    return 1000.0
