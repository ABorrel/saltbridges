import formatCharacter
import loadFile
from re import search
from re import sub


def lineCoords (line):
    """Parsing line of coordinate PDB File
    in: line
    out: dictionnary atom"""

    atom = {}
    try :atom["serial"] = int(line[6:11].replace (" ", ""))
    except :line[6:11].replace (" ", "")
    atom["name"] = line[12:16].replace (" ", "")
    atom["char"] = line[16]
    atom["resName"] = line[17:20].replace (" ", "")
    atom["chainID"] = str(line[21])
    atom["resSeq"] = int (line[22:26].replace (" ", ""))
    atom["iCode"] = str(line[26])
    atom["x"] = float (line[30:38].replace (" ", ""))
    atom["y"] = float (line[38:46].replace (" ", ""))
    atom["z"] = float (line[46:54].replace (" ", ""))
    atom["element"] = line[76:78].replace (" ", "")
    # pqr without element
    if atom["element"] == "" :
        if type (atom["name"][0]) is int :
            atom["element"] = atom["name"][1]
        else :
            atom["element"] = atom["name"][0] 
    
    atom["charge"] = line[78:80].replace (" ", "")
    atom["occupancy"] = line[54:60].replace (" ", "")
    atom["tempFactor"] = line[60:66].replace (" ", "")
    
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



def countH2O (path_file_PDB) : 
    
    count_H20 = 0
    filin = open (path_file_PDB, "r")
    list_lines = filin.readlines()
    filin.close ()
    
    for line_PDB in list_lines : 
        if search("^HETATM", line_PDB) : 
            atom_parsed = lineCoords(line_PDB)
            if atom_parsed["resName"] == "HOH" : 
                count_H20 = count_H20 + 1
    return count_H20
    


def loadCoordSectionPDB (path_PDB_file, section = "", debug = 1):
    """
    Retrieve every atom in cordiante section. If it is NMR complex
    retrieve only first model
    
    """
    
    list_atom = []
    filin = open (path_PDB_file, "r")
    list_line_PDB = filin.readlines()
    filin.close ()
    
    for line_PDB in list_line_PDB :
        #End model
        if search ("^ENDMDL", line_PDB) : 
            break
        
        if section == "" : 
            if search ("^ATOM", line_PDB) or search ("^HETATM", line_PDB) : 
                list_atom.append (lineCoords(line_PDB))
        else : 
            if search ("^" + section, line_PDB)  : 
                list_atom.append (lineCoords(line_PDB))
                
#    if debug : 
#        print "TEST"            
#        print len (list_atom)
#    
    return list_atom
    


