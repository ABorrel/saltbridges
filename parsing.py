import formatCharacter
import loadFile
import calcul
import pathManage

from re import search, sub
from copy import deepcopy
from os import path

def lineCoords (line, remove_H = 0):
    """Parsing line of coordinate PDB File
    in: line
    out: dictionnary atom"""

    line = line.strip ()
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
    if atom["element"] == "H" and remove_H == 1 :
        return None 
    return atom


def EmptyAtom ():


    atom = {}
    atom["serial"] = 0
    atom["name"] = ""
    atom["char"] = ""
    atom["resName"] = ""
    atom["chainID"] = ""
    atom["resSeq"] = 0
    atom["iCode"] = "0"
    atom["x"] = 0.0
    atom["y"] = 0.0
    atom["z"] = 0.0
    atom["element"] = "Z"
    atom["charge"] = "0"
    atom["occupancy"] = "0"
    atom["tempFactor"] = "0"
    
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

def header(namePDB):
    """Retrieve for one PDB file the header
    in : name PDB file
    out : header -> format string lower"""

    l_lines = loadFile.openPdbFile(namePDB)
    return l_lines[0][6:].lower().strip ()


def methodStructure (namePDB):
    l_lines = loadFile.openPdbFile(namePDB)
    
    for l in l_lines : 
        if search("^EXPDTA", l) and search("SOLUTION", l) :
            return "SOLUTION"
    return "XRAY" 
            


def Quality(namePDB):
    """Retrieve by PDB file the Quality if X-ray header
    in : name of pdb file
    out : Quality -> format float, return 1000.00 if have not Quality in file"""

    l_lines = loadFile.openPdbFile(namePDB)
    nb_line = len (l_lines)

    i = 0
    while i < nb_line :
        if "rfactor" in locals ().keys () and "resolution" in locals ().keys () :
            return [resolution, rfactor]
        if search("^REMARK   2 RESOLUTION", l_lines[i]):
            line = sub('[ ]{2,}', ' ', l_lines[i])
            try:
                resolution = formatCharacter.formatFloat(line.split(" ")[3])
            except:
                resolution = 1000.0
        
        elif search ("REMARK   3   R VALUE", l_lines[i]) : 
            rfactor = l_lines[i].strip ().split (":")[-1].replace (" ", "")
            if rfactor == "NULL" : 
                rfactor = 0.0
            else : 
                rfactor =  float (rfactor)
        i = i + 1
    
    return [1000.0, 100.0]




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
    



def loadCoordSectionPDB (p_PDBin, section = "", remove_H = 0, debug = 0):
    """
    Retrieve every atom in coordinate section. If it is NMR complex
    retrieve only first model
    """
    
    l_atom = []
    try : filin = open (p_PDBin, "r")
    except : return []
    list_line_PDB = filin.readlines()
    filin.close ()
    
    for line_PDB in list_line_PDB :
        #End model
        if search ("^ENDMDL", line_PDB) : 
            break
        
        if section == "" : 
            if search ("^ATOM", line_PDB) or search ("^HETATM", line_PDB) : 
                atom = lineCoords(line_PDB, remove_H)
                if atom != None : 
                    l_atom.append (atom)
        else : 
            if search ("^" + section, line_PDB)  : 
                atom = lineCoords(line_PDB, remove_H)
                if atom != None : 
                    l_atom.append (atom)
                
#    if debug : 
#        print "TEST"            
#        print len (l_atom)
#    
    return l_atom
    


def retrieveLigand  (l_atom, name_ligand, extend = 2, debug = 0):
    """
    Retrieve list of ligand in PDB
    args: -> PDB parsed
          -> name ligand
    return: list of ligand with atoms
    """
    
    if debug : 
        print "-control ligand select-"
        print name_ligand, "NAME ligand"
        print "l144 - parse PDB"
        print "****************"
        
    name_ligand = name_ligand.upper ()
    list_atom_ligand = []
    for element in l_atom : 
        if element ["resName"] == name_ligand : 
            list_atom_ligand.append (deepcopy(element))
            
    for atomLigand in list_atom_ligand:
        atomLigand["connect"].append(atomLigand["serial"])
        for atom_ligand_connect in list_atom_ligand:
            distance = calcul.distanceTwoatoms(atomLigand, atom_ligand_connect)
            if distance < extend and distance != 0:
                if not atom_ligand_connect["serial"] in atomLigand["connect"]:
                    atomLigand["connect"].append(atom_ligand_connect["serial"])
    
    # Check exotic bound (SE-C...) with special length of bound
    if checkConnectMatrix(list_atom_ligand) == 0 :
        if debug : 
            print "- exotic ligand -, recursive distance"
            print "l165 - parse PDB"
            print "Threshold: ", extend
            print "--------------------"
        if extend <= 4.0 : 
            return retrieveLigand (l_atom, name_ligand, extend + 0.1, debug = 0)
    
    return separateByLigand (list_atom_ligand)



def checkConnectMatrix (list_atom_ligand):
    
    list_nb_atom = []
    list_atom_serial_by_ligand = retrieveListAtomID (list_atom_ligand)
    for atom_serial in list_atom_serial_by_ligand : 
        list_nb_atom.append(len(atom_serial)), "len of ligand"
        
    list_nb_atom = list(set(list_nb_atom))
    if len(list_nb_atom) != 1 : 
        return 0
    return 1


    
def retrieveListAtomID (list_atom):
    """
    Retrieve list ID atoms
    args: -> list atoms
    return: list serial
    """
    
    list_out = []
    list_atom_temp = deepcopy(list_atom)
    #except : return []
    
    while len (list_atom_temp) != 0 : 
        list_out.append (retrieveSerialLigand(list_atom_temp))
    
    return list_out
    
    
def retrieveSerialLigand (list_atom):
    """
    Retrieve serial atom by ligand
    args: -> list atoms
    return: -> list serial atoms
    """
    
    list_out = list_atom[0]["connect"]
    del list_atom[0]
    nb_atom = len (list_atom)
    
    validate = 100
    while validate != 0 :
        nb_atom_temp = nb_atom 
        i = 0
        while i < nb_atom :
            if list_atom[i]["serial"] in list_out :
                list_out = list_out + list_atom[i]["connect"]
                del list_atom[i]
                nb_atom = nb_atom - 1
            else :
                i = i + 1
        validate = nb_atom - nb_atom_temp 
    return list(set(list_out))
    




def checkLigandHooked (l_atom_complex, l_atom_lig, thresold = 1.65):
    
    l_atom_prot = deepcopy(l_atom_complex)
    
    for atom_lig in l_atom_lig : 
        i = 0
        nb_atom = len(l_atom_prot)
        while i < nb_atom : 
            d = calcul.distanceTwoatoms(l_atom_prot[i], atom_lig)
            if d < thresold : 
                print d
                return 1
            # reduce combination with atom away
            elif d > (thresold * len (l_atom_lig)) : 
                del l_atom_prot[i]
                nb_atom = nb_atom - 1

            i = i + 1
            
    return 0




def separateByLigand (l_atom_ligand, debug = 0) :
    """
    Separate list atoms ligand with same name or ID by atomic position
    args: -> list atoms ligand
    return: -> list of list with several atom by ligand
    """
    if debug : 
        print "--- check ligand ---"
        print  "NB atoms: ", len (l_atom_ligand)
        print "First atoms: ", l_atom_ligand[0]
        print "l209 parsePDB"
        print "--------------------"
    
    ####################################
    # ligand separed by ID and num res #
    ####################################
    if l_atom_ligand == [] :
        return []
    
    
    l_code_ligand = []
    for atom_ligand in l_atom_ligand : 
        code_ligand = atom_ligand["chainID"] + "_" + str(atom_ligand["resSeq"])
        if not code_ligand in l_code_ligand : 
            l_code_ligand.append (code_ligand)
    
    if len (l_code_ligand) == 1 : 
        return [l_atom_ligand]
    
    else : 
        # control number of atom but eg DEX case where same ligand without number of atoms
#         if len (l_atom_ligand) % len (l_code_ligand) != 0 : 
#             print "- ligand several ligand without same number of atoms"
#             print "l229 parsePDB"
#             print "----------------"
        
        l_out = []
        # separe by ligand, list of list
        for code_ligand in l_code_ligand : 
            l_atom_ligand_out = []
            for atom_ligand in l_atom_ligand : 
                code_atom = atom_ligand["chainID"] + "_" + str(atom_ligand["resSeq"])
                if code_atom == code_ligand : 
                    l_atom_ligand_out.append (atom_ligand)
            l_out.append (l_atom_ligand_out)
        
        return l_out




def retrieveListLigand ( p_file_PDB, debug = 0 ):
    """
    Retrieve list name_lig in PDB file, only ligands, remove metals
    args: - path file
    return: list ligands
    """
    
    if not path.exists(p_file_PDB) : 
        return []

    filin = open ( p_file_PDB, "r" )
    l_line_pdb = filin.readlines ()
    filin.close()
    
    l_out = []
    for line_pdb in l_line_pdb : 
        if search ( "^HETATM", line_pdb ) : 
            name_lig = line_pdb[17:21].replace ( " ", "" )
            if debug : print name_lig
            if not name_lig in l_out and name_lig != "HOH"  : # remove water
                if len (name_lig) > 2 : # remove metals just 2 letters
                    l_out.append ( name_lig )
                    
    return l_out 
    


def BuildDicoRes(l_atoms) :
    d_res = {}
    for atom in l_atoms :
        
        res_ID = atom["resSeq"]
        res_name = atom["resName"]
        res_chain = atom["chainID"]
        
        k_in = str(res_name) + "_" + str(res_ID) + "_" + str (res_chain)
        
        if not k_in in d_res.keys () :
            d_res[k_in] = []
        d_res[k_in].append (atom)

    return d_res


def BuildByRes(PDBid):

    rep = pathManage.openPdbFile()
    pfilin = rep + PDBid + ".pdb"


    if not path.exists(pfilin):
        print "ERROR l-417 parsing"
        return {}
    else:
        filin = open(pfilin, "r")
        llines = filin.readlines()
        filin.close()

    dres = {}
    for linefile in llines:
        if search("^ATOM", linefile) or search("^HETATM", linefile):
            atom = lineCoords(linefile)
            if atom != {} and atom["element"] != "H":
                res_ID = atom["resSeq"]
                res_name = atom["resName"]
                res_chain = atom["chainID"]
                # key by residues
                k_in = str(res_name) + "_" + str(res_ID) + "_" + str (res_chain)
                if not k_in in dres.keys():
                    dres[k_in] = []
                dres[k_in].append(atom)
    return dres




