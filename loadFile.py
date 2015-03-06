from re import search, sub
from os import path, listdir

import calcul
import formatCharacter
import parsing
import pathManage
import structure
import checkPDBfile



def openPdbFile(namePDB):
    """Open PDB file, check if there are in file many model structure and chose first model
        in: name file
        out: lines interest in list"""

    rep = pathManage.openPdbFile()
    filin = open(rep + namePDB + ".pdb")
    list_lines = filin.readlines()
    filin.close()

    int_nblines = len(list_lines)
    for i in range(0, int_nblines):
        # retrieve only first model
        if search("^ENDMDL", list_lines[i]):
            lineEnd = i
        if search("^CONECT", list_lines[i]):
            lineStartConect = i
            break

    if "lineStartConect" in locals() and "lineEnd" in locals():
        return list_lines[0: lineEnd] + list_lines[lineStartConect:int_nblines] # concatene lists 
    else:
        return list_lines



def resultFilterPDBLigand (path_file):
    """load result file with ligand and PDB ID associated
    in: name file result
    out: list of ligand"""

    print path.getsize(path_file)
    # list empty
    if path.getsize(path_file) == 0 : 
        return []

    file_open = open(path_file, "r")
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




def ligandInPDB(PDB_ID, ligand_ID):
    """load list_atom_ligand in structure in PDB file on single list_atom_ligand by pdb and remove H
    out : atom in ligands
    out: list of atoms that list_atom_ligand"""


    linesPDB = openPdbFile(PDB_ID)
    list_atom_ligand = []
    for line in linesPDB:
        if search ("^HETATM", line):
            atom = parsing.lineCoords(line)
            if atom["resName"] == ligand_ID and atom["element"] != "H":
                list_atom_ligand.append(atom)
    checkOnlyOneLigand(list_atom_ligand)
    connectMatrix = connectMatrixInPDB(PDB_ID)

    if connectAtom(connectMatrix, list_atom_ligand) == "False":
        print "Construt Matrix"
        calcul.buildConnectMatrix(list_atom_ligand, PDB_ID)

    return list_atom_ligand


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
    checkOnlyOneLigand(ligand)  ####retrieve only first ligand
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




def LigandInPDB(p_file_lig):
    """Load file result, PDB associated at ligand
    in: Name of file
    out: list ligand with PDB associated"""


    fileOpen = open(p_file_lig, "r")

    lineFile = fileOpen.readlines()
    d_out = {}

    for line in lineFile:
        line = line.split("\t")
        PDB = line[0]

        l_ligand = line[1]
        l_ligand = l_ligand.split("\n")[0]
        l_ligand = l_ligand.split(" ")
        if l_ligand == []:
            continue

        for ligand in l_ligand:
            if ligand != "":
                try:
                    d_out[ligand].append(PDB)
                except:
                    d_out[ligand] = []
                    d_out[ligand].append(PDB)

    return d_out




def globalPDB(PDB, ligand = ""):
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
    connectAtom(connectMatrix, out)
        # calcul.buildConnectMatrix(ligand, PDB) # why ????

    return out



def loadCloseStruct (pr_result, option_onePDB_ligand = 0) :
    
    if not path.isdir(pr_result) : 
        return None, None
    
    struct_neighbor = structure.neighborStruct()
    flag_file_empty = 0
    
    l_files = listdir(pr_result)

    
    if option_onePDB_ligand == 0 : 
        for name_file in l_files : 
            if search(".sum", name_file) : 
                if path.getsize(pr_result + name_file) != 0 : 
                    flag_file_empty = flag_file_empty + 1
                sub_struct = name_file.split ("_")[-1].split (".")[0]
                if sub_struct == "one" : 
                    continue
                else : 
                    if path.getsize(pr_result + name_file) != 0 : 
                        flag_file_empty = flag_file_empty + 1
                    if sub_struct == "global" : 
                        struct_global_temp = loadSummary(pr_result + name_file)
                    else : 
                        struct_neighbor[sub_struct] = loadSummary(pr_result + name_file)
        if flag_file_empty < 2 : 
            return None, None
        else : 
            struct_global_neighbor = {}
            struct_global_neighbor["global"] = struct_global_temp
            return struct_neighbor, struct_global_neighbor 
    else : 
        for name_file in l_files : 
            if search("one.sum", name_file) : 
                if path.getsize(pr_result + name_file) != 0 : 
                    flag_file_empty = flag_file_empty + 1
                sub_struct = name_file.split ("_")[-2].split (".")[0]
                if sub_struct == "global" : 
                    struct_global_temp = loadSummary(pr_result + name_file)
                else : 
                    struct_neighbor[sub_struct] = loadSummary(pr_result + name_file)
                    
        if flag_file_empty < 2 : 
            struct_neighbor_all, struct_global_neighbor_all = loadCloseStruct (pr_result, option_onePDB_ligand = 0)
            
            if struct_neighbor_all == None : 
                # run the extraction
                return None, None
            else : 
                # select one PDB by ligand
                for subs in  struct_neighbor.keys () : 
                    loadOnePDBbyLigand (struct_neighbor_all[subs], pr_result + "neighbor_" + str(subs) + "_one.sum")
                loadOnePDBbyLigand (struct_global_neighbor_all["global"], pr_result + "neighbor_global_one.sum")
                
                return struct_neighbor, struct_global_neighbor_all 
        else :
            struct_global_neighbor = {}
            struct_global_neighbor["global"] = struct_global_temp
            return struct_neighbor, struct_global_neighbor
            
        
    
def loadSummary (path_summary) : 
    
    
    l_out = []
    
    filin = open (path_summary, "r")
    l_lines = filin.readlines ()
    filin.close ()
    
    
    for l in l_lines : 
        d_line = {}
        l_s = l.split ("\t")
        d_line["PDB"] = l_s[0]
        d_line["serial"] = l_s[1].split ("/")[0]
        d_line["resName"] = l_s[1].split ("/")[1]
        d_line["x"] = float(l_s[1].split ("/")[2])
        d_line["y"] = float(l_s[1].split ("/")[3])
        d_line["z"] = float(l_s[1].split ("/")[4])
        d_line["neighbors"] = []
        for neigbor in l_s[-1].split ("//")[:-1] : 
            d_n = {}
            element_n = neigbor.split(" ")
            d_n["serial"] = element_n[0]
            d_n["resSeq"] = element_n[1]
            d_n["element"] = element_n[2]
            d_n["name"] = element_n[3]
            d_n["resName"] = element_n[4]
            d_n["distance"] = float(element_n[5])
            d_n["x"] = float(element_n[6])
            d_n["y"] = float(element_n[7])
            d_n["z"] = float(element_n[8])
            
            d_n["angle"] = []
            for angle in element_n [9:] :
                d_n["angle"].append (float(angle))
            d_line["neighbors"].append (d_n)
        l_out.append (d_line)
    return l_out     
    
    
    
def loadOnePDBbyLigand (st_all, p_filout, debug = 0):

    d_PDB = {}
    
    for atom_central in st_all : 
        if not atom_central["resName"] in d_PDB.keys () : 
            d_PDB[atom_central["resName"]] = []
        if not atom_central["PDB"] in  d_PDB[atom_central["resName"]] : 
            d_PDB[atom_central["resName"]].append (atom_central["PDB"])
    
    if debug == 1 : 
        for lig_id in d_PDB.keys () : 
            print lig_id, d_PDB[lig_id], "DEBUG"
    
    # retrieve best PDB
    for lig_ID in d_PDB.keys () : 
        if not len (d_PDB[lig_ID]) == 1 : 
            d_PDB[lig_ID] = checkPDBfile.SelectBestComplexRX (d_PDB[lig_ID])


    
    if debug == 1 :
        print "**********"
        print d_PDB 
        for lig_id in d_PDB.keys () : 
            print lig_id, d_PDB[lig_id], "DEBUG"
    


    # change struct    
    i = 0
    nb_lig = len (st_all)
    while i < nb_lig : 
        if not st_all[i]["PDB"] in d_PDB[st_all[i]["resName"]]: 
            del st_all[i]
            nb_lig = nb_lig - 1
            continue
        else : 
            i = i + 1
    
    # write file
    filout = open (p_filout, "w")
    for atom_central in st_all :
        lineWrite = str(atom_central["PDB"]) + "\t" + str(atom_central["serial"]) + "/" + str(atom_central["resName"]) + "/" + str(atom_central["x"]) + "/" + str(atom_central["y"]) +  "/" + str(atom_central["z"]) + "\t"
        for neighbor in atom_central["neighbors"]:
            lineWrite = lineWrite + str(neighbor["serial"]) + " " + str(neighbor["resSeq"]) + " " + str(neighbor["element"]) + " " + str(neighbor["name"]) + " " + str(neighbor["resName"]) + " " + str("%.2f" % neighbor["distance"]) + " " +str("%.3f" % neighbor["x"]) + " " +str("%.3f" % neighbor["y"]) + " " +str("%.3f" % neighbor["z"])  
            for angle in neighbor["angle"]:
                lineWrite = lineWrite + " " + str("%.2f" % angle)
            lineWrite = lineWrite + "//"
        lineWrite = lineWrite + "\n"
        filout.write (lineWrite)
    filout.close ()
    




