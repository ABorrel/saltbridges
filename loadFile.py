from re import search, sub
from os import path, listdir

import calcul
import parsing
import pathManage
import checkPDBfile



def openPdbFile(PDBin):
    """Open PDB file, check if there are in file many model structure and chose first model
        in: name file
        out: lines interest in list"""

    if path.exists(PDBin) : 
        filin = open(PDBin, "r")
    else : 
        rep = pathManage.openPdbFile()
        filin = open(rep + PDBin + ".pdb", "r")

    l_lines = filin.readlines()
    filin.close()

    int_nblines = len(l_lines)
    for i in range(0, int_nblines):
        # retrieve only first model
        if search("^ENDMDL", l_lines[i]):
            lineEnd = i
        if search("^CONECT", l_lines[i]):
            lineStartConect = i
            break

    if "lineStartConect" in locals() and "lineEnd" in locals():
        return l_lines[0: lineEnd] + l_lines[lineStartConect:int_nblines] # concatene lists 
    else:
        return l_lines



def resultFilterPDBLigand (path_file):
    """load result file with ligand and PDB ID associated
    in: name file result
    out: list of ligand"""

#     print path.getsize(path_file)
    # list empty
    if path.getsize(path_file) == 0 : 
        return []

    file_open = open(path_file, "r")
    lines = file_open.readlines()
    file_open.close()
    len_file = len(lines)

    l_out = []
    for i in range(0, len_file):
        d_lig = {}
        line = lines[i].split("\t")
        d_lig["name"] = line[0]
        l_pdb = line[1].split("\n")[0]
        l_pdb = l_pdb.split(" ")
        d_lig["PDB"] = []

        for pdb in l_pdb:
            if pdb != "":
                d_lig["PDB"].append(pdb)

        l_out.append(d_lig)

    return l_out




def ligandInPDB(PDBin, ligand_ID):
    """load l_atom_lig in structure in PDB file on single l_atom_lig by pdb and remove H
    out : atom in ligands
    out: list of atoms that l_atom_lig"""

    linesPDB = openPdbFile(PDBin)
    
    l_atom_lig = []
    for line in linesPDB:
        if search ("^HETATM", line):
            atom = parsing.lineCoords(line)
            if atom["resName"] == ligand_ID and atom["element"] != "H":
                l_atom_lig.append(atom)
    checkOnlyOneLigand(l_atom_lig)
    connectMatrix = connectMatrixInPDB(PDBin)

    if connectAtom(connectMatrix, l_atom_lig) == "False":
        print "Construt Matrix"
        calcul.buildConnectMatrix(l_atom_lig, PDBin)

    return l_atom_lig


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


def ligandInPDBConnectMatrixLigand(PDB_ID, ligand):
    """load ligand in structure in PDB file
    out : list atom in ligands with connect matrix calculated"""

    l_lines = openPdbFile(PDB_ID)
    l_atom = []
    l_serial = []
    for line in l_lines:
        if search ("^HETATM", line):
            atom = parsing.lineCoords(line)
            if atom["resName"] == ligand and atom["element"] != "H":
                if not atom["serial"] in l_serial : 
                    l_atom.append(atom)
                    l_serial.append(atom["serial"])
    checkOnlyOneLigand(l_atom)  ####retrieve only first ligand
    calcul.buildConnectMatrix(l_atom, PDB_ID)

    return l_atom



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

    for line in lineFile: ##### !!!!!! 
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




def loadCloseStruct (pr_result, control_empty_file = 1) :
    
    if not path.isdir(pr_result) : 
        return None
    
    d_summarize = {}
    flag_file_empty = 0
    
    l_files = listdir(pr_result)

    
    for name_file in l_files : 
        if search(".sum", name_file) : 
            if path.getsize(pr_result + name_file) == 0 : 
                flag_file_empty = flag_file_empty + 1
            else : 
                sub_struct = name_file.split ("_")[-1].split (".")[0]
                d_summarize[sub_struct] = loadSummary(pr_result + name_file)
                    
    if control_empty_file == 1 and flag_file_empty > 2 : # case file empty -> need control 
            # run the extraction
            
            print "== ERROR 2 files summarize empty"
            return None
    else : 
        if d_summarize == {} : 
            return None
    
    
    return d_summarize
            
        
    
def loadSummary (p_summary) : 
    
    
    l_out = []
    
    filin = open (p_summary, "r")
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
            
            d_n["angleSubs"] = []
            for angleSubs in element_n [9:] :
                if angleSubs == "NA" :                     
                    d_n["angleSubs"].append ("NA")
                else : 
                    d_n["angleSubs"].append (float(angleSubs))
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
            for angleSubs in neighbor["angleSubs"]:
                lineWrite = lineWrite + " " + str("%.2f" % angleSubs)
            lineWrite = lineWrite + "//"
        lineWrite = lineWrite + "\n"
        filout.write (lineWrite)
    filout.close ()
    



def ExtractInfoPDBID(PDB_ID) : 

    # control PDB exist in the folder where the PDB is included
    p_PDBfile = pathManage.pathDitrectoryPDB() + PDB_ID.lower() + ".pdb"
    if not path.exists(p_PDBfile) :
        print "ERROR load PDB ID -> ", PDB_ID
        return {} 
    
    # initialisation of the output
    d_out = {}
    d_out["protein"] = []
    d_out["RX"] = 100.0
    d_out["RFree"] = 100.0
    
    
    filin = open (p_PDBfile, "r")
    l_linesPDB = filin.readlines ()
    filin.close ()
    
    d_out["Header"] = l_linesPDB[0][6:].lower().strip ()
    
    for linePDB in l_linesPDB : 
        # Resolution
        if search("^REMARK   2 RESOLUTION", linePDB):
            
            lineRX = sub('[ ]{2,}', ' ', linePDB)
            try : d_out["RX"] = float (lineRX.split(" ")[3].replace (" ", ""))
            except : pass
            
        # Rfree
        elif search ("REMARK   3   R VALUE", linePDB) : 
            rfactor = linePDB.strip ().split (":")[-1].replace (" ", "")
            if rfactor == "NULL" : 
                rfactor = 0.0
            else : 
                rfactor =  float (rfactor)
            d_out["RFree"] =  rfactor
        
        # protein
        elif search ("^ATOM", linePDB) :
            atom_prot = parsing.lineCoords (linePDB, remove_H = 1)
            if atom_prot != None : 
                d_out["protein"].append (atom_prot)
        
        elif search ("^HETATM", linePDB) : 
            atom_HET = parsing.lineCoords (linePDB, remove_H = 1)
            if atom_HET != None : 
                name_lig = atom_HET["resName"]
                if not name_lig in d_out.keys () : 
                    d_out[name_lig] = []
                d_out[name_lig].append (atom_HET)
        
        # kept only first model in the protein in case of RMN structure
        elif search ("^ENDMDL", linePDB) : 
            break
        
        
    # separate the ligand in double
    for k in d_out.keys() : 
        if k != "protein" and k != "RX" and k != "RFree" and k != "Header" : 
            d_out[k] = parsing.separateByLigand (d_out[k], debug = 0)
            
    return d_out
                
