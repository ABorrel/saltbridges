from os import listdir, path 
from re import search
from copy import deepcopy
import calcul
import loadFile
import log
import parsing
import retrieveAtom
import structure
import tool
import toolSubstructure
import writeFile
import statistic
import loadFile
import managePDB



#############global search######################
def interestStructure (listAtomLigand, debug = 0):
    """For one serial_nitrogen atom search substructure
    in : list atom ligand
    out : list of substructures found (list string)"""

    l_serial_N = listAtomType (listAtomLigand, "N")
    l_serial_O = listAtomType (listAtomLigand, "O")
    
    if debug : 
        print l_serial_O
        print l_serial_N
    
    structureFound = []

    for serial_nitrogen in l_serial_N:
        listAtomConnectNitrogen, connect = retrieveAtom.atomConnect(listAtomLigand, serial_nitrogen)
        if imidazole(listAtomConnectNitrogen, listAtomLigand)[0] == 1:
            structureFound.append("Imidazole")
        elif guanidium(listAtomConnectNitrogen, listAtomLigand)[0] == 1 : 
            structureFound.append("Guanidium")    
#         elif diAmine(listAtomConnectNitrogen, listAtomLigand) == 1 : 
#             structureFound.append("Diamine")     
#         elif pyridine(listAtomConnectNitrogen, listAtomLigand) == 1 : 
#             structureFound.append("Pyridine")    
        elif cn (listAtomConnectNitrogen, listAtomLigand) == 1:
            structureFound.append("Primary")
        elif cnc(listAtomConnectNitrogen, listAtomLigand) == 1:
            structureFound.append("Secondary")
        elif cncc(listAtomConnectNitrogen, listAtomLigand) == 1:
            structureFound.append("Tertiary")
       

    for serial_oxygen in l_serial_O:
        listAtomConnectOxygen, connect = retrieveAtom.atomConnect(listAtomLigand, serial_oxygen)
        if acidCarboxylic(listAtomConnectOxygen, listAtomLigand)[0] >= 1:
            structureFound.append("AcidCarboxylic")
            
    return structureFound


def ligands(name_database, pr_init):
    '''search ligands in PDB database
    out : list of ligands with PDB files associated'''
    
    
    # control file exist
    if path.exists(pr_init + "resultLigandInPDB") and path.getsize(pr_init + "resultLigandInPDB") != 0: 
        return pr_init + "resultLigandInPDB"
    
    # import list PBD from file .dat
    listPDB = managePDB.retriveListPDB(name_database)
    
    start, fileLog = log.initAction("Search ligands in PDB")

    listPDBLigand = []

    for PDBFile in listPDB:
        namePDB = PDBFile.split(".")[0]
        fileLog.write(namePDB + "\n")
        ligandInPDB = structure.ligandPDB()
        try : linesPDB = loadFile.openPdbFile(namePDB)
        except : continue
        ligandInPDB["name"] = namePDB

        for linePDB in linesPDB:
            if(search ("^HETATM", linePDB)):
                atom = parsing.lineCoords(linePDB)
                if not atom["resName"] in ligandInPDB["ligands"]:
                    ligandInPDB["ligands"].append(atom["resName"])

        listPDBLigand.append(ligandInPDB)

    path_file = writeFile.resultLigandInPDB(listPDBLigand, pr_init)
    log.endAction("Search ligands in PDB", start, fileLog)
    return path_file




################substructure search -> stAtom###############

def cn(listAtomConnectNitrogen, listAtomLigan):
    """search primary stAtom 
    in: list atom connected of nitrogen, list atom ligand
    out: 1 or 0"""

    stAtom = toolSubstructure.matrixElement(listAtomConnectNitrogen)
    if stAtom == ["N", "C"] or stAtom == ["N", "C"]:
        if toolSubstructure.checkSingleBond(listAtomConnectNitrogen[0], listAtomConnectNitrogen[1]) == 1:
            if toolSubstructure.checkConectOnlyC(listAtomConnectNitrogen[1], listAtomLigan) == 1 : 
                return 1
    return 0


def cnc(listAtomConnectNitrogen, listAtomLigand):
    """search secondary stAtom 
    in: list atom connected of nitrogen, list atom ligand
    out: 1 or 0"""

    stAtom = toolSubstructure.matrixElement(listAtomConnectNitrogen)

    if stAtom == ["N", "C", "C"]:
        if toolSubstructure.checkConectOnlyC(listAtomConnectNitrogen[1], listAtomLigand) == 1 and toolSubstructure.checkConectOnlyC(listAtomConnectNitrogen[2], listAtomLigand) == 1:
            if toolSubstructure.checkSingleBond(listAtomConnectNitrogen[0], listAtomConnectNitrogen[1]) == 1 and toolSubstructure.checkSingleBond(listAtomConnectNitrogen[0], listAtomConnectNitrogen[2]) == 1:
                return 1
    return 0


def cncc(listAtomConnectNitrogen, listAtomligand):
    """search tertiary stAtom 
    in: list atom connected of nitrogen, list atom ligand
    out: 1 or 0"""

    stAtom = toolSubstructure.matrixElement(listAtomConnectNitrogen)

    if stAtom == ["N", "C", "C", "C"]:
        if toolSubstructure.checkConectOnlyC(listAtomConnectNitrogen[1], listAtomligand) == 1 and toolSubstructure.checkConectOnlyC(listAtomConnectNitrogen[2], listAtomligand) == 1 and toolSubstructure.checkConectOnlyC(listAtomConnectNitrogen[3], listAtomligand) == 1:
            if toolSubstructure.checkCoplanar(listAtomConnectNitrogen[0], listAtomligand) == 1:
                if toolSubstructure.checkSingleBond(listAtomConnectNitrogen[0], listAtomConnectNitrogen[1]) == 1 and toolSubstructure.checkSingleBond(listAtomConnectNitrogen[0], listAtomConnectNitrogen[2]) == 1 and toolSubstructure.checkSingleBond(listAtomConnectNitrogen[0], listAtomConnectNitrogen[3]) == 1:
                    return 1
    return 0


def guanidium(l_at_connect_N, listAtomLigand):
    """search guanidium 
    in: list atom connected of nitrogen, list atom ligand
    out: 1 or 0"""
    
    stAtom = toolSubstructure.matrixElement(l_at_connect_N)
    findNH = 0
    # atom lateral
    if stAtom == ["N", "C"] :
        l_at_connect_c, connectC = retrieveAtom.atomConnect(listAtomLigand , l_at_connect_N[1]["serial"])
        if connectC == ["C", "N", "N", "N"] :
            groupAtomN1, conect_N1 = retrieveAtom.atomConnect(listAtomLigand , int (l_at_connect_c[0]["connect"][1]))
            groupAtomN2, conect_N2 = retrieveAtom.atomConnect(listAtomLigand , int (l_at_connect_c[0]["connect"][2]))
            groupAtomN3, conect_N3 = retrieveAtom.atomConnect(listAtomLigand , int (l_at_connect_c[0]["connect"][3]))
            l_conect = [conect_N1, conect_N2, conect_N3]
            l_group_atom =  [groupAtomN1, groupAtomN2, groupAtomN3]
            
            i = 0
            while i < 3 :  
                if l_conect[i] == ["N", "C"] : 
                    findNH = findNH + 1
                    i = i + 1
                elif l_conect[i] == ["N", "C", "C"]:  
                    if l_group_atom[i][1]["serial"] != l_at_connect_c[0]["serial"] : 
                        serial_out = l_group_atom[i][1]["serial"]
                    elif l_group_atom[i][2]["serial"] != l_at_connect_c[0]["serial"] : 
                        serial_out = l_group_atom[i][2]["serial"]
                    else : 
                        print "ERROR l158"
                        return [0, []]
                    i = i + 1
                else :
                    return [0, []]
                
            # check number primary stAtom -> case GAI not take, change ?
            if findNH == 2  : 
                return [1, [l_at_connect_c[0]["serial"], groupAtomN1[0]["serial"], groupAtomN2[0]["serial"], groupAtomN3[0]["serial"], serial_out]]
                    
            
        else :
            return [0, []]
    
    # atom central structure
    elif stAtom == ["N", "C", "C"] : 
        for at_conect in l_at_connect_N[1:] : 
            l_group_at, conect_N = retrieveAtom.atomConnect(listAtomLigand , at_conect["serial"])
            if conect_N == ["C", "N", "N", "N"] :
                l_c_central = l_group_at
                for group_at in l_group_at[1:] : 
                    l_goup_N, connect_N = retrieveAtom.atomConnect(listAtomLigand , group_at["serial"])
                    if connect_N == ["N", "C"] : 
                        findNH = findNH + 1
            else : 
                serial_c = l_group_at[0]["serial"]

        if findNH >= 2 and "serial_c" in locals() : 
            return [1,[ serial_c, l_c_central[0]["serial"], l_c_central[1]["serial"], l_c_central[2]["serial"], l_c_central[3]["serial"]]]
        else : 
            return [0, []]
        
    else : 
        return [0, []]
    
    return [0, []]
        
        
#         
#         
#         for atomConnectNitrogen in l_at_connect_N : 
#             if atomConnectNitrogen ["element"] == "C" : 
#                 listAtomC, connectElementC = retrieveAtom.atomConnect(listAtomLigand, atomConnectNitrogen["serial"])
#                 if connectElementC == ["C", "N", "N", "N"] : 
#                     for atomConnectC in listAtomC[1:] : 
#                         nitrogenConnect, connectNitrogen = retrieveAtom.atomConnect(listAtomLigand, atomConnectC["serial"])
#                         if connectNitrogen == ["N", "C"] : 
#                             findNH = findNH + 1
#         
#         if findNH == 2 : 
#             return 1
#     
#     return [0, []]


# def guanidiumNitrogenCentral(l_at_connect_N, listAtomLigand):
#     """Search guanidium 
#     in: list Atom connected of nitrogen, list global atoms
#     out: boolean"""
#     
#     findNH = 0
#     for atomConnectNitrogen in l_at_connect_N : 
#         if atomConnectNitrogen ["element"] == "C" : 
#             listAtomC, connectElementC = retrieveAtom.atomConnect(listAtomLigand, atomConnectNitrogen["serial"])
#             if connectElementC == ["C", "N", "N", "N"] : 
#                 for atomConnectC in listAtomC[1:] : 
#                     nitrogenConnect, connectNitrogen = retrieveAtom.atomConnect(listAtomLigand, atomConnectC["serial"])
#                     if connectNitrogen == ["N", "C"] : 
#                         findNH = findNH + 1
#     
#     if findNH == 2 : 
#         return 1
#     
#     return 0
                
    

def pyridine(listAtomConnectNitrogen, listAtomLigand):
    """search pyridine 
    in: list atom connected of nitrogen, list atom ligand
    out: boolean"""
        
    stAtom = toolSubstructure.matrixElement(listAtomConnectNitrogen)
    
    if stAtom == ["N", "C", "C"] : 
        nitrogenInit = listAtomConnectNitrogen[0]
        if  cycleOnlyTestCarbon(nitrogenInit["serial"], nitrogenInit["serial"], nitrogenInit["serial"], listAtomLigand, 6, 0) == 1 : 
            return 1
        
         
def diAmine (l_at_connect_N, listAtomLigand):
    """search diamine 
    in: list atom connected of nitrogen, list atom ligand
    out: boolean"""
        
    connect_element = toolSubstructure.matrixElement(l_at_connect_N)
#     print connect_element, "l199"
    if connect_element != ["N", "C"] :
        return 0
    else : 
        l_connect_C1, connect_C1 = retrieveAtom.atomConnect(listAtomLigand, l_at_connect_N[0]["connect"][1]) 
#         print connect_C1, "l204"
        if connect_C1 == ["C", "C", "N", "N"] or connect_C1 == ["C", "N", "C", "N"] or connect_C1 == ["C", "N", "N", "C"] : 
            for atom_connect in l_connect_C1[:1] :
                l_connect_atom, connect_atom = retrieveAtom.atomConnect(listAtomLigand, atom_connect ["serial"])
                
                if atom_connect ["element"] == "N" : 
                    if connect_atom != ["N", "C"]  : 
#                         print "l211"
                        return 0
                else  : 
                    if connect_atom == ["C", "C", "C", "C"] :
#                         print "l214" 
                        return 0
        else : 
#             print connect_C1
            return 0
    return 1
            


def acidCarboxylic(listAtomConnectOx, listAtomLigand) : 
    
    l_out = []
    connect_element = toolSubstructure.matrixElement(listAtomConnectOx)
    if connect_element != ["O", "C"] : 
        return [0, []]
    else : 
        l_atom_connect_C, connect_matrix_C = retrieveAtom.atomConnect(listAtomLigand, listAtomConnectOx[1]["serial"])
        l_out.append (l_atom_connect_C[0]["serial"])
        if connect_matrix_C == ["C", "O", "O", "C"] or connect_matrix_C == ["C", "O", "C", "O"] or  connect_matrix_C == ["C", "C", "O", "O"] : 
            for atom_connect_C in l_atom_connect_C[1:] :
                l_out.append (atom_connect_C ["serial"])
                if atom_connect_C["element"] == "O" : 
                    l_atom_connect_ox, connect_matrix_ox = retrieveAtom.atomConnect(listAtomLigand, atom_connect_C["serial"])
                    if connect_matrix_ox != ["O", "C"] :
                        return [0, []]
                    else : 
                        if not atom_connect_C["serial"] in l_out : 
                            l_out.append (atom_connect_C["serial"])
                            
        else : 
            return [0, []]
    
    return [1,l_out]


    
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

def imidazole(listAtomConnectNitrogen,list_atom_ligand):
    """search imidazole global
    in: list atom connected of nitrogen, list atom ligand
    out: boolean"""

    stAtom = toolSubstructure.matrixElement(listAtomConnectNitrogen)
    
    if stAtom == ["N", "C", "C"]:
        l_atom_check = [listAtomConnectNitrogen[0]["serial"]]
        groupAtomC1, conect_C1 = retrieveAtom.atomConnect(list_atom_ligand, int (listAtomConnectNitrogen[0]["connect"][1]))
        groupAtomC2, conect_C2 = retrieveAtom.atomConnect(list_atom_ligand, int (listAtomConnectNitrogen[0]["connect"][2]))
#         print conect_C1, "C1"
#         print conect_C2, "C2"
        if conect_C1 == ["C", "C", "N", "N"] or conect_C1 == ["C", "N", "N", "C"] or conect_C1 == ["C", "N", "C", "N"] or conect_C1 == ["C", "N", "N"]:
            l_atom_check.append (groupAtomC1[0]["serial"])
            if imidazoleATOM3(groupAtomC1, l_atom_check, list_atom_ligand)[0] == 1:
#                 print "l344"
                return [1, l_atom_check]

        elif conect_C2 == ["C", "C", "N", "N"] or conect_C2 == ["C", "N", "N", "C"] or conect_C2 == ["C", "N", "C", "N"] or conect_C2 == ["C", "N", "N"]:
#             print "IN C2"
#             print l_atom_check, "l348"
            l_atom_check.append (groupAtomC2[0]["serial"])
            if imidazoleATOM3(groupAtomC2, l_atom_check, list_atom_ligand)[0] == 1:
#                 print "l350"
                return [1, l_atom_check]

    return [0, []]


def imidazoleATOM3(listAtom, l_atom_check, listAtomLigand):
    """Check the atom 3 in circle of imidazole
    in : atoms ligand in list, serial of first nitrogen, name of ligand
    out : boolean"""

    for atom in listAtom:
        if atom["element"] == "N":
            if not atom["serial"] in l_atom_check : 
                groupAtomN1, conect_N1 = retrieveAtom.atomConnect(listAtomLigand, atom["serial"])



#     print conect_N1, "CN1, l367"
    if conect_N1 == ['N', 'C', 'C']:
        l_atom_check.append (groupAtomN1[0]["serial"])
#         print l_atom_check, "l370"
        if imidazoleATOM4(groupAtomN1, l_atom_check, listAtomLigand)[0] == 1:
            return [1, l_atom_check]

    return [0, l_atom_check]


def imidazoleATOM4(listAtomNitrogen, l_atom_check, listAtomLigand):
    """Check the atom 4 in cercle of imidazole
    in : atoms ligand in list, serial of first nitrogen, name of ligand
    out : boolean"""
    
#     print l_atom_check, "l387"
    for atom in listAtomNitrogen:
        if atom["element"] == "C":
            if not atom["serial"] in l_atom_check : 
                groupAtomC3, conect_C3 = retrieveAtom.atomConnect(listAtomLigand, atom["serial"])
                l_atom_check.append (groupAtomC3[0]["serial"])
#                 print l_atom_check, "l392"
                if  imidazoleATOM5(groupAtomC3, l_atom_check, listAtomLigand)[0] == 1:
                    return [1, l_atom_check]

    return [0, l_atom_check]
            

def imidazoleATOM5(l_at_c4, l_atom_check, listAtomLigand):
    """Check the atom 5 in cercle of imidazole
    in : atoms ligand in list, serial of first nitrogen, name of ligand
    out : boolean"""
    
    
    for at_c4 in l_at_c4 : 
        if at_c4["element"] == "C" and not at_c4["serial"] in l_atom_check:
            l_conect_c4, conect_c4 =  retrieveAtom.atomConnect(listAtomLigand, at_c4["serial"])
            
            for at_connect in l_conect_c4 : 
                if at_connect["element"] == "N" and at_connect["serial"] in l_atom_check :
                    l_atom_check.append (at_c4["serial"])
                    return [1, l_atom_check]
    return [0, l_atom_check]
    

######################################################################################################################


def globalSearch (max_distance, p_file_dataset,  pr_result, option_one_PDB = 0):
    
    
    pr_summary = pr_result + "Sum" +  p_file_dataset.split ("/")[-1][8:] + "/"
    
    # load structure in summary ---> if use need place option one PDB by ligand
    struct_neighbor, struct_global_neighbor = loadFile.loadCloseStruct (pr_summary, option_one_PDB)
    if struct_neighbor != None : 
        print "********INNNNNN*************"
        print type (struct_global_neighbor["global"])
        return struct_neighbor, struct_global_neighbor
     
    
    #start, logFile = log.initAction("search neighbors in " +path.basename(p_file_dataset))
     
    list_ligands_in_PDB = loadFile.resultFilterPDBLigand(p_file_dataset)
    nbLigand = len(list_ligands_in_PDB)
    
    # ##Count Structure
    struct_neighbor = structure.neighborStruct()
    struct_global_neighbor = []

    
    # ##Write summary file
    files_summary = writeFile.openFileSummary(pr_summary)# sumary result
    
    # inialization    
    i = 0
    while i < nbLigand :
        print "Ligand: " + str(list_ligands_in_PDB[i]["name"]) + " " + str(i) + " " + str(nbLigand)
        nbPDB = len(list_ligands_in_PDB[i]["PDB"])
        
        # take only one PDB by ligand
            
        j = 0
        while j < nbPDB : 
            name_PDB = list_ligands_in_PDB[i]["PDB"][j]
            list_atom_ligand = loadFile.ligandInPDB(list_ligands_in_PDB[i]["PDB"][j], list_ligands_in_PDB[i]["name"])
            #print list_atom_ligand
            
            # search neighbor for every atom in ligand selected
            globalNeighbors(max_distance, list_atom_ligand, name_PDB, struct_global_neighbor)
            
            # search neighbor for interest 
            interestGroup(max_distance, list_atom_ligand, name_PDB, struct_neighbor)
            
            j = j + 1
        i = i + 1
            
    
    writeFile.neighborStruct(struct_neighbor, struct_global_neighbor, files_summary)
    writeFile.closeFileSummary(files_summary)
    
    if option_one_PDB == 1 : 
        struct_neighbor, struct_global_neighbor = loadFile.loadCloseStruct (pr_summary,  option_onePDB_ligand = option_one_PDB)
    
    # case where load directly substructure
    if not "global" in struct_global_neighbor.keys () : 
        struct_global = {}
        struct_global["global"] = struct_global_neighbor
        return struct_neighbor, struct_global
    else :
        struct_global = struct_global_neighbor
    
    return struct_neighbor, struct_global
            


def interestGroup (max_distance, list_atom_ligand, name_PDB, struct_neighbor):
    """Search different groups
    in : ligands in namePDB
    out : major nitrogen in the group of different structures"""
    
    listSerialNitrogen = listAtomType(list_atom_ligand, "N")
    listSerialOxygen = listAtomType(list_atom_ligand, "O")
    
    d_imd_temp = {}
    d_imd_temp["Imidazole"] = []
    d_gua_temp = {}
    d_gua_temp["Guanidium"] = []   
    d_dia_temp = {}
    d_dia_temp["Diamine"] = []  
    d_acd_temp = {}
    d_acd_temp["AcidCarboxylic"] = []  
    
    # different struct_neighbor
    for serialNitrogen in listSerialNitrogen:
        listAtomConnectNitrogen, conect = retrieveAtom.atomConnect(list_atom_ligand, serialNitrogen)
        # check every substructure
        if imidazole(listAtomConnectNitrogen, list_atom_ligand)[0] == 1:
            implementNeighborStruct (max_distance, listAtomConnectNitrogen, name_PDB, list_atom_ligand, "Imidazole", d_imd_temp)
        elif guanidium(listAtomConnectNitrogen, list_atom_ligand)[0] == 1:
            implementNeighborStruct (max_distance, listAtomConnectNitrogen, name_PDB, list_atom_ligand, "Guanidium", d_gua_temp)
#         elif diAmine(listAtomConnectNitrogen, list_atom_ligand) == 1:
#             implementNeighborStruct (max_distance, listAtomConnectNitrogen, name_PDB, list_atom_ligand, "Diamine", d_dia_temp)
#         elif pyridine(listAtomConnectNitrogen, list_atom_ligand) == 1:
#             implementNeighborStruct (max_distance, listAtomConnectNitrogen, name_PDB, list_atom_ligand, "Pyridine", struct_neighbor)
        elif cn(listAtomConnectNitrogen, list_atom_ligand) == 1:
            implementNeighborStruct (max_distance, listAtomConnectNitrogen, name_PDB, list_atom_ligand, "Primary", struct_neighbor)
        elif cnc(listAtomConnectNitrogen, list_atom_ligand) == 1:
            implementNeighborStruct (max_distance, listAtomConnectNitrogen, name_PDB, list_atom_ligand, "Secondary", struct_neighbor)
        elif cncc(listAtomConnectNitrogen, list_atom_ligand) == 1:
            implementNeighborStruct (max_distance, listAtomConnectNitrogen, name_PDB, list_atom_ligand, "Tertiary", struct_neighbor)
            
               
    for serial_oxygen in listSerialOxygen :
        listAtomConnectOx, conect = retrieveAtom.atomConnect(list_atom_ligand, serial_oxygen)
        if acidCarboxylic(listAtomConnectOx, list_atom_ligand)[0]== 1:
            implementNeighborStruct (max_distance, listAtomConnectOx, name_PDB, list_atom_ligand, "AcidCarboxylic", d_acd_temp)
            
    # group neighbor
    if len(d_acd_temp["AcidCarboxylic"]) > 0 : 
        regroupAtomNeighbor(d_acd_temp["AcidCarboxylic"], struct_neighbor["AcidCarboxylic"], list_atom_ligand)
    if len(d_gua_temp["Guanidium"]) > 0 : 
        regroupAtomNeighborGuanidium(d_gua_temp["Guanidium"], struct_neighbor["Guanidium"], list_atom_ligand)
    if len(d_imd_temp["Imidazole"]) > 0 : 
        regroupAtomNeighbor(d_imd_temp["Imidazole"], struct_neighbor["Imidazole"], list_atom_ligand)
    if len(d_dia_temp["Diamine"]) > 0 : 
        regroupAtomNeighbor(d_dia_temp["Diamine"], struct_neighbor["Diamine"], list_atom_ligand)

 

#######regroup neighbors case of imidazole, guanidium and diamine###########


def implementNeighborStruct (max_distance, list_atom_connect_central, name_PDB, l_atom_ligand, chemical_found, struct_neighbor):
    
    
    atom_neighbors = buildAtom(max_distance, list_atom_connect_central[0], name_PDB, chemical_found, l_atom_ligand)
    
    # !!!!!!!!!!!!!!!!!!!!!!!!
    # !!!!!!!!!
#     if option_angle == 1 : -> check with count
#         checkAngleInSearchNeighbor(atom_neighbors, chemical_found)
    #if len(atom_neighbors["neighbors"]) != 0 : 
    
    
    struct_neighbor[chemical_found].append(atom_neighbors)


def regroupAtomNeighborGuanidium(listAtomRegroup, stAtom, listAtomLigand):  # #A revoir tres tres lourds en temps -> BUG
    """Regroup neighbors for 3 nitrogen atoms guanidium
    in: list atom finds in search neighbors, count structure stAtom, list atom of ligand in pdb
    out: modification stAtom structure"""
    
    
    for atomRegroup in listAtomRegroup : 
        atomRegroupConnect, connectAtomRegroup = retrieveAtom.atomConnect(listAtomLigand, atomRegroup["serial"])
        if connectAtomRegroup == ["N", "C"] : 
            atom1 = atomRegroup 
            
    outAtom = {}
    outAtom["PDB"] = atom1["PDB"]
    outAtom["resName"] = atom1["resName"]
    outAtom["serial"] = atom1["serial"]
    outAtom["x"] = atom1["x"]
    outAtom["y"] = atom1["y"]
    outAtom["z"] = atom1["z"]
    outAtom["neighbors"] = deepcopy(atom1["neighbors"]) 
    
    listAtomConnectN, connectMatrixN = retrieveAtom.atomConnect(listAtomLigand, atom1["serial"])
    listAtomConnectC, connectMatrixC = retrieveAtom.atomConnect(listAtomLigand, listAtomConnectN[1]["serial"])
    
    for atomConnectC in listAtomConnectC : 
        for atom in listAtomRegroup :
            if atom["serial"] == atomConnectC["serial"] : 
                for neighbor in atom["neighbors"] : 
                    if not neighbor in outAtom["neighbors"] : 
                        outAtom["neighbors"].append(neighbor)
    
    stAtom.append(outAtom)
    



def regroupAtomNeighbor (atomRegroup, struct_neighbor, listAtomLigand):
    """regroup equivalent atom
    in: atom no regroup, struct_neighbor structure neighbor, list atom ligand
    out: struct_neighbor modified"""
    
    
    if len(atomRegroup) == 2 : 
        struct_neighbor.append(regroupNeighbor(atomRegroup[0]["serial"], atomRegroup[1]["serial"], atomRegroup))
        
        
    else : 
        l_atom_serial = []
        for atom in atomRegroup : 
            l_atom_serial.append(atom["serial"])
        nbSerial = len(l_atom_serial)
      
        i = 0
        while i < nbSerial : 
            l_atom_bond, conect = retrieveAtom.atomConnect(listAtomLigand, l_atom_serial[i])
            
            for atom_bond in l_atom_bond[1:] : 
                l_intersect = list(set(l_atom_serial) & set(atom_bond["connect"]))# intersect list
                if len(l_intersect) > 2 : 
                    print "ERROR regroup (searchPDB line 625)"
                elif len(l_intersect) == 2 : 
                    struct_neighbor.append(regroupNeighbor(atomRegroup[0]["serial"], atomRegroup[1]["serial"], atomRegroup))
                    l_atom_serial.remove (l_intersect[0])
                    l_atom_serial.remove (l_intersect[1])
                    nbSerial = nbSerial -2
                    continue
                else : 
                    pass
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
    outAtom["x"] = atom1["x"]
    outAtom["y"] = atom1["y"]
    outAtom["z"] = atom1["z"]
    
    outAtom["neighbors"] = deepcopy(atom1["neighbors"])
    
    for neighbor in atom2["neighbors"] : 
        outAtom["neighbors"].append(neighbor)
    
#    print len(outAtom["neighbors"])
    return outAtom
    

######################################################################################


def buildAtom(rayon, at_central, namePDB, typeStudy, listAtomLigand):
    """Build count structure for stAtom structure
    in: distance, serial atom nitrogen, type study, list atom ligand
    out: count strucutre"""
    
    atom = {}
    atom["PDB"] = namePDB
    atom["resName"] = at_central["resName"]
    atom["serial"] = at_central["serial"]
    atom["x"] = at_central["x"]
    atom["y"] = at_central["y"]
    atom["z"] = at_central["z"]
    
    atom["neighbors"] = neighbors(rayon, at_central, namePDB, typeStudy, listAtomLigand)

    return atom



def neighbors(rayon, central_atom, pdb, typeStructure = "global", ligandPDB = [] ):
    """Search neighbors for all ligand
    in : rayon where is atoms, central atom, pdb file
    out : list atoms found"""

    listAtom = []
    linesPDB = loadFile.openPdbFile(pdb)
    for line in linesPDB:
        if search("^ATOM", line) or search("^HETATM", line): 
            atom = parsing.lineCoords(line)
            if atom != {} and atom["element"] != "H":
                distance = calcul.distanceTwoatoms(central_atom, atom)
                if distance <= rayon and distance != 0.0:
                    if central_atom["resSeq"] != atom["resSeq"]: # check if variation
                        if tool.atomInList(listAtom, atom) == 0:
                            atom["distance"] = distance
                            atom["angle"] = calcul.angle(central_atom, atom, ligandPDB, typeStructure)
                            atom["classification"] = structure.classificationATOM(atom)
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


def globalNeighbors(distance, l_atom_ligand, pdb, struct_global_neighbor):
    """Retrieve neighbor in sphere
    in : distance arroud  atom, atom in ligand, name of pdb
    out : list of atoms"""
    

    for atom in l_atom_ligand:
        groupAtom, conect = retrieveAtom.atomConnect(l_atom_ligand, atom["serial"])
        atom = buildAtom(distance, groupAtom[0], pdb, "global", l_atom_ligand)
        struct_global_neighbor.append(atom)



def listAtomType(groupAtom, type_atom):
    """Search N in the atoms list
    return a list of number atom"""

    list_out = []
    for i in groupAtom:
        if i["element"] == type_atom:
            serial = i["serial"]
            if not serial in list_out:
                list_out.append(i["serial"])
    return  list_out






def checkAngleInSearchNeighbor(atomRetrieve, typeStructStudy):
    """Check for every neighbors of atom retrieve the angle
    in: atomRetrieve
    out: atomRetrieve modified"""
    
    conditionAngle = structure.selectionAngle()
    nbNeighbors = len(atomRetrieve["neighbors"])
        
    i = 0
    while i < nbNeighbors :
        # print atomRetrieve["neighbors"]
        
        if checkListAngle(atomRetrieve["neighbors"][i]["angle"], conditionAngle[typeStructStudy]) == 0 : 
            del atomRetrieve["neighbors"][i]
            nbNeighbors = nbNeighbors - 1
            continue
        else : 
            i = i + 1
    

def checkListAngle (list_angle, condition): 
    
    
    for angle in list_angle : 
        if angle < condition["INF"] : 
            return 0
        if angle > condition["SUP"] : 
            return 0
    
    return 1
        
        
    
    
    
      



##################################Serine Protease #############################################

# def cycle(atomInitial, firstPosition, atomLigand, nbTest, listRetrieve, serialTest):
#
#    #print nbTest
#    if atomInitial == 0 or len(atomInitial["connect"]) < 3 or atomInitial["serial"] == serialTest : 
# #        print "OUT1"
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
# #            print "OUT2"
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
# def cycle2(atomInitial, firstPosition, atomLigand, nbTest, strSerialList, serialTest):
#
#    #print nbTest
#    if atomInitial == 0 or len(atomInitial["connect"]) < 3 or atomInitial["serial"] == serialTest : 
# #        print "OUT1"
#        return 0, ""
#    else : 
#        strSerialList = strSerialList + "_" + str(atomInitial["serial"])
#    
#    if nbTest == 0 : 
#        if firstPosition in atomInitial["connect"] [1:] and firstPosition != serialTest:
#            return 1 , strSerialList
#        else : 
# #            print "OUT2"
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


# def checkCycle (atom, firstPosition, lastPosition, n, atomLigand):
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



