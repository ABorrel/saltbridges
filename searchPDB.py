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
import pathManage



#############global search######################
def interestStructure (l_atom_lig, more_flex = 0, debug = 0):
    """For one serial_nitrogen atom search substructure
    in : list atom ligand
    out : list of substructures found (list string)"""

    l_serial_N = ListSerialElement (l_atom_lig, "N")
    l_serial_O = ListSerialElement (l_atom_lig, "O")
    
    if debug : 
        print l_serial_O
        print l_serial_N
    
    l_substruct = []

    for serial_nitrogen in l_serial_N:
        l_atom_connectN, connect = retrieveAtom.atomConnect(l_atom_lig, serial_nitrogen)
        if imidazole(l_atom_connectN, l_atom_lig)[0] == 1:
            l_substruct.append("IMD")
        elif Guanidium(l_atom_connectN, l_atom_lig)[0] == 1 : 
            l_substruct.append("GAI")    
#         elif diAmine(l_atom_connectN, l_atom_lig) == 1 : 
#             l_substruct.append("Diamine")     
#         elif pyridine(l_atom_connectN, l_atom_lig) == 1 : 
#             l_substruct.append("Pyridine")    
        elif cn (l_atom_connectN, l_atom_lig) == 1:
            l_substruct.append("I")
        elif cnc(l_atom_connectN, l_atom_lig, more_flex = more_flex) == 1:
            l_substruct.append("II")
        elif cncc(l_atom_connectN, l_atom_lig, more_flex = more_flex) == 1:
            l_substruct.append("III")
       

    for serial_oxygen in l_serial_O:
        l_atom_connectO, connect = retrieveAtom.atomConnect(l_atom_lig, serial_oxygen)
        if acidCarboxylic(l_atom_connectO, l_atom_lig)[0] == 1:
            l_substruct.append("COO")
    
    return l_substruct


def ligands(name_database, pr_init):
    '''search ligands in PDB database
    out : list of ligands with PDB files associated'''
    
    print "Start Search Ligand In PDB file"
    # control file exist
    if path.exists(pr_init + "resultLigandInPDB") and path.getsize(pr_init + "resultLigandInPDB") != 0: 
        return pr_init + "resultLigandInPDB"
    
    # import list PBD from file .dat
    # http://www.rcsb.org/pdb/rest/representatives?cluster=50
    l_PDB = managePDB.retriveListPDB(name_database)
    
    l_d_lig = []

    for PDB_ID in l_PDB:
        d_lig_PDB = structure.ligandPDB()
        d_lig_PDB["name"] = PDB_ID
        l_lig_in = parsing.retrieveListLigand(pathManage.pathDitrectoryPDB() + PDB_ID.lower() + ".pdb")
        d_lig_PDB["ligands"] = l_lig_in
        l_d_lig.append(d_lig_PDB)
    
    # write result file
    p_out = writeFile.resultLigandInPDB(l_d_lig, pr_init)

    print "END Search Ligand In PDB file"
    return p_out




################substructure search -> stAtom###############

def cn(listAtomConnectNitrogen, l_atom_lig):
    """search primary stAtom 
    in: list atom connected of nitrogen, list atom ligand
    out: 1 or 0"""

    stAtom = toolSubstructure.matrixElement(listAtomConnectNitrogen)
    if stAtom == ["N", "C"] or stAtom == ["N", "C"]:
        if toolSubstructure.checkSingleBond(listAtomConnectNitrogen[0], listAtomConnectNitrogen[1]) == 1:
            if toolSubstructure.checkConectOnlyC(listAtomConnectNitrogen[1], l_atom_lig) == 1 : 
                return 1
    return 0


def cnc(l_atom_connectN, l_atom_lig, more_flex = 0):
    """search secondary connect_element 
    in: list atom connected of nitrogen, list atom ligand
    out: 1 or 0"""

    connect_element = toolSubstructure.matrixElement(l_atom_connectN)

    if connect_element == ["N", "C", "C"]:
#         print "IN - NH2"
        if more_flex == 1 : # more flexible control just connectivity
            #if toolSubstructure.checkConectOnlyC(l_atom_connectN[1], l_atom_lig) == 1 and toolSubstructure.checkConectOnlyC(l_atom_connectN[2], l_atom_lig) == 1:
            #if toolSubstructure.checkSingleBond(l_atom_connectN[0], l_atom_connectN[1]) == 1 and toolSubstructure.checkSingleBond(l_atom_connectN[0], l_atom_connectN[2]) == 1:
            return 1
        else : 
            if toolSubstructure.checkConectOnlyC(l_atom_connectN[1], l_atom_lig) == 1 and toolSubstructure.checkConectOnlyC(l_atom_connectN[2], l_atom_lig) == 1:
                if toolSubstructure.checkSingleBond(l_atom_connectN[0], l_atom_connectN[1]) == 1 and toolSubstructure.checkSingleBond(l_atom_connectN[0], l_atom_connectN[2]) == 1:
                    return 1
    return 0


def cncc(l_atom_connectN, l_atom_lig, more_flex = 0):
    """
    Search tertiary amine in list of atom lig
    - Append option 
    """

    connect_element = toolSubstructure.matrixElement(l_atom_connectN)
    
    if connect_element == ["N", "C", "C", "C"]:
        
        if more_flex == 1 : 
            #if toolSubstructure.checkCoplanar(l_atom_connectN[0], l_atom_lig) == 1:
            #    if toolSubstructure.checkSingleBond(l_atom_connectN[0], l_atom_connectN[1], d_min = 1.34) == 1 and toolSubstructure.checkSingleBond(l_atom_connectN[0], l_atom_connectN[2], d_min = 1.34) == 1 and toolSubstructure.checkSingleBond(l_atom_connectN[0], l_atom_connectN[3], d_min = 1.34) == 1:
            return 1
        else :     
            if toolSubstructure.checkConectOnlyC(l_atom_connectN[1], l_atom_lig) == 1 and toolSubstructure.checkConectOnlyC(l_atom_connectN[2], l_atom_lig) == 1 and toolSubstructure.checkConectOnlyC(l_atom_connectN[3], l_atom_lig) == 1:
                if toolSubstructure.checkCoplanar(l_atom_connectN[0], l_atom_lig) == 1:
                    if toolSubstructure.checkSingleBond(l_atom_connectN[0], l_atom_connectN[1]) == 1 and toolSubstructure.checkSingleBond(l_atom_connectN[0], l_atom_connectN[2]) == 1 and toolSubstructure.checkSingleBond(l_atom_connectN[0], l_atom_connectN[3]) == 1:
                        return 1
    return 0


def Guanidium(l_at_connect_N, l_atom_lig):
    """search Guanidium 
    in: list atom connected of nitrogen, list atom ligand
    out: 1 or 0"""
    
    connectN1 = toolSubstructure.matrixElement(l_at_connect_N)
    nb_NH = 0
    l_atom_out = []
    # atom lateral
    if connectN1 == ["N", "C"] :
        l_at_connect_c, connectC = retrieveAtom.atomConnect(l_atom_lig , l_at_connect_N[1]["serial"])
        if connectC == ["C", "N", "N", "N"] :
            l_atom_out.append (deepcopy(l_at_connect_c[0])) # append C
            groupAtomN1, conect_N1 = retrieveAtom.atomConnect(l_atom_lig , int (l_at_connect_c[0]["connect"][1]))
            groupAtomN2, conect_N2 = retrieveAtom.atomConnect(l_atom_lig , int (l_at_connect_c[0]["connect"][2]))
            groupAtomN3, conect_N3 = retrieveAtom.atomConnect(l_atom_lig , int (l_at_connect_c[0]["connect"][3]))
            l_conect = [conect_N1, conect_N2, conect_N3]
            l_group_atom =  [groupAtomN1, groupAtomN2, groupAtomN3]
            
            i = 0
            while i < 3 :  
                if l_conect[i] == ["N", "C"] : 
                    nb_NH = nb_NH + 1
                    i = i + 1
                elif l_conect[i] == ["N", "C", "C"]:  
                    if l_group_atom[i][1]["serial"] != l_at_connect_c[0]["serial"] : 
                        l_atom_out.append (deepcopy(l_group_atom[i][1]))
                    elif l_group_atom[i][2]["serial"] != l_at_connect_c[0]["serial"] : 
                        l_atom_out.append (deepcopy(l_group_atom[i][2]))
                    else : 
                        print "ERROR l158"
                        return [0, []]
                    i = i + 1
                else :
                    return [0, []]
                
            # check number primary connectN1 -> case GAI not take, change ?
            if nb_NH == 2  : 
                l_atom_out.append (deepcopy(groupAtomN1[0]))
                l_atom_out.append (deepcopy(groupAtomN2[0]))
                l_atom_out.append (deepcopy(groupAtomN3[0]))
                
        else :
            return [0, []]
    
    # atom central structure
    elif connectN1 == ["N", "C", "C"] : 
        for at_conect in l_at_connect_N[1:] : 
            l_group_at, conect_N = retrieveAtom.atomConnect(l_atom_lig , at_conect["serial"])
            if conect_N == ["C", "N", "N", "N"] :
                l_c_central = l_group_at
                l_atom_out.append (deepcopy(l_c_central[0]))
                for group_at in l_group_at[1:] : 
                    l_goup_N, connect_N = retrieveAtom.atomConnect(l_atom_lig , group_at["serial"])
                    if connect_N == ["N", "C"] : 
                        nb_NH = nb_NH + 1
            else : 
                l_atom_out.append (deepcopy(l_group_at[0]))

        if nb_NH >= 2 and len (l_atom_out) >= 2 :
            
            l_atom_out.append (deepcopy(l_c_central[1]))
            l_atom_out.append (deepcopy(l_c_central[2]))
            l_atom_out.append (deepcopy(l_c_central[3]))
        else : 
            return [0, []]
        
    
    
    if len (l_atom_out) != 5 : 
        return [0,[]]
    else : 
        # append 2 carbons
        l_C, connectC = retrieveAtom.atomConnect(l_atom_lig , l_atom_out[0]["serial"])
        if connectC == ["C", "N", "N", "N"] : 
            l_serial = [l_atom_out[1]["serial"], l_atom_out[0]["serial"]]
            atomC1 = l_atom_out[1]
        else : 
            l_serial = [l_atom_out[0]["serial"], l_atom_out[1]["serial"]]
            atomC1 = l_atom_out[0]
        
        # distance between C2 and N
        d1 = calcul.distanceTwoatoms(atomC1, l_atom_out[2])
        d2 = calcul.distanceTwoatoms(atomC1, l_atom_out[3])
        d3 = calcul.distanceTwoatoms(atomC1, l_atom_out[4])
        l_dist = [d1, d2, d3]
        l_dist.sort ()
        
        for dist in l_dist : 
            if dist == d1 : 
                l_serial.append (l_atom_out[2]["serial"])
            elif dist == d2 : 
                l_serial.append (l_atom_out[3]["serial"])
            else : 
                l_serial.append (l_atom_out[4]["serial"])
        
        return [1,l_serial]
    
    
    return [0, []]        
# # 
# # def pyridine(listAtomConnectNitrogen, listAtomLigand):
# #     """search pyridine 
# #     in: list atom connected of nitrogen, list atom ligand
# #     out: boolean"""
# #         
# #     connectN1 = toolSubstructure.matrixElement(listAtomConnectNitrogen)
# #     
# #     if connectN1 == ["N", "C", "C"] : 
# #         nitrogenInit = listAtomConnectNitrogen[0]
# #         if  cycleOnlyTestCarbon(nitrogenInit["serial"], nitrogenInit["serial"], nitrogenInit["serial"], listAtomLigand, 6, 0) == 1 : 
# #             return 1
# #         
# #          
# # def diAmine (l_at_connect_N, listAtomLigand):
# #     """search diamine 
# #     in: list atom connected of nitrogen, list atom ligand
# #     out: boolean"""
# #         
# #     connect_element = toolSubstructure.matrixElement(l_at_connect_N)
# # #     print connect_element, "l199"
# #     if connect_element != ["N", "C"] :
# #         return 0
# #     else : 
# #         l_connect_C1, connect_C1 = retrieveAtom.atomConnect(listAtomLigand, l_at_connect_N[0]["connect"][1]) 
# # #         print connect_C1, "l204"
# #         if connect_C1 == ["C", "C", "N", "N"] or connect_C1 == ["C", "N", "C", "N"] or connect_C1 == ["C", "N", "N", "C"] : 
# #             for atom_connect in l_connect_C1[:1] :
# #                 l_connect_atom, connect_atom = retrieveAtom.atomConnect(listAtomLigand, atom_connect ["serial"])
# #                 
# #                 if atom_connect ["element"] == "N" : 
# #                     if connect_atom != ["N", "C"]  : 
# # #                         print "l211"
# #                         return 0
# #                 else  : 
# #                     if connect_atom == ["C", "C", "C", "C"] :
# # #                         print "l214" 
# #                         return 0
# #         else : 
# # #             print connect_C1
# #             return 0
# #     return 1
            


def acidCarboxylic(l_C2, l_atom_lig) : 
    
    l_serial = []
    connectO1 = toolSubstructure.matrixElement(l_C2)
    if connectO1 != ["O", "C"] : 
        return [0, []]
    else :
        l_serial.append (l_C2[0]["serial"]) 
        
        l_connect3, connectC2 = retrieveAtom.atomConnect(l_atom_lig, l_C2[1]["serial"])
        connectC2.sort()
        if connectC2 == ["C", "C", "O", "O"] :  
            l_serial.append (l_connect3[0]["serial"])
            
            for atom3 in l_connect3[1:] :
                if atom3["element"] == "O" and not atom3["serial"] in l_serial: 
                        l_O2, connectO2 = retrieveAtom.atomConnect(l_atom_lig, atom3["serial"])
                        if connectO2 == ["O", "C"] :
                            l_serial.append (atom3["serial"])
                        else : 
                            return [0, []]
                elif atom3["element"] == "C" and not atom3["serial"] in l_serial : 
                    l_atomC4, connectC4 = retrieveAtom.atomConnect(l_atom_lig, atom3["serial"])
#                     print connectC4, "Check"
                    connectC4_unique = sorted(set(connectC4),key=connectC4.index)
#                     print connectC4_unique
                    if connectC4_unique != ["C"] : 
                        return [0, []]
                    else :
                        for atom_conex in l_atomC4[1:] : 
                            if toolSubstructure.checkSingleBond(l_atomC4[0], atom_conex) == 0 : 
                                return [0, []]
                    
                    
                        
    # check 3 atom and not a hydroxyl group
    if len (l_serial) == 3 :
        return [1, l_serial]
    
    return [0, []]


    
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

def imidazole(l_atom2, l_atom_lig):
    """search imidazole global
    in: list atom connected of nitrogen, list atom ligand
    out: boolean"""
    
    
    connect1 = toolSubstructure.matrixElement(l_atom2)
#     print connect1
    # append first atom
    
    if connect1 == ["N", "C", "C"]:
        l_atom3, conect2 = retrieveAtom.atomConnect(l_atom_lig, l_atom2[1]["serial"])
        l_atom32, conect22 = retrieveAtom.atomConnect(l_atom_lig, l_atom2[2]["serial"])

        # sort list to reduce if condition
        conect2.sort ()
        conect22.sort ()
        
        if conect2 == ["C", "C", "N", "N"] or conect2 == ["C", "N", "N"]:
            l_serial_check = [l_atom2[0]["serial"], l_atom3[0]["serial"]]
            if imidazoleATOM3(l_atom3, l_serial_check, l_atom_lig)[0] == 1:
#                 print "l344"
                return [1, l_serial_check]

        elif conect22 == ["C", "C", "N", "N"] or conect22 == ["C", "N", "N"]:
            l_serial_check = [l_atom2[0]["serial"], l_atom32[0]["serial"]]
            if imidazoleATOM3(l_atom32, l_serial_check, l_atom_lig)[0] == 1:
#                 print "l350"
                return [1, l_serial_check]

    return [0, []]


def imidazoleATOM3(l_atom3, l_serial_check, l_atom_lig):
    """Check the atom3 3 in circle of imidazole
    in : atoms ligand in list, serial of first nitrogen, name of ligand
    out : boolean"""

    for atom3 in l_atom3:
        if atom3["element"] == "N":
            if not atom3["serial"] in l_serial_check : 
                l_atom4, conect3 = retrieveAtom.atomConnect(l_atom_lig, atom3["serial"])
                break

#     print conect3, "CN1, l367"
    if conect3 == ['N', 'C', 'C']:
        l_serial_check.append (l_atom4[0]["serial"])
#         print l_serial_check, "l370"
        if imidazoleATOM4(l_atom4, l_serial_check, l_atom_lig)[0] == 1:
            return [1, l_serial_check]

    return [0, l_serial_check]


def imidazoleATOM4(l_atom4, l_serial_check, l_atom_lig):
    """Check the atom4 4 in cercle of imidazole
    in : atoms ligand in list, serial of first nitrogen, name of ligand
    out : boolean"""
    
#     print l_serial_check, "l387"
    for atom4 in l_atom4:
        if atom4["element"] == "C":
            if not atom4["serial"] in l_serial_check : 
                l_atom5, conect4 = retrieveAtom.atomConnect(l_atom_lig, atom4["serial"])
                conect4.sort ()
                if conect4 == ["C", "C", "C", "N"] or conect4 == ["C", "C", "N"] : 
                    l_serial_check.append (l_atom5[0]["serial"])
                    imd, l_serial_check = imidazoleATOM5(l_atom5, l_serial_check, l_atom_lig)
#                     print l_serial_check, imd
                    if  imd == 1:
                        return [1, l_serial_check]
                    else :
                        l_serial_check.remove(l_atom5[0]["serial"])

    return [0, l_serial_check]
            

def imidazoleATOM5(l_atom5, l_serial_check, l_atom_lig):
    """Check the atom 5 in cercle of imidazole
    in : atoms ligand in list, serial of first nitrogen, name of ligand
    out : boolean"""
    
    
    for atom5 in l_atom5 : 
        if atom5["element"] == "C" and not atom5["serial"] in l_serial_check:
            l_atom6, conect5 =  retrieveAtom.atomConnect(l_atom_lig, atom5["serial"])
            # check connect atom4 IMD
            conect5.sort()
            if conect5 == ["C","C", "C", "N"] or conect5 == ["C","C", "N"] : 
                for atom6 in l_atom6 : 
                    if atom6["element"] == "N" and atom6["serial"] in l_serial_check :
                        l_serial_check.append (l_atom6[0]["serial"])
                        return [1, l_serial_check]
    return [0, l_serial_check]
    

######################################################################################################################


def globalSearch (dist_thresold, p_file_dataset,  pr_result, debug = 1):
    
    
    pr_summary = pr_result + "Sum/"
    
    if debug : print "Directory summary", pr_summary
    
    
    # load structure in summary ---> if use need place option one PDB by ligand
    d_neighbor = loadFile.loadCloseStruct (pr_summary)
    
    if d_neighbor != None : 
        if debug : 
            print "Type of structure stock", type (d_neighbor["global"])
        return d_neighbor
     
    
    #start, logFile = log.initAction("search neighbors in " +path.basename(p_file_dataset))
    l_lig = loadFile.resultFilterPDBLigand(p_file_dataset)
    nb_lig = len(l_lig)
    
    # ##Count Structure
    d_neighbor = {}
    l_neighbor_global = []

    
    # ##Write summary file
    d_files_summary = writeFile.openFileSummary(pr_summary)# sumary result
    
    # inialization  !!!!!!! 
     
    i = 0
    while i < nb_lig :
        if debug: print "Ligand: " + str(l_lig[i]["name"]) + " " + str(i) + " " + str(nb_lig)
        nb_PDB = len(l_lig[i]["PDB"])
        
        # take only one PDB by ligand
            
        j = 0
        while j < nb_PDB : 
            name_PDB = l_lig[i]["PDB"][j]
            l_atom_ligand = loadFile.ExtractInfoPDBID(l_lig[i]["PDB"][j])[l_lig[i]["name"]][0] # change not tested
            
            # search neighbor for every atom in ligand selected
            globalNeighbors(dist_thresold, l_atom_ligand, name_PDB, l_neighbor_global)
            # search neighbor for interest 
            interestGroup(dist_thresold, l_atom_ligand, name_PDB, d_neighbor)
            
            j = j + 1
        i = i + 1
            
    
    
    writeFile.neighborStruct(d_neighbor, l_neighbor_global, d_files_summary)
    writeFile.closeFileSummary(d_files_summary)
    
    # case where load directly substructure => why do not load directly in dictionnary
    d_neighbor["global"] = l_neighbor_global
    return d_neighbor



def SearchEnvironmentSaltBridgeProt(pr_result, l_PDB, max_dist, nb_lines = 150000, debug = 1):
    """
    This is for the intra-protein environment
    """
    
    pr_summary = pr_result + "Sum/"
    
    if debug == 1 : print "Directory summary", pr_summary
    
    
    # load structure in summary ---> if use need place option one PDB by ligand
    d_neighbor = loadFile.loadCloseStruct (pr_summary, nb_lines)
    
    if d_neighbor != None : 
        if debug : 
            print "Type of structure stock -> run building"
        return d_neighbor
     
    
    nb_PDB = len (l_PDB)
    
    # ##Write summary file
    d_files_summary = writeFile.openFileSummary(pr_summary, type_sum = "prot")# summary result
    
    i = 0
    while i < nb_PDB :
        print l_PDB[i], i
        InterestResForSaltBridge (d_files_summary, l_PDB[i], max_dist)
        
        i = i + 1
            
    writeFile.closeFileSummary(d_files_summary)
    
    # return non because files were created but not loaded in tempory memory
    return SearchEnvironmentSaltBridgeProt(pr_result, l_PDB, max_dist)
    
    

def InterestResForSaltBridge (d_filout, PDB_ID, max_distance):
    
    d_atom_interest = structure.DProtAtomSaltBridge()
    
    # double parsing - list of atom and dictionary
    l_atom_PDB = parsing.loadCoordSectionPDB(pathManage.pathDitrectoryPDB() + PDB_ID.lower() + ".pdb", "ATOM")
    d_res_PDB = parsing.BuildDicoRes(l_atom_PDB)
    
    
    for res_PDB in d_res_PDB.keys () :
        res = res_PDB.split ("_")[0]
        if  res in d_atom_interest.keys () :
            atom_and_neighbor = {} # control for write
            atom_and_neighbor[d_atom_interest[res]["subs"]] = {}
            l_atom_res = calcul.buildConnectMatrix(d_res_PDB[res_PDB])
            l_temp = []
            for atom in d_res_PDB[res_PDB] : 
                if atom["name"] in d_atom_interest[res]["atom"] : 
                    if len (d_atom_interest[res]["atom"]) == 1 : 
                        atom_and_neighbor[d_atom_interest[res]["subs"]] = [buildAtom(max_distance + structure.CalibrateDistanceNeighbor()[d_atom_interest[res]["subs"]], atom, PDB_ID, d_atom_interest[res]["subs"], l_atom_res)]
                    else : 
                        l_temp.append (atom)
                        if len (l_temp) == 2 : 
                            atom_and_neighbor[d_atom_interest[res]["subs"]] = [buildAtom(max_distance + structure.CalibrateDistanceNeighbor()[d_atom_interest[res]["subs"]], calcul.CenterPoint(l_temp[0], l_temp[1]), PDB_ID, d_atom_interest[res]["subs"], l_atom_res)]
            
              
            writeFile.neighborStruct(atom_and_neighbor, [], d_filout)        
                    
    
def interestGroup (max_distance, l_atom_lig, name_PDB, d_stock_neighbor, more_flex = 0):
    """Search different groups
    in : ligands in namePDB
    out : major nitrogen in the group of different structures
    change distance max for Guanidium + 1.5
    imidazole + 1.1
    and acid carboxylic + 1.3
    append more flex to GPCR study"""
    
    l_serialN = ListSerialElement(l_atom_lig, "N")
    l_serialO = ListSerialElement(l_atom_lig, "O")
    
    
    # different d_stock_neighbor
    for serialN in l_serialN:
        l_atom_connectN, conect = retrieveAtom.atomConnect(l_atom_lig, serialN)
        # check every substructure
        if imidazole(l_atom_connectN, l_atom_lig)[0] == 1:
            implementNeighborStruct (max_distance + structure.CalibrateDistanceNeighbor()["IMD"], l_atom_connectN, name_PDB, l_atom_lig, "IMD", d_stock_neighbor)
        elif Guanidium(l_atom_connectN, l_atom_lig)[0] == 1:
            implementNeighborStruct (max_distance + structure.CalibrateDistanceNeighbor()["GAI"], l_atom_connectN, name_PDB, l_atom_lig, "GAI", d_stock_neighbor)
#         elif diAmine(l_atom_connectN, l_atom_lig) == 1:
#             implementNeighborStruct (max_distance, l_atom_connectN, name_PDB, l_atom_lig, "Diamine", d_dia_temp)
#         elif pyridine(l_atom_connectN, l_atom_lig) == 1:
#             implementNeighborStruct (max_distance, l_atom_connectN, name_PDB, l_atom_lig, "Pyridine", d_stock_neighbor)
        elif cn(l_atom_connectN, l_atom_lig) == 1:
            implementNeighborStruct (max_distance, l_atom_connectN, name_PDB, l_atom_lig, "I", d_stock_neighbor)
        elif cnc(l_atom_connectN, l_atom_lig, more_flex = more_flex) == 1:
#             print "l562 -> search PDB"
            implementNeighborStruct (max_distance, l_atom_connectN, name_PDB, l_atom_lig, "II", d_stock_neighbor)
        elif cncc(l_atom_connectN, l_atom_lig, more_flex = more_flex) == 1:
            implementNeighborStruct (max_distance, l_atom_connectN, name_PDB, l_atom_lig, "III", d_stock_neighbor)
            
               
    for serialO in l_serialO :
        l_atom_connectO, conect = retrieveAtom.atomConnect(l_atom_lig, serialO)
        if acidCarboxylic(l_atom_connectO, l_atom_lig)[0] == 1:
            
            implementNeighborStruct (max_distance + structure.CalibrateDistanceNeighbor()["COO"], l_atom_connectO, name_PDB, l_atom_lig, "COO", d_stock_neighbor)
            

#######regroup neighbors case of imidazole, Guanidium and diamine###########


def implementNeighborStruct (max_distance, l_atom_connect_central, name_PDB, l_atom_lig, subs, st_neighbor):
    

    if subs == "IMD" : 
        atom_central = calcul.CenterImidazole (l_atom_connect_central, l_atom_lig)
    elif subs == "GAI" : 
        atom_central = calcul.CenterGuanidium (l_atom_connect_central, l_atom_lig)
    elif subs == "COO" : 
        atom_central = calcul.CenterAcidCarboxylic (l_atom_connect_central, l_atom_lig)
    else : 
        atom_central = l_atom_connect_central[0]
    
#     print "Test atom considered", atom_central
    atom_neighbors = buildAtom(max_distance, atom_central, name_PDB, subs, l_atom_lig)
    
    # dynamic implementation    
    if not subs in st_neighbor.keys () :
        st_neighbor[subs] = []
    if not atom_neighbors in st_neighbor[subs] : 
        st_neighbor[subs].append(atom_neighbors)



######################################################################################


def buildAtom(rayon, at_central, PDB, subs, l_atom_lig):
    """Build count structure for stAtom structure
    in: distance, serial atom nitrogen, type study, list atom ligand
    out: count strucutre"""
    
    atom = {}
    atom["PDB"] = PDB
    atom["resName"] = at_central["resName"]
    atom["serial"] = at_central["serial"]
    atom["x"] = at_central["x"]
    atom["y"] = at_central["y"]
    atom["z"] = at_central["z"]
    
    atom["neighbors"] = neighbors(rayon, at_central, PDB, subs, l_atom_lig)

    return atom



def neighbors(rayon, atom_central, pdb, subs = "global", l_atom_lig = [] ): # change the name because in the same time function and variable
    """Search neighbors for all ligand
    in : rayon where is atoms, central atom, pdb file
    out : list atoms found"""

    l_atom = []
    linesPDB = loadFile.openPdbFile(pdb)
    for line in linesPDB:
        if search("^ATOM", line) or search("^HETATM", line): 
            atom = parsing.lineCoords(line)
            if atom != {} and atom["element"] != "H":
                distance = calcul.distanceTwoatoms(atom_central, atom)
                if distance <= rayon and distance != 0.0:
                    if atom_central["resSeq"] != atom["resSeq"]: # check if variation
                        if tool.atomInList(l_atom, atom) == 0:
                            atom["distance"] = distance
                            atom["angleSubs"] = calcul.angleSubs(atom_central, atom, l_atom_lig, subs)
                            atom["classification"] = structure.classificationATOM(atom)
                            l_atom.append(atom)

    return l_atom


def repetitionInPDB(ligand, pdb):
    """Count number of same ligand in pdb file
    in : name of ligand, name of pdb file
    out : number of ligand in pdb"""

    groupAtom = loadFile.ExtractInfoPDBID(pdb)[ligand][0] # change not tested

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



def ListSerialElement(l_atom, element):
    """Search N in the atoms list
    return a list of number atom"""

    l_out = []
    for atom in l_atom:
        if atom["element"] == element:
            serial = atom["serial"]
            if not serial in l_out:
                l_out.append(atom["serial"])
    return  l_out






def checkAngleInSearchNeighbor(atomRetrieve, subs):
    """Check for every neighbors of atom retrieve the angleSubs
    in: atomRetrieve
    out: atomRetrieve modified"""
    
    d_angle_limit = structure.criteraAngle(subs)
    nbNeighbors = len(atomRetrieve["neighbors"])
        
    i = 0
    while i < nbNeighbors :
        # print atomRetrieve["neighbors"]
        
        if checkListAngle(atomRetrieve["neighbors"][i]["angleSubs"], d_angle_limit) == 0 : 
            del atomRetrieve["neighbors"][i]
            nbNeighbors = nbNeighbors - 1
            continue
        else : 
            i = i + 1
    

def checkListAngle (l_angle, d_limit): 
    
    
    for angleSubs in l_angle : 
        if angleSubs < d_limit["INF"] : 
            return 0
        if angleSubs > d_limit["SUP"] : 
            return 0
    
    return 1
        
        
    
    
        
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

################ Cyles #######################
    
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


#################### other fonction search ################"


def Nclose (atom, l_at_lig) : 
    
    d_min = 100
    for at_lig in l_at_lig : 
        if at_lig["element"] != "N" : 
            continue
        else : 
            d = calcul.distanceTwoatoms(atom, at_lig)
            if d < d_min : 
                d_min = d
                atom_temp = deepcopy(at_lig)
    
    if "atom_temp" in locals() : 
        return atom_temp
    else : 
        return {}



       
            
