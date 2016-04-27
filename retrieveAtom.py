import tool
from copy import deepcopy
import searchPDB


def serial(at_serial, list_atom) :
    """For a atom id search atom dictionnary
    in: at_serial atom and list atoms
    out: type of atom"""

    #print list_atom
    for atom in list_atom :
        if int(atom["serial"]) == int(at_serial) :
            return atom

    return 0



def atomConnect (l_atom_lig, serial_atom):
    """For a serial atom found the connect atoms
    in: list atom ligand, serial atom
    out: list atoms connect and connect matrix element"""

    atom_in = serial(serial_atom, l_atom_lig)
    conect = []
    l_atom_out = []

    if atom_in == 0:
        return l_atom_out, conect

    matrixConect = atom_in["connect"]
       

    for elementMatrixConect in matrixConect:
        atom = (serial(int(elementMatrixConect), l_atom_lig))
        if atom != 0:
            l_atom_out.append(atom)
            conect.append(atom["element"])
        # remove else condition because H are removed before and give -> out
#         else :
#             conect.append("out")
            
    tool.dellH(conect)
    return l_atom_out, conect

###################serine protease################

def cycle(serialAtomInit, listAtomLigand):
    """retrieve cycle in list of atom
    in: serial atom in cycle
    out: list atoms"""
    
    listOUT = []
    listOUT.append(serialAtomInit)
    
    flagCycle = 0
    atomList = [serialAtomInit]
    
    while flagCycle != 1 :
        flagCycle = 1
        atomList2 = []
        for atomcheck in atomList :
            for connect in atomcheck["connect"] :
                atom = serial(connect, listAtomLigand)
                if atom != 0 : 
                    if not atom in listOUT : 
                        if atom["cycle"] == 1 :
                            if not atom in listOUT :  
                                listOUT.append(atom)
                            flagCycle = 0
                            for serialAtomConnectCycle in atom["connect"][1:] : 
                                atomConnectCycle = serial(serialAtomConnectCycle, listAtomLigand)
                                if atomConnectCycle != 0 :  
                                    if atomConnectCycle ["cycle"] == 1 : 
                                        if not atomConnectCycle in atomList2 : 
                                            atomList2.append(atomConnectCycle)
                                    else : 
                                        if not atomConnectCycle in listOUT : 
                                            listOUT.append(atomConnectCycle)
                        else : 
                            if not atom in listOUT : 
                                listOUT.append(atom)
        
        atomList = deepcopy(atomList2)
        
    return listOUT
    
###################################################################
###################################################################


def substructure (subs, at_central, l_at_lig) : 
    
    if subs == "I" or subs == "II" or subs == "III" : 
        serial_atom = at_central["serial"]
        l_sub = atomConnect(l_at_lig, serial_atom)
        return l_sub[0]
    
    elif subs == "IMD" : 
        atom_N = searchPDB.Nclose (at_central, l_at_lig)
        l_connectN, conect = atomConnect(l_at_lig, atom_N["serial"])
        out_imd = searchPDB.imidazole(l_connectN, l_at_lig)
        
        if out_imd[0] == 0 : 
            print "ERROR, retrieveAtom: l-104"
            return []
        else : 
            l_atom_IMD = []
            for serial_atom in out_imd[1] : 
                l_atom_IMD.append (serial(serial_atom, l_at_lig))
            return l_atom_IMD
    
    elif subs == "GAI" : 
        l_out = []
        l_atomC, connectC = atomConnect(l_at_lig, at_central["serial"])
        l_atom_N, connectN = atomConnect(l_at_lig, l_atomC[1]["serial"])
        out_gua = searchPDB.Guanidium(l_atom_N, l_at_lig)
        
        if out_gua[0] == 0 : 
            print "ERROR"
            return []
        out_search = out_gua[1]
        
        for serial_at in out_search : 
#             print serial_at
            l_out.append (serial(serial_at, l_at_lig))
        return l_out
    
    elif subs == "COO" : 
         
        l_atomC, connectC = atomConnect(l_at_lig, at_central["serial"])
        if connectC[1] == "O" : 
            l_atom_O, connectO = atomConnect(l_at_lig, l_atomC[1]["serial"])
        else : 
            l_atom_O, connectO = atomConnect(l_at_lig, l_atomC[2]["serial"])  
         
        l_out = []
        out_search = searchPDB.acidCarboxylic(l_atom_O, l_at_lig)
        
#         print out_search[0], len (out_search[1])
        if out_search[0] == 0 : 
            print "ERROR"
            return []
        l_serial = out_search[1]
        
        for serial_at in l_serial : 
            l_out.append (serial(serial_at, l_at_lig))
        return l_out
    
    


 
