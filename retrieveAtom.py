import tool
from copy import deepcopy
import calcul


def serial(at_serial, list_atom) :
    """For a atom id search atom dictionnary
    in: at_serial atom and list atoms
    out: type of atom"""

    #print list_atom
    for atom in list_atom :
        if int(atom["serial"]) == int(at_serial) :
            return atom

    return 0



def atomConnect (list_atom_ligand, serial_atom):
    """For a serial atom found the connect atoms
    in: list atom ligand, serial atom
    out: list atoms connect and connect matrix element"""

    atomN = serial(serial_atom, list_atom_ligand)
    conect = []
    listAtom = []

    if atomN == 0:
        return listAtom, conect

    matrixConect = atomN["connect"]
       

    for elementMatrixConect in matrixConect:
        atom = (serial(int(elementMatrixConect), list_atom_ligand))
        if atom != 0:
            listAtom.append(atom)
            conect.append(atom["element"])
        else :
            conect.append("out")  ########################################### case out ligand
            
    tool.dellH(conect)
    return listAtom, conect

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


def substructure (substruct, serial_at_central, l_at) : 
    
    if substruct == "Primary" or substruct == "Secondary" or substruct == "Tertiary" : 
        l_out = atomConnect(l_at, serial_at_central)
        return l_out[0]
    
    if substruct == "Imidazole" : 
        l_out = []
        l_n = atomConnect(l_at, serial_at_central)
        
        
    
    
    


 
