import retrieveAtom
import calcul


def checkConectOnlyC(atom, l_atom_lig) :
    """Check if atom is connect N is only H and C 
    in: atom checked, list of atom
    out: 0 or 1"""

    connectMatrixElement = connectMatrixSerialToConnectMatrixElement(atom["connect"], l_atom_lig)
    flag = 0    
    
    for element in connectMatrixElement [1:]: 
        if element != "C" : 
            if element != "H" : 
                flag = flag + 1 
                
    if flag > 1 : 
        return 0
    return 1 


def connectMatrixSerialToConnectMatrixElement(connectMatrixSerial, listAtom):
    '''transform Connect matrix with serial 
    in: connect matrix serial, list atom
    out: list -> matrix connect element'''

    link = []
    for serial in connectMatrixSerial:
        try : link.append(retrieveAtom.serial(int(serial), listAtom)["element"])
        except : pass
    return link


def checkSingleBond(atom1, atom2, d_min = 1.42):
    """check distance between two atoms (1.42)
    in: atom1 and atom2
    out: 0 or 1"""
    
    distance = calcul.distanceTwoatoms(atom1, atom2)
    if distance >= d_min:
        return 1
    else:
        return 0


def checkCoplanar (Natom, l_atom, debug = 0):
    """Check if tertiary amine is coplanar
    in: atom nitrogen
    out: list atom"""
    
    try: distance = calcul.coplanar(Natom, l_atom)
    except: return 0
    
    if debug == 1 : print distance, "planarity"

    if distance >= 1.0:
        return 1
    else:
        return 0


def retrieveAtomConnected(atom, listAtom):
    """retrieve list atom connected
    in: atom, list atoms
    out: list atoms connected"""

    listAtom = []
    atom = retrieveAtom.serial(int(atom), listAtom)
    connectMatrix = atom['connect']
    for i in connectMatrix:
        listAtom.append(retrieveAtom.serial(i, listAtom))
    return listAtom

def matrixElement(l_atom):
    """For list of atom retrieve list element
    in: list Atom
    out: list Element"""
    
    out = []
    for atom in l_atom : 
        try :
            out.append(atom["element"])
        except :
            continue
    
    return out
