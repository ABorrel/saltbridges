import retrieveAtom
import calcul


def checkConectOnlyC(atom, listAtom) :
    """Check if atom is connect N is only H and C 
    in: atom checked, list of atom
    out: 0 or 1"""

    connectMatrixElement = connectMatrixSerialToConnectMatrixElement(atom["connect"], listAtom)
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
#     print distance
    if distance >= d_min:
        return 1
    else:
        return 0


def checkCoplanar (NAtom, l_atom):
    """Check if tertiary amine is coplanar
    in: atom nitrogen
    out: list atom"""
    
    try: distance = calcul.coplanar(NAtom, l_atom)
    except: return 0
    
#     print distance

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

def matrixElement(listAtom):
    """For list of atom retrieve list element
    in: list Atom
    out: list Element"""
    
    out = []
    for atom in listAtom : 
        try :
            out.append(atom["element"])
        except :
            continue
    
    return out
