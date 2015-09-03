from math import acos, asin, cos, sin, degrees
from math import sqrt
from re import search
from copy import deepcopy

import loadFile
import parsing
import retrieveAtom
from sympy import solve
from sympy import symbols


def distanceTwoatoms(atom1, atom2):  ##############to review
    '''calculate distance of 2 atoms
    in : - atom1 structure
         - atom2 structure
    out : distance -> float
          100 if impossible calcul'''

    try:
        x1 = float(atom1['x'])
        x2 = float(atom2['x'])
        xd = x2 - x1

        y1 = float(atom1['y'])
        y2 = float(atom2['y'])
        yd = y2 - y1

        z1 = float(atom1['z'])
        z2 = float(atom2['z'])
        zd = z2 - z1

        return sqrt(xd * xd + yd * yd + zd * zd)
    except:
        return 100


def coplanar (atom, l_atom_ligand):
    """Calculate the orthogonal distance between nitrogen atom of tertiary amine and plan with 3 carbons connect
    in : - atom of nitrogen -> dictionnary
         - all atom of ligand -> list of atom dictionnary
    out : - distance -> float"""

    matrix = atom["connect"]
    if len(matrix) < 4:
        print "Atom does not 3 bonds !!"
        return

    else:
        point = [float(retrieveAtom.serial(matrix[0], l_atom_ligand)["x"]), float(retrieveAtom.serial(matrix[0], l_atom_ligand)["y"]), float(retrieveAtom.serial(matrix[0], l_atom_ligand)["z"])]

        d = symbols("d")
        
        atom1 = retrieveAtom.serial(matrix[1], l_atom_ligand)
        atom2 = retrieveAtom.serial(matrix[2], l_atom_ligand)
        atom3 = retrieveAtom.serial(matrix[3], l_atom_ligand)
        
        point1 = [float(atom1["x"]), float(atom1["y"]), float(atom1["z"])]
        point2 = [float(atom2["x"]), float(atom2["y"]), float(atom2["z"])]
        point3 = [float(atom3["x"]), float(atom3["y"]), float(atom3["z"])]


        v1 = [point3[0] - point1[0], point3[1] - point1[1], point3[2] - point1[2]]
        v2 = [point2[0] - point1[0], point2[1] - point1[1], point2[2] - point1[2]]
        a = v1[1] * v2[2] - v1[2] * v2[1]
        b = v1[2] * v2[0] - v1[0] * v2[2]
        c = v1[0] * v2[1] - v1[1] * v2[0]

        d = float(solve(a * point1[0] + b * point1[1] + c * point1[2] - d, d)[0])

        nominator = abs(a * point[0] + b * point[1] + c * point[2] - d)
        denominator = (a ** 2 + b ** 2 + c ** 2) ** (1 / 2)

        distance = nominator / denominator

        return distance


def coplanarPoint (atom, l_atom):
    """Calculate the orthogonal distance between nitrogen atom of tertiary amine and plan with 3 carbons connect
    in : - atom of nitrogen -> dictionnary
         - all atom of ligand -> list of atom dictionnary
    out : - distance -> float"""

    if len(l_atom) != 3:
        print "Atom does not 3 bonds !!"
        return

    else:
        point = [float(atom["x"]), float(atom["y"]), float(atom["z"])]

        d = symbols("d")
        
        
        point1 = [float(l_atom[0]["x"]), float(l_atom[0]["y"]), float(l_atom[0]["z"])]
        point2 = [float(l_atom[1]["x"]), float(l_atom[1]["y"]), float(l_atom[1]["z"])]
        point3 = [float(l_atom[2]["x"]), float(l_atom[2]["y"]), float(l_atom[2]["z"])]


        v1 = [point3[0] - point1[0], point3[1] - point1[1], point3[2] - point1[2]]
        v2 = [point2[0] - point1[0], point2[1] - point1[1], point2[2] - point1[2]]
        a = v1[1] * v2[2] - v1[2] * v2[1]
        b = v1[2] * v2[0] - v1[0] * v2[2]
        c = v1[0] * v2[1] - v1[1] * v2[0]

        d = float(solve(a * point1[0] + b * point1[1] + c * point1[2] - d, d)[0])

        nominator = abs(a * point[0] + b * point[1] + c * point[2] - d)
        denominator = (a ** 2 + b ** 2 + c ** 2) ** (1 / 2)

        distance = nominator / denominator

        return distance



def buildConnectMatrix(l_atom_lig, namePDB = ""):
    """Building the matrix connect. threshold -> 1.6 angstroms
    in : - Atoms of ligands -> list of atoms -> dictionnary
         - PDB file associted of ligand
    out : - Ligand with connect matrix"""


#     linesPDB = loadFile.openPdbFile(namePDB)
#     groupAtomPDB = []
#     for line in linesPDB:
#         if search("^ATOM", line) or search("^HETATM", line):
#             groupAtomPDB.append(parsing.lineCoords(line))

    # initialization of connect matrix
    for atomLigand in l_atom_lig:
        atomLigand["connect"] = []

    for atomLigand in l_atom_lig:
        atomLigand["connect"].append(atomLigand["serial"])
        for atomPDB in l_atom_lig:
            distance = distanceTwoatoms(atomLigand, atomPDB)
            if distance < 1.8 and distance != 0.1 and distance != "ERROR":
                if not atomPDB["serial"] in atomLigand["connect"]:
                    atomLigand["connect"].append(atomPDB["serial"])

    return l_atom_lig


def equationPlan (A, B, C):
    
#     print A['x'], B['x'], C['x']
#     print A['y'], B['y'], C['y']
#     print A['z'], B['z'], C['z']
    a = (B['y'] - A['y']) * (C['z'] - A['z']) - (B ['z'] - A['z']) * (C['y'] - A['y'])
    b = -((B['x'] - A['x']) * (C['z'] - A['z']) - (B ['z'] - A['z']) * (C ['x'] - A['x']))
    c = (B ['x'] - A['x']) * (C['y'] - A['y']) - (B['y'] - A['y']) * (C['x'] - A['x'])
    d = -(a * A['x'] + b * A['y'] + c * A['z'])
    
    return [a, b, c, d]
    

def vectoriel (N, C1):
    
    
    NC1x = C1["x"] - N["x"]
    NC1y = C1["y"] - N["y"]
    NC1z = C1["z"] - N["z"]
    
    return [NC1x, NC1y, NC1z]


def normeVectoriel (listVectoriel1, listVectoriel2):
    
    terme1 = pow(listVectoriel1[1] * listVectoriel2[2] - listVectoriel1[2] * listVectoriel2[1], 2)
    terme2 = pow(listVectoriel1[2] * listVectoriel2[0] - listVectoriel1[0] * listVectoriel2[2], 2)
    terme3 = pow(listVectoriel1[0] * listVectoriel2[1] - listVectoriel1[1] * listVectoriel2[0], 2)
    
    return sqrt (terme1 + terme2 + terme3)

def normeVector (point1, point2):
    
    terme1 = pow(point1["x"] - point2["x"], 2)
    terme2 = pow(point1["y"] - point2["y"], 2)
    terme3 = pow(point1["z"] - point2["z"], 2)
    
    return sqrt(terme1 + terme2 + terme3)


def scalar(vector1, vector2):
    
    x = vector1[0] * vector2[0]
    y = vector1[1] * vector2[1]
    z = vector1[2] * vector2[2]

    return x + y + z


def angleVector(pointD, pointCentral, pointG):
    
    vectorNC1 = vectoriel(pointCentral, pointD)
    vectorNC2 = vectoriel(pointCentral, pointG)
    normeNC1 = normeVector(pointCentral, pointD)
    normeNC2 = normeVector(pointCentral, pointG)

    scalarNC1NC2 = scalar(vectorNC1, vectorNC2) 

    try :
        alpha = degrees(acos(scalarNC1NC2 / (normeNC1 * normeNC2)))
        return alpha
    except :
        return "ERROR"
    

    

def checkPoint(pointTest, pointN, pointC1, pointC2, plan, alpha):
    Eqplan = pointTest["x"] * plan[0] + pointTest["y"] * plan[1] + pointTest["z"] * plan[2] + plan[3]
    angleRef = (360 - alpha) / 2 - 5
    
    
    
    angleC1NTest = angleVector(pointC1, pointN, pointTest)
    angleC2NTest = angleVector(pointC2, pointN, pointTest)
    
    if angleC1NTest > angleRef and angleC2NTest > angleRef : 
        diffAngle = abs(angleC1NTest - angleC2NTest)
        return  diffAngle + abs(Eqplan) * 10
    else :
        return 1000






def angleTertiaryAmine (atomNitrogen, atomCounterIon, listAtomLigand):
    """Calcul for nitrogen and counter ion angles
    in: atom nitrogen, atom Counter Ions
    out: list of 3 angles"""
    
    matrix = atomNitrogen["connect"]
    if len(matrix) < 4:
        print "Atom does not 3 bonds !!"
        return

    else:

        carbon1 = retrieveAtom.serial(matrix[1], listAtomLigand)
        carbon2 = retrieveAtom.serial(matrix[2], listAtomLigand)
        carbon3 = retrieveAtom.serial(matrix[3], listAtomLigand)

        distanceNitrogenCounterIon = distanceTwoatoms(atomNitrogen, atomCounterIon)

        distanceCarbon1CounterIon = distanceTwoatoms(atomCounterIon, carbon1)
        distanceCarbon2CounterIon = distanceTwoatoms(atomCounterIon, carbon2)
        distanceCarbon3CounterIon = distanceTwoatoms(atomCounterIon, carbon3)

        distanceNitrogenCarbon1 = distanceTwoatoms(atomNitrogen, carbon1)
        distanceNitrogenCarbon2 = distanceTwoatoms(atomNitrogen, carbon2)
        distanceNitrogenCarbon3 = distanceTwoatoms(atomNitrogen, carbon3)


        angle1 = degrees(acos((distanceNitrogenCarbon1 * distanceNitrogenCarbon1 + distanceNitrogenCounterIon * distanceNitrogenCounterIon - distanceCarbon1CounterIon * distanceCarbon1CounterIon) / (2 * distanceNitrogenCarbon1 * distanceNitrogenCounterIon)))
        angle2 = degrees(acos((distanceNitrogenCarbon2 * distanceNitrogenCarbon2 + distanceNitrogenCounterIon * distanceNitrogenCounterIon - distanceCarbon2CounterIon * distanceCarbon2CounterIon) / (2 * distanceNitrogenCarbon2 * distanceNitrogenCounterIon)))
        angle3 = degrees(acos((distanceNitrogenCarbon3 * distanceNitrogenCarbon3 + distanceNitrogenCounterIon * distanceNitrogenCounterIon - distanceCarbon3CounterIon * distanceCarbon3CounterIon) / (2 * distanceNitrogenCarbon3 * distanceNitrogenCounterIon)))

        return [angle1, angle2, angle3]


def AngleSecondaryFromLig(atomN, atom_neighbor, l_atom_lig):
    """Calcul for nitrogen and neighbor
    in: atom nitrogen, atom neighbor
    out: angle between muddle of 2 cabones and N and neigbbor"""
    
    l_connect = atomN["connect"]
    if len(l_connect) < 3:
        print "Atom does not 2 bonds !!"
        return

    else:

        atomC1 = retrieveAtom.serial(l_connect[1], l_atom_lig)
        atomC2 = retrieveAtom.serial(l_connect[2], l_atom_lig)
        
        angle1 = Angle3Atoms (atomC1, atomN, atom_neighbor)
        angle2 = Angle3Atoms (atomC2, atomN, atom_neighbor)
        
        return [angle1, angle2]


def Angle3Atoms (atom1, atom_center, atom3):
    
    dist_central_atom3 = distanceTwoatoms(atom_center, atom3)
    dist_central_atom1 = distanceTwoatoms(atom_center, atom1)
    dist_atom1_atom3 = distanceTwoatoms(atom3, atom1)
        
    angle_out = degrees(acos((dist_central_atom1 * dist_central_atom1 + dist_central_atom3 * dist_central_atom3 - dist_atom1_atom3 * dist_atom1_atom3) / (2 * dist_central_atom1 * dist_central_atom3)))
        
    return angle_out
    


def anglePrimaryAmine(atomNitrogen, atomCounterIon, atomLigands):
    """Calcul for nitrogen and counter ion angles
    in: atom nitrogen, atom Counter Ions
    out: list of 2 angles"""
    
    matrix = atomNitrogen["connect"]
    if len(matrix) < 2:
        print "Atom does not 1 bonds !!"
        return

    else:

        carbon1 = retrieveAtom.serial(matrix[1], atomLigands)
        distanceNitrogenCounterIon = distanceTwoatoms(atomNitrogen, atomCounterIon)
        distanceCarbon1CounterIon = distanceTwoatoms(atomCounterIon, carbon1)
        distanceNitrogenCarbon1 = distanceTwoatoms(atomNitrogen, carbon1)

        angle1 = degrees(acos((distanceNitrogenCarbon1 * distanceNitrogenCarbon1 + distanceNitrogenCounterIon * distanceNitrogenCounterIon - distanceCarbon1CounterIon * distanceCarbon1CounterIon) / (2 * distanceNitrogenCarbon1 * distanceNitrogenCounterIon)))

        return [angle1]



def angleImidazole(atom_central, atom_check, l_atom_lig, debug = 0):
    """
    based one central atom, take C between the two N
    """
    
    for atom_lig in l_atom_lig : 
        if atom_lig["element"] == "C" :
            d_temp = distanceTwoatoms(atom_central, atom_lig)
            if d_temp < 2.0 : 
                d_Cconsidered = deepcopy(atom_lig)
                return [angleVector(d_Cconsidered, atom_central, atom_check)]
    return ["NA"]
    




def angleSubs(central_atom, atom_found, l_atom_lig, subs):
    """calcul angleSubs for each structure study
    in: atom nitrogen, list atoms ligand, structure study
    out: list of angles"""
    
    if subs == "Primary" :
        return anglePrimaryAmine(central_atom, atom_found, l_atom_lig)
    elif subs == "Secondary" :
        return AngleSecondaryFromLig(central_atom, atom_found, l_atom_lig)
    elif subs == "Tertiary" :
        return angleTertiaryAmine(central_atom, atom_found, l_atom_lig)
    elif subs == "Imidazole" :
        return angleImidazole(central_atom, atom_found, l_atom_lig)
#     elif subs == "Diamine" :
#         return anglePrimaryAmine(central_atom, atom_found, l_atom_lig)
#     elif subs == "Pyridine" :
#         return AngleSecondaryFromLig(central_atom, atom_found, l_atom_lig)
    elif subs == "Guanidium" :
        return angleGuanidium(central_atom, atom_found, l_atom_lig)
    elif subs == "AcidCarboxylic" :
        return angleAcidCarboxylic(central_atom, atom_found, l_atom_lig)
    
    else :
        return []
    



def angleAcidCarboxylic(central_atom, atom_check, l_atom_lig) : 
    
    l_atom_connect, matrix_connect = retrieveAtom.atomConnect(l_atom_lig, central_atom["serial"])
    l_O_temp = []
    
    for atom_connect in l_atom_connect[1:] : 
        l_temp_connect, connect_matrix = retrieveAtom.atomConnect(l_atom_lig, atom_connect["serial"])
        if connect_matrix == ["O", "C"] : 
            l_O_temp.append (atom_connect)
    
    
    at_considered = CenterPoint(l_O_temp[0], l_O_temp[1])    
        
    
    return [angleVector(central_atom, at_considered, atom_check)]



def angleGuanidium(central_atom, atom_check, l_atom_lig) : 
    
    l_N_temp = []
    
    l_atom_connect_central, l_atom_element = retrieveAtom.atomConnect(l_atom_lig, central_atom["serial"])
    
    for atom_N in l_atom_connect_central[1:] : 
        l_atom_connect_N, l_element_N = retrieveAtom.atomConnect(l_atom_lig, atom_N["serial"])
        
        
        if l_element_N == ["N", "C"] : 
            if not l_atom_connect_N[0] in l_N_temp : 
                l_N_temp.append (deepcopy(l_atom_connect_N[0]))


    if len(l_N_temp) == 2 :
        at_considered = CenterPoint(l_N_temp[0], l_N_temp[-1])

    elif len(l_N_temp) > 2 : 
        for N_temp in l_N_temp : 
            if distanceTwoatoms(central_atom, N_temp) > 3.0 : 
                del N_temp 
        if len(l_N_temp) == 2 :
            at_considered = CenterPoint(l_N_temp[0], l_N_temp[-1])
    else : 
        return "ERROR"
        
    return [angleVector(central_atom, at_considered, atom_check)]
    
    



def anglePrimaryAmineCalculVol(atomN, atomC, atomTest):
    """calcul for three atoms angle for primary amine
    in: atom nitrogen, atom carbon and atom test
    out: list angle"""
    
    distanceNitrogenCounterIon = distanceTwoatoms(atomN, atomTest)
    distanceCarbon1CounterIon = distanceTwoatoms(atomTest, atomC)
    distanceNitrogenCarbon1 = distanceTwoatoms(atomN, atomC)

    try : 
        angle1 = degrees(acos((distanceNitrogenCarbon1 * distanceNitrogenCarbon1 + distanceNitrogenCounterIon * distanceNitrogenCounterIon - distanceCarbon1CounterIon * distanceCarbon1CounterIon) / (2 * distanceNitrogenCarbon1 * distanceNitrogenCounterIon)))
        return [angle1]
    except : return [0.00]


# def angleSecondaryAmineCalculVol(atomN, atomC1, atomC2, atom_test):
#     """calcul for three atoms angleVector for secondary amine
#     in: atom nitrogen, atoms carbons and atom test
#     out: list angles"""
# 
#     distanceNitrogenTest = distanceTwoatoms(atomN, atom_test)
# 
#     distanceCarbon1Test = distanceTwoatoms(atom_test, atomC1)
#     distanceCarbon2Test = distanceTwoatoms(atom_test, atomC2)
# 
#     distanceNitrogenCarbon1 = distanceTwoatoms(atomN, atomC1)
#     distanceNitrogenCarbon2 = distanceTwoatoms(atomN, atomC2)
# 
#     angle1 = degrees(acos((distanceNitrogenCarbon1 * distanceNitrogenCarbon1 + distanceNitrogenTest * distanceNitrogenTest - distanceCarbon1Test * distanceCarbon1Test) / (2 * distanceNitrogenCarbon1 * distanceNitrogenTest)))
#     angle2 = degrees(acos((distanceNitrogenCarbon2 * distanceNitrogenCarbon2 + distanceNitrogenTest * distanceNitrogenTest - distanceCarbon2Test * distanceCarbon2Test) / (2 * distanceNitrogenCarbon2 * distanceNitrogenTest)))
# 
# 
#     return [angle1, angle2]    
    
    
    
def angleTertiaryAmineCalculVol (atomNitrogen, atomTest, atomC1, atomC2, atomC3):
    """calcul for three atoms angleVector for tertiary amine
    in: atom nitrogen, atoms carbons and atom test
    out: list angles"""
       
    distanceNitrogenCounterIon = distanceTwoatoms(atomNitrogen, atomTest)

    distanceCarbon1CounterIon = distanceTwoatoms(atomTest, atomC1)
    distanceCarbon2CounterIon = distanceTwoatoms(atomTest, atomC2)
    distanceCarbon3CounterIon = distanceTwoatoms(atomTest, atomC3)

    distanceNitrogenCarbon1 = distanceTwoatoms(atomNitrogen, atomC1)
    distanceNitrogenCarbon2 = distanceTwoatoms(atomNitrogen, atomC2)
    distanceNitrogenCarbon3 = distanceTwoatoms(atomNitrogen, atomC3)


    angle1 = degrees(acos((distanceNitrogenCarbon1 * distanceNitrogenCarbon1 + distanceNitrogenCounterIon * distanceNitrogenCounterIon - distanceCarbon1CounterIon * distanceCarbon1CounterIon) / (2 * distanceNitrogenCarbon1 * distanceNitrogenCounterIon)))
    angle2 = degrees(acos((distanceNitrogenCarbon2 * distanceNitrogenCarbon2 + distanceNitrogenCounterIon * distanceNitrogenCounterIon - distanceCarbon2CounterIon * distanceCarbon2CounterIon) / (2 * distanceNitrogenCarbon2 * distanceNitrogenCounterIon)))
    angle3 = degrees(acos((distanceNitrogenCarbon3 * distanceNitrogenCarbon3 + distanceNitrogenCounterIon * distanceNitrogenCounterIon - distanceCarbon3CounterIon * distanceCarbon3CounterIon) / (2 * distanceNitrogenCarbon3 * distanceNitrogenCounterIon)))

    return [angle1, angle2, angle3]



def angleImidazoleCalculVol(atomN1, atomN3, atom_test):


    atom_center = {}
    atom_center["x"] = (atomN1["x"] + atomN3["x"]) / 2
    atom_center["y"] = (atomN1["y"] + atomN3["y"]) / 2
    atom_center["z"] = (atomN1["z"] + atomN3["z"]) / 2
    
    
    return anglePrimaryAmineCalculVol(atom_center, atomN1, atom_test)
    
  
    
def CenterImidazole (l_atom_connectN, l_atom_lig) : 
    
    atomN1 = l_atom_connectN[0]
    
    
    for atom_connectN in l_atom_connectN[1:] :
        l_atom_connect2, l_serial = retrieveAtom.atomConnect(l_atom_lig, atom_connectN["serial"])
        
        for atom_connect2 in l_atom_connect2[1:] : 
            if atom_connect2["element"] == "N" : 
                d_N = distanceTwoatoms(atomN1, atom_connect2)
                if d_N < 2.4 : 
                    return CenterPoint(atomN1, atom_connect2)
                
    print "ERROR"            
    return "ERROR"   
    

def CenterPoint (atom1, atom2):

    a_out = parsing.EmptyAtom()
    
    a_out["serial"] = 0
    a_out["name"] = atom1["name"]
    a_out["char"] = atom1["char"]
    a_out["resName"] = atom1["resName"]
    try : a_out["chainID"] = atom1["chainID"]
    except : a_out["chainID"] = ""
    try : a_out["resSeq"] = atom1["resSeq"]
    except : a_out["resSeq"] = "1"
    try : a_out["iCode"] = atom1["iCode"]
    except : a_out["iCode"] = ""
    a_out["element"] = "Z"
    
    a_out["x"] = (atom1["x"] + atom2["x"])/2
    a_out["y"] = (atom1["y"] + atom2["y"])/2
    a_out["z"] = (atom1["z"] + atom2["z"])/2
    
    return a_out

    
def CenterGuanidium (l_atom_connectN, l_atom_lig) : 
    
    
    l_atom_connect = []
    for atom_connectN in l_atom_connectN : 
        l_atom_connect.append (atom_connectN["element"])
    
    if l_atom_connect == ['N', 'C', 'C'] : 
        for C_atom in l_atom_connectN[1:] : 
            l_C_atom, l_connect = retrieveAtom.atomConnect(l_atom_lig, C_atom["serial"])
            if l_connect == ["C", "N", "N", "N"] : 
                return deepcopy(l_C_atom[0])
    
    
    if l_atom_connect == ["N", "C"] : 
        return deepcopy(l_atom_connectN[-1])
 
 
    print "ERROR"
    return "ERROR"
        

def CenterAcidCarboxylic (l_atom_connectO, l_atom_lig):
    
    
    
    if l_atom_connectO[1]["element"] == "C" :
         
        return deepcopy(l_atom_connectO[1])
    
    print "ERROR"
    return "ERROR"
    
    


def anglePlanImidazole(atomN1, atomC1, atomC2, atomN3):
    
    plan = equationPlan(atomC1, atomN1, atomC2)   
    
    xN = atomN1["x"]
    yN = atomN1["y"]
    zN = atomN1["z"]

    minDiff = 1000
    pointTest = {}
    pointTest["x"] = xN
    pointTest["y"] = yN
    pointTest["z"] = zN
    pointTest["element"] = "S"
    pointRef = pointTest
    alpha = angleVector(atomC1, atomN1, atomC2)

    for x in range (-5, 5) :
        pointTest["x"] = xN + x * 0.2
        for y in range (-5, 5) :
            pointTest["y"] = yN + y * 0.2
            for z in range (-5, 5) :
                pointTest["z"] = zN + z * 0.2
                Diff = checkPoint(pointTest, atomN1, atomC1, atomC2, plan, alpha)
                if Diff < minDiff :
                    pointRef = deepcopy(pointTest)
                    minDiff = Diff
                    
    # print minDiff           
    return [angleVector(atomN3, atomN1, pointRef)], pointRef
    
