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


def coplanar (atom, listAtom):
    """Calculate the orthogonal distance between nitrogen atom of tertiary amine and plan with 3 carbons connect
    in : - atom of nitrogen -> dictionnary
         - all atom of ligand -> list of atom dictionnary
    out : - distance -> float"""

    matrix = atom["connect"]
    if len(matrix) < 4:
        print "Atom does not 3 bonds !!"
        return

    else:
        point = [float(retrieveAtom.serial(matrix[0], listAtom)["x"]), float(retrieveAtom.serial(matrix[0], listAtom)["y"]), float(retrieveAtom.serial(matrix[0], listAtom)["z"])]

        d = symbols("d")
        
        atom1 = retrieveAtom.serial(matrix[1], listAtom)
        atom2 = retrieveAtom.serial(matrix[2], listAtom)
        atom3 = retrieveAtom.serial(matrix[3], listAtom)
        
        
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




def buildConnectMatrix(listAtomLigand, namePDB):
    """Building the matrix connect. threshold -> 1.6 angstroms
    in : - Atoms of ligands -> list of atoms -> dictionnary
         - PDB file associted of ligand
    out : - Ligand with connect matrix"""


    linesPDB = loadFile.openPdbFile(namePDB)
    groupAtomPDB = []
    for line in linesPDB:
        if search("^ATOM", line) or search("^HETATM", line):
            groupAtomPDB.append(parsing.lineCoords(line))


    for atomLigand in listAtomLigand:
        atomLigand["connect"].append(atomLigand["serial"])
        for atomPDB in groupAtomPDB:
            distance = distanceTwoatoms(atomLigand, atomPDB)
            if distance < 1.7 and distance != 0:
                if not atomPDB["serial"] in atomLigand["connect"]:
                    atomLigand["connect"].append(atomPDB["serial"])

    return listAtomLigand


def equationPlan (A, B, C):
    
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
        return 0
    
    

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


def angleSecondaryAmine(atomNitrogen, atomCounterIon, listAtomLigand):
    """Calcul for nitrogen and counter ion angles
    in: atom nitrogen, atom Counter Ions
    out: list of 2 angles"""
    
    matrix = atomNitrogen["connect"]
    if len(matrix) < 3:
        print "Atom does not 2 bonds !!"
        return

    else:

        carbon1 = retrieveAtom.serial(matrix[1], listAtomLigand)
        carbon2 = retrieveAtom.serial(matrix[2], listAtomLigand)

        distanceNitrogenCounterIon = distanceTwoatoms(atomNitrogen, atomCounterIon)

        distanceCarbon1CounterIon = distanceTwoatoms(atomCounterIon, carbon1)
        distanceCarbon2CounterIon = distanceTwoatoms(atomCounterIon, carbon2)

        distanceNitrogenCarbon1 = distanceTwoatoms(atomNitrogen, carbon1)
        distanceNitrogenCarbon2 = distanceTwoatoms(atomNitrogen, carbon2)

        angle1 = degrees(acos((distanceNitrogenCarbon1 * distanceNitrogenCarbon1 + distanceNitrogenCounterIon * distanceNitrogenCounterIon - distanceCarbon1CounterIon * distanceCarbon1CounterIon) / (2 * distanceNitrogenCarbon1 * distanceNitrogenCounterIon)))
        angle2 = degrees(acos((distanceNitrogenCarbon2 * distanceNitrogenCarbon2 + distanceNitrogenCounterIon * distanceNitrogenCounterIon - distanceCarbon2CounterIon * distanceCarbon2CounterIon) / (2 * distanceNitrogenCarbon2 * distanceNitrogenCounterIon)))


        return [angle1, angle2]


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


def angleImidazolePyridine(atomNitrogen, atomFound, listAtomLigand):
    
    matrix = atomNitrogen["connect"]
    atomC1 = retrieveAtom.serial(matrix[1], listAtomLigand)
    atomC2 = retrieveAtom.serial(matrix[2], listAtomLigand)
    
    plan = equationPlan(atomC1, atomNitrogen, atomC2)    
    
    xN = atomNitrogen["x"]
    yN = atomNitrogen["y"]
    zN = atomNitrogen["z"]

    minDiff = 100
    pointTest = {}
    pointTest["x"] = xN
    pointTest["y"] = yN
    pointTest["z"] = zN
    pointTest["element"] = "O"

    alpha = angleVector(atomC1, atomNitrogen, atomC2)

    for x in range (-5, 5) :
        pointTest["x"] = xN + x * 0.2
        for y in range (-5, 5) :
            pointTest["y"] = yN + y * 0.2
            for z in range (-5, 5) :
                pointTest["z"] = zN + z * 0.2
                Diff = checkPoint(pointTest, atomNitrogen, atomC1, atomC2, plan, alpha)
                if Diff < minDiff :
                    pointRef = deepcopy(pointTest)
                    minDiff = Diff

    try : return [angleVector(atomFound, atomNitrogen, pointRef)]
    except : return []



def angle(atomNitrogen, atomFound, listAtomLigand, typeStructure):
    """calcul angle for each structure study
    in: atom nitrogen, list atoms ligand, structure study
    out: list of angles"""
    
    if typeStructure == "Primary" :
        return anglePrimaryAmine(atomNitrogen, atomFound, listAtomLigand)
    elif typeStructure == "Secondary" :
        return angleSecondaryAmine(atomNitrogen, atomFound, listAtomLigand)
    elif typeStructure == "Tertiary" :
        return angleTertiaryAmine(atomNitrogen, atomFound, listAtomLigand)
    elif typeStructure == "Imidazole" :
        return angleImidazolePyridine(atomNitrogen, atomFound, listAtomLigand)
    elif typeStructure == "Diamine" :
        return anglePrimaryAmine(atomNitrogen, atomFound, listAtomLigand)
    elif typeStructure == "Pyridine" :
        return angleSecondaryAmine(atomNitrogen, atomFound, listAtomLigand)
    elif typeStructure == "Guanidium" :
        return anglePrimaryAmine(atomNitrogen, atomFound, listAtomLigand)
    else :
        return []
    


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


def angleSecondaryAmineCalculVol(atomNitrogen, atomC1, atomC2, atomTest):
    """calcul for three atoms angleVector for secondary amine
    in: atom nitrogen, atoms carbons and atom test
    out: list angles"""

    distanceNitrogenTest = distanceTwoatoms(atomNitrogen, atomTest)

    distanceCarbon1Test = distanceTwoatoms(atomTest, atomC1)
    distanceCarbon2Test = distanceTwoatoms(atomTest, atomC2)

    distanceNitrogenCarbon1 = distanceTwoatoms(atomNitrogen, atomC1)
    distanceNitrogenCarbon2 = distanceTwoatoms(atomNitrogen, atomC2)

    angle1 = degrees(acos((distanceNitrogenCarbon1 * distanceNitrogenCarbon1 + distanceNitrogenTest * distanceNitrogenTest - distanceCarbon1Test * distanceCarbon1Test) / (2 * distanceNitrogenCarbon1 * distanceNitrogenTest)))
    angle2 = degrees(acos((distanceNitrogenCarbon2 * distanceNitrogenCarbon2 + distanceNitrogenTest * distanceNitrogenTest - distanceCarbon2Test * distanceCarbon2Test) / (2 * distanceNitrogenCarbon2 * distanceNitrogenTest)))


    return [angle1, angle2]    
    
    
    
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



def angleImidazolePyridineCalculVol(atomNitrogen, atomC1, atomC2, atomTest):
    
    plan = equationPlan(atomC1, atomNitrogen, atomC2)   
    
    xN = atomNitrogen["x"]
    yN = atomNitrogen["y"]
    zN = atomNitrogen["z"]

    minDiff = 1000
    pointTest = {}
    pointTest["x"] = xN
    pointTest["y"] = yN
    pointTest["z"] = zN
    pointTest["element"] = "S"
    pointRef = pointTest
    alpha = angleVector(atomC1, atomNitrogen, atomC2)

    for x in range (-5, 5) :
        pointTest["x"] = xN + x * 0.2
        for y in range (-5, 5) :
            pointTest["y"] = yN + y * 0.2
            for z in range (-5, 5) :
                pointTest["z"] = zN + z * 0.2
                Diff = checkPoint(pointTest, atomNitrogen, atomC1, atomC2, plan, alpha)
                if Diff < minDiff :
                    pointRef = deepcopy(pointTest)
                    minDiff = Diff
                    
    # print minDiff           
    return [angleVector(atomTest, atomNitrogen, pointRef)], pointRef
    
