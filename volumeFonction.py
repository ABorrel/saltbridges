from re import search
from re import sub
import parsing
import calcul
import writePDBfile
from os import path


def primaryAmine(path_filePDB, inferiorLimit, superiorLimit, study):
    """calculation of volume around nitrogen of primary amine
    in: filePDB with only primary amine, extreme value of angleVector, structure study
    out: file format filePDB with water"""
    
    
    nameFileOut = "volume" + study + str(inferiorLimit) + "_" + str(superiorLimit) + ".pdb"
    superiorLimit = float(superiorLimit)
    inferiorLimit = float(inferiorLimit)
    
    filin = open (path_filePDB, "r")
    linePDB = filin.readlines()
    
    for line in linePDB :
        line = sub("\n", "", line)
        if search("^HETATM", line) : 
            atom = parsing.lineCoords(line)
            if atom["element"] == "N" : 
                atomN = atom
            elif atom["element"] == "C" : 
                atomC = atom

    fileWrite = path.dirname(path_filePDB) + "/" + nameFileOut
    fileWrite = open(fileWrite, "w")
    
    
    atomN["x"] = atomN["x"] + 10
    atomN["y"] = atomN["y"] + 10
    atomN["z"] = atomN["z"] + 10
    atomC["x"] = atomC["x"] + 10
    atomC["y"] = atomC["y"] + 10
    atomC["z"] = atomC["z"] + 10
    
    writePDBfile.coordinateStructure(atomC, "HETATM", fileWrite)
    writePDBfile.coordinateStructure(atomN, "HETATM", fileWrite)

    atomTest = {}
    atomTest["element"] = "O"
    atomTest["resName"] = "HOH"
    atomTest["name"] = "O"
    serial = 2
    
    listVal = []
    i = -5
    while i <= 5 : 
        listVal.append(i)
        i = i + 0.4
    
    count = 0
    for x in listVal:
        atomTest["x"] = x + atomN["x"]
        for y in listVal : 
            atomTest["y"] = y + atomN["y"]
            for z in listVal :
                atomTest["z"] = z + atomN["z"]
                distance = calcul.distanceTwoatoms(atomTest, atomN)
                if distance < 5 and distance > 2: 
                    angleVector = calcul.anglePrimaryAmineCalculVol(atomN, atomC, atomTest)
                    if angleVector[0] > inferiorLimit :
                        if angleVector[0] < superiorLimit : 
                            serial = serial + 1
                            atomTest["serial"] = serial
                            writePDBfile.coordinateStructure(atomTest, "HETATM", fileWrite)
                            count = count + 1
                            
    fileWrite.write("END\n")
    fileWrite.close()
    print count, study
    
    
def secondaryAmine(path_filePDB, inferiorLimit, superiorLimit, study):
    """calculation of volume arround nitrogen of secondary amine
    in: filePDB with only secondary amine and extreme value of angleVector, structure study
    out: file format filePDB with water"""

    nameFileOut = "volume" + str(study) + "_" + str(inferiorLimit) + "_" + str(superiorLimit) + ".pdb"
    superiorLimit = float(superiorLimit)
    inferiorLimit = float(inferiorLimit)

    filin = open (path_filePDB, "r")
    linePDB = filin.readlines()
    
    listAtom = []
    for line in linePDB :
        line = sub("\n", "", line)
        if search("^ATOM", line) : 
            atom = parsing.lineCoords(line)
            listAtom.append(atom)

    
    count = 0        
    for atom in listAtom :         
            if atom["element"] == "N" : 
                atomN = atom
            elif atom["element"] == "C" : 
                if "atomC1" in locals() :
                    atomC2 = atom
                else : 
                    atomC1 = atom

    fileWrite = path.dirname(path_filePDB) + "/" + nameFileOut
    fileWrite = open(fileWrite, "w")
    
    
    atomN["x"] = atomN["x"] + 10
    atomN["y"] = atomN["y"] + 10
    atomN["z"] = atomN["z"] + 10
    atomC1["x"] = atomC1["x"] + 10
    atomC1["y"] = atomC1["y"] + 10
    atomC1["z"] = atomC1["z"] + 10
    atomC2["x"] = atomC2["x"] + 10
    atomC2["y"] = atomC2["y"] + 10
    atomC2["z"] = atomC2["z"] + 10
    
    writePDBfile.coordinateStructure(atomC1, "HETATM", fileWrite)
    writePDBfile.coordinateStructure(atomC2, "HETATM", fileWrite)
    writePDBfile.coordinateStructure(atomN, "HETATM", fileWrite)

    atomTest = {}
    atomTest["element"] = "O"
    atomTest["resName"] = "HOH"
    atomTest["name"] = "O"
    serial = 3
    
    listVal = []
    i = -5
    while i <= 5 : 
        listVal.append(i)
        i = i + 0.3
    
    for x in listVal:
        atomTest["x"] = x + atomN["x"]
        for y in listVal : 
            atomTest["y"] = y + atomN["y"]
            for z in listVal :
                atomTest["z"] = z + atomN["z"]
                distance = calcul.distanceTwoatoms(atomTest, atomN)
                if distance < 5 and distance > 2: 
                    
                    angles = calcul.angleSecondaryAmineCalculVol(atomN, atomC1, atomC2, atomTest)
                    if angles[0] > inferiorLimit and angles[1] > inferiorLimit:
                        if angles[0] < superiorLimit and angles[1] < superiorLimit: 
                            serial = serial + 1
                            atomTest["serial"] = serial
                            writePDBfile.coordinateStructure(atomTest, "HETATM", fileWrite)
                            count = count + 1
                            
    fileWrite.write("END\n")
    fileWrite.close()
    print count , study
    
    
def tertiaryAmine(path_filePDB, inferiorLimit, superiorLimit, study):
    """calculation of volume around nitrogen of amine
    in: filePDB with only secondary amine and extreme value of angleVector
    out: file format filePDB with water"""   
    
    nameFileOut = "volume" + str(study) + "_" + str(inferiorLimit) + "_" + str(superiorLimit) + ".pdb"
    superiorLimit = float(superiorLimit)
    inferiorLimit = float(inferiorLimit)   
    
    filin = open (path_filePDB, "r")
    linePDB = filin.readlines()
    count = 0
    for line in linePDB :
        line = sub("\n", "", line)
        if search("^ATOM", line) : 
            atom = parsing.lineCoords(line)
            if atom["element"] == "N" : 
                atomN = atom
            elif atom["element"] == "C" : 
                if "atomC1" in locals() :
                    if "atomC2" in locals():
                        atomC3 = atom
                    else : 
                        atomC2 = atom
                else : 
                    atomC1 = atom

    fileWrite = path.dirname(path_filePDB) + "/" + nameFileOut
    fileWrite = open(fileWrite, "w")
    
    atomN["x"] = atomN["x"] + 10
    atomN["y"] = atomN["y"] + 10
    atomN["z"] = atomN["z"] + 10
    atomC1["x"] = atomC1["x"] + 10
    atomC1["y"] = atomC1["y"] + 10
    atomC1["z"] = atomC1["z"] + 10
    atomC2["x"] = atomC2["x"] + 10
    atomC2["y"] = atomC2["y"] + 10
    atomC2["z"] = atomC2["z"] + 10
    atomC3["x"] = atomC3["x"] + 10
    atomC3["y"] = atomC3["y"] + 10
    atomC3["z"] = atomC3["z"] + 10
    
    writePDBfile.coordinateStructure(atomC1, "HETATM", fileWrite)
    writePDBfile.coordinateStructure(atomC2, "HETATM", fileWrite)
    writePDBfile.coordinateStructure(atomC3, "HETATM", fileWrite)
    writePDBfile.coordinateStructure(atomN, "HETATM", fileWrite)

    atomTest = {}
    atomTest["element"] = "O"
    atomTest["resName"] = "HOH"
    atomTest["name"] = "O"
    serial = 4
    
    listVal = []
    i = -5
    while i <= 5 : 
        listVal.append(i)
        i = i + 0.3
    
    for x in listVal:
        atomTest["x"] = x + atomN["x"]
        for y in listVal : 
            atomTest["y"] = y + atomN["y"]
            for z in listVal :
                atomTest["z"] = z + atomN["z"]
                distance = calcul.distanceTwoatoms(atomTest, atomN)
                if distance < 5 and distance > 2: 
                    
                    angles = calcul.angleTertiaryAmineCalculVol(atomN, atomTest, atomC1, atomC2, atomC3)
                    if angles[0] > inferiorLimit and angles[1] > inferiorLimit and angles[2] > inferiorLimit:
                        if angles[0] < superiorLimit and angles[1] < superiorLimit and angles[2] < superiorLimit:
                            serial = serial + 1
                            atomTest["serial"] = serial
                            writePDBfile.coordinateStructure(atomTest, "HETATM", fileWrite)
                            count = count + 1
                            
    fileWrite.write("END\n")
    fileWrite.close()
    print count, study

def guanidium(path_filePDB, inferiorLimit, superiorLimit, infSecondary, supSecondary, study):
    """calculation of volume around nitrogen of primary amine
    in: filePDB with only primary amine, extreme value of angleVector, structure study
    out: file format filePDB with water"""
    
    
    nameFileOut = "volume" + study + str(inferiorLimit) + "_" + str(superiorLimit) + "_" + str(infSecondary) + "_" + str(supSecondary) + ".pdb"
    superiorLimit = float(superiorLimit)
    inferiorLimit = float(inferiorLimit)
    
    filin = open (path_filePDB, "r")
    linePDB = filin.readlines()
    
    listAtom = []
    for line in linePDB :
        line = sub("\n", "", line)
        if search("^ATOM", line) or search("^HETATM", line) : 
            atom = parsing.lineCoords(line)
            listAtom.append(atom)
    
    for atom in listAtom :
        if not "atomN1" in locals() : 
            if atom["element"] == "N" : 
                atomN1 = atom
                continue
        if not "atomN2" in locals() : 
            if atom["element"] == "N" : 
                atomN2 = atom
                continue
        if not "atomN3" in locals() : 
            if atom["element"] == "N" : 
                atomN3 = atom
                continue
            
        if atom["element"] == "C" : 
            if not "atomC" in locals():
                atomC = atom
            else : 
                atomC1 = atom
                    
    fileWrite = path.dirname(path_filePDB) + "/" + nameFileOut
    fileWrite = open(fileWrite, "w")
    
    
    writePDBfile.coordinateStructure(atomC, "HETATM", fileWrite)
    writePDBfile.coordinateStructure(atomC1, "HETATM", fileWrite)
    writePDBfile.coordinateStructure(atomN1, "HETATM", fileWrite)
    writePDBfile.coordinateStructure(atomN2, "HETATM", fileWrite)
    writePDBfile.coordinateStructure(atomN3, "HETATM", fileWrite)
   
    atomTest = {}
    atomTest["element"] = "O"
    atomTest["resName"] = "HOH"
    atomTest["name"] = "O"
    serial = 2
    
    listVal = []
    i = -5
    while i <= 15 : 
        listVal.append(i)
        i = i + 0.4
    
    for x in listVal:
        atomTest["x"] = x + atomN2["x"]
        for y in listVal : 
            atomTest["y"] = y + atomN2["y"]
            for z in listVal :
                atomTest["z"] = z + atomN2["z"]
                distance = calcul.distanceTwoatoms(atomTest, atomN2)
                if distance < 5 and distance > 2: 
                    angleVector = calcul.anglePrimaryAmineCalculVol(atomN2, atomC, atomTest)
                    if angleVector[0] > inferiorLimit :
                        if angleVector[0] < superiorLimit : 
                            serial = serial + 1
                            atomTest["serial"] = serial
                            writePDBfile.coordinateStructure(atomTest, "HETATM", fileWrite)
    
    for x in listVal:
        atomTest["x"] = x + atomN3["x"]
        for y in listVal : 
            atomTest["y"] = y + atomN3["y"]
            for z in listVal :
                atomTest["z"] = z + atomN3["z"]
                distance = calcul.distanceTwoatoms(atomTest, atomN3)
                if distance < 5 and distance > 2: 
                    angleVector = calcul.anglePrimaryAmineCalculVol(atomN3, atomC, atomTest)
                    if angleVector[0] > inferiorLimit :
                        if angleVector[0] < superiorLimit : 
                            serial = serial + 1
                            atomTest["serial"] = serial
                            writePDBfile.coordinateStructure(atomTest, "HETATM", fileWrite)
                            
                            
    for x in listVal:
        atomTest["x"] = x + atomN1["x"]
        for y in listVal : 
            atomTest["y"] = y + atomN1["y"]
            for z in listVal :
                atomTest["z"] = z + atomN1["z"]
                distance = calcul.distanceTwoatoms(atomTest, atomN1)
                if distance < 5 and distance > 2: 
                    
                    angles, atomRef = calcul.angleImidazolePyridineCalculVol(atomN1, atomC1, atomC, atomTest)
                    # print angles
                    if angles[0] > infSecondary :
                        if angles[0] < supSecondary : 
#                             print "write"
                            serial = serial + 1
                            atomTest["serial"] = serial
                            writePDBfile.coordinateStructure(atomTest, "HETATM", fileWrite)
    
    
    
    fileWrite.write("END\n")
    fileWrite.close()
    
    
def pyridine(path_filePDB, inferiorLimit, superiorLimit, study):
    """calculation of volume arround nitrogen of secondary amine
    in: filePDB with only secondary amine and extreme value of angleVector, structure study
    out: file format filePDB with water"""

    nameFileOut = "volume" + str(study) + "_" + str(inferiorLimit) + "_" + str(superiorLimit) + ".pdb"
    superiorLimit = float(superiorLimit)
    inferiorLimit = float(inferiorLimit)

    filin = open (path_filePDB, "r")
    linePDB = filin.readlines()
    
    listAtom = []
    for line in linePDB :
        line = sub("\n", "", line)
        if search("^HETATM", line) : 
            atom = parsing.lineCoords(line)
            listAtom.append(atom)
            
    for atom in listAtom : 
        if atom["element"] == "N" : 
            atomN = atom
        elif atom["name"] == "C08" : 
            atomC2 = atom
        elif  atom["name"] == "C06" :   
            atomC1 = atom

    fileWrite = path.dirname(path_filePDB) + "/" + nameFileOut
    fileWrite = open(fileWrite, "w")
    
    for atom in listAtom : 
        writePDBfile.coordinateStructure(atom, "HETATM", fileWrite)

    atomTest = {}
    atomTest["element"] = "O"
    atomTest["resName"] = "HOH"
    atomTest["name"] = "O"
    serial = 3
    
    listVal = []
    i = -4
    while i <= 4 : 
        listVal.append(i)
        i = i + 0.5
    
    for x in listVal:
        atomTest["x"] = x + atomN["x"]
        for y in listVal : 
            atomTest["y"] = y + atomN["y"]
            for z in listVal :
                atomTest["z"] = z + atomN["z"]
                distance = calcul.distanceTwoatoms(atomTest, atomN)
                if distance < 5 and distance > 2: 
                    
                        angles, pointRef = calcul.angleImidazolePyridineCalculVol(atomN, atomC1, atomC2, atomTest)
                    

                        if angles != None :
                            if angles[0] > inferiorLimit :
                                if angles[0] < superiorLimit : 
                                    serial = serial + 1
                                    atomTest["serial"] = serial
                                    writePDBfile.coordinateStructure(atomTest, "HETATM", fileWrite)
                                
    # writePDBfile.coordinateStructure(pointRef, "HETATM", fileWrite) 
               
    fileWrite.write("END\n")
    fileWrite.close()
    


def diamine(path_filePDB, inferiorLimit, superiorLimit, study):
    """calculation of volume around nitrogen of primary amine
    in: filePDB with only primary amine, extreme value of angleVector, structure study
    out: file format filePDB with water"""
    
    
    nameFileOut = "volume" + study + str(inferiorLimit) + "_" + str(superiorLimit) + ".pdb"
    superiorLimit = float(superiorLimit)
    inferiorLimit = float(inferiorLimit)
    
    filin = open (path_filePDB, "r")
    linePDB = filin.readlines()
    listAtom = []
    
    for line in linePDB :
        line = sub("\n", "", line)
        if search("^HETATM", line) or search("^ATOM", line) : 
            atom = parsing.lineCoords(line)
            listAtom.append(atom)
            if atom["element"] == "N" :
                if "atomN1" in locals() : 
                    atomN2 = atom
                else : 
                    atomN1 = atom
                    
            elif atom["name"] == "C01" : 
                atomC = atom

    fileWrite = path.dirname(path_filePDB) + "/" + nameFileOut
    fileWrite = open(fileWrite, "w")
    
    for atom in listAtom : 
        writePDBfile.coordinateStructure(atom, "HETATM", fileWrite)

    atomTest = {}
    atomTest["element"] = "O"
    atomTest["resName"] = "HOH"
    atomTest["name"] = "O"
    serial = 2
    
    listVal = []
    i = -5
    while i <= 5 : 
        listVal.append(i)
        i = i + 0.4
    
    for x in listVal:
        atomTest["x"] = x + atomN1["x"]
        for y in listVal : 
            atomTest["y"] = y + atomN1["y"]
            for z in listVal :
                atomTest["z"] = z + atomN1["z"]
                distance = calcul.distanceTwoatoms(atomTest, atomN1)
                if distance < 5 and distance > 2: 
                    angleVector = calcul.anglePrimaryAmineCalculVol(atomN1, atomC, atomTest)
                    if angleVector[0] > inferiorLimit :
                        if angleVector[0] < superiorLimit : 
                            serial = serial + 1
                            atomTest["serial"] = serial
                            writePDBfile.coordinateStructure(atomTest, "HETATM", fileWrite)
                            
    for x in listVal:
        atomTest["x"] = x + atomN2["x"]
        for y in listVal : 
            atomTest["y"] = y + atomN2["y"]
            for z in listVal :
                atomTest["z"] = z + atomN2["z"]
                distance = calcul.distanceTwoatoms(atomTest, atomN2)
                if distance < 5 and distance > 2: 
                    angleVector = calcul.anglePrimaryAmineCalculVol(atomN2, atomC, atomTest)
                    if angleVector[0] > inferiorLimit :
                        if angleVector[0] < superiorLimit : 
                            serial = serial + 1
                            atomTest["serial"] = serial
                            writePDBfile.coordinateStructure(atomTest, "HETATM", fileWrite)
                            
    fileWrite.write("END\n")
    fileWrite.close()


def imidazole(path_filePDB, inferiorLimit, superiorLimit, study):
    """calculation of volume arround nitrogen of secondary amine
    in: filePDB with only secondary amine and extreme value of angleVector, structure study
    out: file format filePDB with water"""

    nameFileOut = "volume" + str(study) + "_" + str(inferiorLimit) + "_" + str(superiorLimit) + ".pdb"
    superiorLimit = float(superiorLimit)
    inferiorLimit = float(inferiorLimit)

    filin = open (path_filePDB, "r")
    linePDB = filin.readlines()
    
    listAtom = []
    for line in linePDB :
        line = sub("\n", "", line)
        if search("^HETATM", line) : 
            atom = parsing.lineCoords(line)
            listAtom.append(atom)
            
    for atom in listAtom : 
        if atom["name"] == "N1" : 
            atomN1 = atom
        if atom["name"] == "N3" : 
            atomN2 = atom
        elif atom["name"] == "C2" : 
            atomC2 = atom
        elif  atom["name"] == "C4" :   
            atomC4 = atom
        elif  atom["name"] == "C5" :   
            atomC5 = atom    

    fileWrite = path.dirname(path_filePDB) + "/" + nameFileOut
    fileWrite = open(fileWrite, "w")
    
    for atom in listAtom : 
        writePDBfile.coordinateStructure(atom, "HETATM", fileWrite)

    atomTest = {}
    atomTest["element"] = "O"
    atomTest["resName"] = "HOH"
    atomTest["name"] = "O"
    serial = 3
    
    listVal = []
    i = -4
    while i <= 4 : 
        listVal.append(i)
        i = i + 0.5
    
    for x in listVal:
        atomTest["x"] = x + atomN1["x"]
        for y in listVal : 
            atomTest["y"] = y + atomN1["y"]
            for z in listVal :
                atomTest["z"] = z + atomN1["z"]
                distance = calcul.distanceTwoatoms(atomTest, atomN1)
                if distance < 5 and distance > 2: 
                    
                        angles, pointRef = calcul.angleImidazolePyridineCalculVol(atomN1, atomC2, atomC5, atomTest)
                    

                        if angles != None :
                            if angles[0] > inferiorLimit :
                                if angles[0] < superiorLimit : 
                                    serial = serial + 1
                                    atomTest["serial"] = serial
                                    writePDBfile.coordinateStructure(atomTest, "HETATM", fileWrite)
    
    
    
    for x in listVal:
        atomTest["x"] = x + atomN2["x"]
        for y in listVal : 
            atomTest["y"] = y + atomN2["y"]
            for z in listVal :
                atomTest["z"] = z + atomN2["z"]
                distance = calcul.distanceTwoatoms(atomTest, atomN2)
                if distance < 5 and distance > 2: 
                    
                        angles, pointRef = calcul.angleImidazolePyridineCalculVol(atomN2, atomC2, atomC4, atomTest)
                    

                        if angles != None :
                            if angles[0] > inferiorLimit :
                                if angles[0] < superiorLimit : 
                                    serial = serial + 1
                                    atomTest["serial"] = serial
                                    writePDBfile.coordinateStructure(atomTest, "HETATM", fileWrite) 
    
                                
               
    fileWrite.write("END\n")
    fileWrite.close()
    

