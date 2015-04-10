from re import search, sub
from os import path

import parsing
import calcul
import pathManage
import writePDBfile
import structure




def defVolumePrimary (pr_init, subs):
    """calculation of volume around nitrogen of primary amine
    in: filePDB with only primary amine, extreme value of l_angle, structure subs
    out: file format filePDB with water"""
    
    pr_volume = pathManage.CreatePathDir(pr_init + "Volume/")
    filout = open (pr_volume + "volume_" + subs + ".pdb", "w")
    l_atom_sub = structure.substructureCoord(subs)
    def_volume = structure.criteraAngle(subs)
    writePDBfile.coordinateSection(filout, l_atom_sub, "HETATM")

    
    angle_inf = def_volume["angle"][0]
    angle_sup = def_volume["angle"][1]
    
    d_inf = def_volume["distance"][0]
    d_sup = def_volume["distance"][1]
    
    
    for atom_sub in l_atom_sub :
        if atom_sub["element"] == "N" : 
            atomN = atom_sub
        elif atom_sub["element"] == "C" : 
            atomC = atom_sub

    serial = 0
    for x_test in [atomN["x"] + x * 0.1 for x in range (-50,60)] : 
        for y_test in [atomN["y"] + y * 0.1 for y in range (-50,60)] : 
            for z_test in [atomN["z"] + z * 0.1 for z in range (-50,60)] :  
                atom_test = structure.genericAtom(x_test, y_test, z_test)
                distance = calcul.distanceTwoatoms(atom_test, atomN)
                if distance < d_sup and distance > d_inf: 
                    l_angle = calcul.anglePrimaryAmineCalculVol(atomN, atomC, atom_test)
                    if l_angle[0] > angle_inf :
                        if l_angle[0] < angle_sup : 
                            serial = serial + 1
                            atom_test["serial"] = serial
                            atom_test["resSeq"] = serial
                            writePDBfile.coordinateStructure(atom_test, "HETATM", filout)
                            
    filout.close()
    
    WriteParameter (pr_volume + subs + ".param", subs, def_volume, serial)
    
    
    
    
    
    
def defVolumeSecondary (pr_init, subs):
    """calculation of volume around nitrogen of primary amine
    in: filePDB with only primary amine, extreme value of l_angle, structure subs
    out: file format filePDB with water"""
    
    pr_volume = pathManage.CreatePathDir(pr_init + "Volume/")
    filout = open (pr_volume + "volume_" + subs + ".pdb", "w")
    l_atom_sub = structure.substructureCoord(subs)
    def_volume = structure.criteraAngle(subs)
    writePDBfile.coordinateSection(filout, l_atom_sub, "HETATM")

    
    angle_inf = def_volume["angle"][0]
    angle_sup = def_volume["angle"][1]
    
    d_inf = def_volume["distance"][0]
    d_sup = def_volume["distance"][1]
    
    
    for atom_sub in l_atom_sub :
        if atom_sub["element"] == "N" : 
            atomN = atom_sub
        elif atom_sub["element"] == "C" : 
            if "atomC1" in locals() :
                atomC2 = atom_sub
            else : 
                atomC1 = atom_sub

    serial = 0
    for x_test in [atomN["x"] + x * 0.1 for x in range (-50,60)] : 
        for y_test in [atomN["y"] + y * 0.1 for y in range (-50,60)] : 
            for z_test in [atomN["z"] + z * 0.1 for z in range (-50,60)] :  
                atom_test = structure.genericAtom(x_test, y_test, z_test)
                distance = calcul.distanceTwoatoms(atom_test, atomN)
                if distance < d_sup and distance > d_inf: 
                    l_angle = calcul.angleSecondaryAmineCalculVol(atomN, atomC1, atomC2, atom_test)
                    if l_angle[0] > angle_inf and l_angle[1] > angle_inf:
                        if l_angle[0] < angle_sup and l_angle[1] < angle_sup: 
                            serial = serial + 1
                            atom_test["serial"] = serial
                            atom_test["resSeq"] = serial
                            writePDBfile.coordinateStructure(atom_test, "HETATM", filout)
                            
    filout.close()
    
    WriteParameter (pr_volume + subs + ".param", subs, def_volume, serial)
        




def defVolumeTertiary (pr_init, subs):
    """calculation of volume around nitrogen of primary amine
    in: filePDB with only primary amine, extreme value of l_angle, structure subs
    out: file format filePDB with water"""
    
    pr_volume = pathManage.CreatePathDir(pr_init + "Volume/")
    filout = open (pr_volume + "volume_" + subs + ".pdb", "w")
    l_atom_sub = structure.substructureCoord(subs)
    def_volume = structure.criteraAngle(subs)
    writePDBfile.coordinateSection(filout, l_atom_sub, "HETATM")

    
    angle_inf = def_volume["angle"][0]
    angle_sup = def_volume["angle"][1]
    
    d_inf = def_volume["distance"][0]
    d_sup = def_volume["distance"][1]
    
    
    for atom_sub in l_atom_sub :
        if atom_sub["element"] == "N" : 
            atomN = atom_sub
        elif atom_sub["element"] == "C" : 
            if "atomC1" in locals() :
                if "atomC2" in locals():
                    atomC3 = atom_sub
                else : 
                    atomC2 = atom_sub
            else : 
                atomC1 = atom_sub
        

    serial = 0
    for x_test in [atomN["x"] + x * 0.1 for x in range (-50,60)] : 
        for y_test in [atomN["y"] + y * 0.1 for y in range (-50,60)] : 
            for z_test in [atomN["z"] + z * 0.1 for z in range (-50,60)] :  
                atom_test = structure.genericAtom(x_test, y_test, z_test)
                distance = calcul.distanceTwoatoms(atom_test, atomN)
                if distance < d_sup and distance > d_inf: 
                    
                    l_angles = calcul.angleTertiaryAmineCalculVol(atomN, atom_test, atomC1, atomC2, atomC3)
                    
                    if l_angles[0] > angle_inf and l_angles[1] > angle_inf and l_angles[2] > angle_inf:
                        if l_angles[0] < angle_sup and l_angles[1] < angle_sup and l_angles[2] < angle_sup:
                            serial = serial + 1
                            atom_test["serial"] = serial
                            atom_test["resSeq"] = serial
                            writePDBfile.coordinateStructure(atom_test, "HETATM", filout)
                    
    filout.close()
    WriteParameter (pr_volume + subs + ".param", subs, def_volume, serial)


  
  
  
def defVolumeImidazole (pr_init, subs):
    """calculation of volume around nitrogen of primary amine
    in: filePDB with only primary amine, extreme value of l_angle, structure subs
    out: file format filePDB with water"""
    
    pr_volume = pathManage.CreatePathDir(pr_init + "Volume/")
    filout = open (pr_volume + "volume_" + subs + ".pdb", "w")
    l_atom_sub = structure.substructureCoord(subs)
    def_volume = structure.criteraAngle(subs)
    writePDBfile.coordinateSection(filout, l_atom_sub, "HETATM")

    
    angle_inf = def_volume["angle"][0]
    angle_sup = def_volume["angle"][1]
    
    d_inf = def_volume["distance"][0]
    d_sup = def_volume["distance"][1]
    
    
    for atom_sub in l_atom_sub : 
        if atom_sub["name"] == "N1" : 
            atomN1 = atom_sub
        if atom_sub["name"] == "N3" : 
            atomN3 = atom_sub
        
    atom_center = {}
    atom_center["x"] = (atomN1["x"] + atomN3["x"]) / 2
    atom_center["y"] = (atomN1["y"] + atomN3["y"]) / 2
    atom_center["z"] = (atomN1["z"] + atomN3["z"]) / 2
        
    
    serial = 0
    for x_test in [atom_center["x"] + x * 0.1 for x in range (-100,100)] : 
        for y_test in [atom_center["y"] + y * 0.1 for y in range (-100,100)] : 
            for z_test in [atom_center["z"] + z * 0.1 for z in range (-100,100)] :  
                atom_test = structure.genericAtom(x_test, y_test, z_test)
                distance1 = calcul.distanceTwoatoms(atom_test, atomN1)
                distance2 = calcul.distanceTwoatoms(atom_test, atomN3)
                l_angleN1 = calcul.angleImidazoleCalculVol(atomN1, atomN3, atom_test)
                l_angleN3 = calcul.angleImidazoleCalculVol(atomN3, atomN1, atom_test)
                
                if distance1 < d_sup and distance1 > d_inf: 
                    if l_angleN1[0] > angle_inf and l_angleN1[0] < angle_sup :
                        serial = serial + 1
                        atom_test["serial"] = serial
                        atom_test["resSeq"] = serial
                        writePDBfile.coordinateStructure(atom_test, "HETATM", filout)
                        continue
                            
                if distance2 < d_sup and distance2 > d_inf: 
                    if l_angleN3[0] > angle_inf and l_angleN3[0] < angle_sup :
                        serial = serial + 1
                        atom_test["serial"] = serial
                        atom_test["resSeq"] = serial
                        writePDBfile.coordinateStructure(atom_test, "HETATM", filout)
                    
                    
    filout.close()
    WriteParameter (pr_volume + subs + ".param", subs, def_volume, serial) 
  


  
def defVolumeGuanidium (pr_init, subs):
    """calculation of volume around nitrogen of primary amine
    in: filePDB with only primary amine, extreme value of l_angle, structure subs
    out: file format filePDB with water"""
    
    pr_volume = pathManage.CreatePathDir(pr_init + "Volume/")
    filout = open (pr_volume + "volume_" + subs + ".pdb", "w")
    l_atom_sub = structure.substructureCoord(subs)
    def_volume = structure.criteraAngle(subs)
    writePDBfile.coordinateSection(filout, l_atom_sub, "HETATM")

    
    angle_inf = def_volume["angle"][0]
    angle_sup = def_volume["angle"][1]
    
    d_inf = def_volume["distance"][0]
    d_sup = def_volume["distance"][1]
    
    
    for atom_sub in l_atom_sub : 
        if atom_sub["name"] == "N01" : 
            atomN1 = atom_sub
        if atom_sub["name"] == "N02" : 
            atomN2 = atom_sub
        
    atom_center = {}
    atom_center["x"] = (atomN1["x"] + atomN2["x"]) / 2
    atom_center["y"] = (atomN1["y"] + atomN2["y"]) / 2
    atom_center["z"] = (atomN1["z"] + atomN2["z"]) / 2
        
    
    serial = 0
    for x_test in [atom_center["x"] + x * 0.1 for x in range (-100,100)] : 
        for y_test in [atom_center["y"] + y * 0.1 for y in range (-100,100)] : 
            for z_test in [atom_center["z"] + z * 0.1 for z in range (-100,100)] :  
                atom_test = structure.genericAtom(x_test, y_test, z_test)
                distance1 = calcul.distanceTwoatoms(atom_test, atomN1)
                distance2 = calcul.distanceTwoatoms(atom_test, atomN3)
                l_angleN1 = calcul.angleImidazoleCalculVol(atomN1, atomN3, atom_test)
                l_angleN3 = calcul.angleImidazoleCalculVol(atomN3, atomN1, atom_test)
                
                if distance1 < d_sup and distance1 > d_inf: 
                    if l_angleN1[0] > angle_inf and l_angleN1[0] < angle_sup :
                        serial = serial + 1
                        atom_test["serial"] = serial
                        atom_test["resSeq"] = serial
                        writePDBfile.coordinateStructure(atom_test, "HETATM", filout)
                        continue
                            
                if distance2 < d_sup and distance2 > d_inf: 
                    if l_angleN3[0] > angle_inf and l_angleN3[0] < angle_sup :
                        serial = serial + 1
                        atom_test["serial"] = serial
                        atom_test["resSeq"] = serial
                        writePDBfile.coordinateStructure(atom_test, "HETATM", filout)
                    
                    
    filout.close()
    WriteParameter (pr_volume + subs + ".param", subs, def_volume, serial) 
  
    


    
    
    

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
                    
                    angles, atomRef = calcul.anglePlanImidazoleCalculVol(atomN1, atomC1, atomC, atomTest)
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
                    
                        angles, pointRef = calcul.anglePlanImidazoleCalculVol(atomN, atomC1, atomC2, atomTest)
                    

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




def WriteParameter (p_file_parameter, subs, def_volume, nb_atom):
    
        
    # file of parameter
    file_parameter = open (p_file_parameter, "w")
    file_parameter.write ("===" + subs + "===\n")
    file_parameter.write ("-> ANGLES\n")
    file_parameter.write ("Inferior: " + str (def_volume["angle"][0]) + "\n")
    file_parameter.write ("Superior: " + str (def_volume["angle"][1]) + "\n")
    file_parameter.write ("-> DISTANCE\n")
    file_parameter.write ("Inferior: " + str (def_volume["distance"][0]) + "\n")
    file_parameter.write ("Superior: " + str (def_volume["distance"][1]) + "\n")
    file_parameter.write ("==> NB atom positionned:" + str (nb_atom) + "\n")
    file_parameter.close ()     

