from re import search, sub
from os import path

import parsing
import calcul
import pathManage
import writePDBfile
import structure




def AreaPrimary (pr_init):
    """calculation of volume around nitrogen of primary amine
    in: filePDB with only primary amine, extreme value of l_angle, structure subs
    out: file format filePDB with water"""
    
    subs = "I"
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
    for x_test in [atomN["x"] + x * 0.2 for x in range (-50,60)] : 
        for y_test in [atomN["y"] + y * 0.2 for y in range (-50,60)] : 
            for z_test in [atomN["z"] + z * 0.2 for z in range (-50,60)] :  
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
    
    
    
    
    
    
def AreaSecondary (pr_init):
    """calculation of volume around nitrogen of primary amine
    in: filePDB with only primary amine, extreme value of l_angle, structure subs
    out: file format filePDB with water"""
    
    subs = "II"
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
    for x_test in [atomN["x"] + x * 0.2 for x in range (-50,60)] : 
        for y_test in [atomN["y"] + y * 0.2 for y in range (-50,60)] : 
            for z_test in [atomN["z"] + z * 0.2 for z in range (-50,60)] :  
                atom_test = structure.genericAtom(x_test, y_test, z_test)
                dist = calcul.distanceTwoatoms(atom_test, atomN)
                if dist < d_sup and dist > d_inf: 
                    angle1 = calcul.Angle3Atoms(atomC1, atomN, atom_test)
                    angle2 = calcul.Angle3Atoms(atomC2, atomN, atom_test)
                    if angle1 >= angle_inf and angle2 >= angle_inf:
                        if angle1 <= angle_sup and angle2 <= angle_sup : 
                            serial = serial + 1
                            atom_test["serial"] = serial
                            atom_test["resSeq"] = serial
                            writePDBfile.coordinateStructure(atom_test, "HETATM", filout)
                            
    filout.close()
    
    WriteParameter (pr_volume + subs + ".param", subs, def_volume, serial)
        




def AeraTertiary (pr_init):
    """calculation of volume around nitrogen of primary amine
    in: filePDB with only primary amine, extreme value of l_angle, structure subs
    out: file format filePDB with water"""
    
    subs = "III"
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
    for x_test in [atomN["x"] + x * 0.2 for x in range (-50,60)] : 
        for y_test in [atomN["y"] + y * 0.2 for y in range (-50,60)] : 
            for z_test in [atomN["z"] + z * 0.2 for z in range (-50,60)] :  
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


  
  
  
def AreaImidazole (pr_init):
    """calculation of volume around nitrogen of primary amine
    in: filePDB with only primary amine, extreme value of l_angle, structure subs
    out: file format filePDB with water
    only one half of imidazole are inspected"""
    
    subs = "IMD"
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
        
    atom_center = calcul.CenterPoint(atomN1, atomN3)
    
    serial = 0
    for x_test in [atom_center["x"] + x * 0.2 for x in range (-50,60)] : 
        for y_test in [atom_center["y"] + y * 0.2 for y in range (-50,60)] : 
            for z_test in [atom_center["z"] + z * 0.2 for z in range (-50,60)] :  
                
                atom_test = structure.genericAtom(x_test, y_test, z_test)
                dist = calcul.distanceTwoatoms(atom_test, atom_center)
                
                # case where the calculation of angle is impossible
                try : angle = calcul.Angle3Atoms(atomN1, atom_center, atom_test)
                except : continue
                
                if dist <= d_sup and dist >= d_inf: 
                    if angle >= angle_inf and angle <= angle_sup :
                        serial = serial + 1
                        atom_test["serial"] = serial
                        atom_test["resSeq"] = serial
                        
                        writePDBfile.coordinateStructure(atom_test, "HETATM", filout)
                    
                    
    filout.close()
    WriteParameter (pr_volume + subs + ".param", subs, def_volume, serial) 




def AreaCOO (pr_init):
    """calculation of volume around nitrogen of primary amine
    in: filePDB with only primary amine, extreme value of l_angle, structure subs
    out: file format filePDB with water
    only one half of imidazole are inspected"""
    
    subs = "COO"
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
        if atom_sub["name"] == "O01" : 
            atomO1 = atom_sub
        if atom_sub["name"] == "O02" : 
            atomO2 = atom_sub
        if atom_sub["name"] == "C02" : 
            atomC2 = atom_sub
        if atom_sub["name"] == "C01" : 
            atomC1 = atom_sub
        
    atom_center = calcul.CenterPoint(atomO1, atomO2)
    
    serial = 0
    for x_test in [atom_center["x"] + x * 0.2 for x in range (-50,60)] : 
        for y_test in [atom_center["y"] + y * 0.2 for y in range (-50,60)] : 
            for z_test in [atom_center["z"] + z * 0.2 for z in range (-50,60)] :  
                
                atom_test = structure.genericAtom(x_test, y_test, z_test)
                dist = calcul.distanceTwoatoms(atom_test, atom_center)
                
                # case where the calculation of angle is impossible
                try : angle = calcul.Angle3Atoms(atomC2, atom_center, atom_test)
                except : continue
                
                if dist <= d_sup and dist >= d_inf: 
                    if angle >= angle_inf and angle <= angle_sup :
                        serial = serial + 1
                        atom_test["serial"] = serial
                        atom_test["resSeq"] = serial
                        
                        writePDBfile.coordinateStructure(atom_test, "HETATM", filout)
                    
                    
    filout.close()
    WriteParameter (pr_volume + subs + ".param", subs, def_volume, serial) 


  
def AreaGuanidium (pr_init):
    """calculation of volume around nitrogen of primary amine
    in: filePDB with only primary amine, extreme value of l_angle, structure subs
    out: file format filePDB with water"""
    
    subs = "GAI"
    pr_volume = pathManage.CreatePathDir(pr_init + "Volume/")
    filout = open (pr_volume + "volume_" + subs + ".pdb", "w")
    l_atom_sub = structure.substructureCoord(subs)
    def_volume = structure.criteraAngle(subs)
    writePDBfile.coordinateSection(filout, l_atom_sub, "HETATM")
    #filout.close ()
    #ddd

    
    angle_inf = def_volume["angle"][0]
    angle_sup = def_volume["angle"][1]
    
    
    d_inf = def_volume["distance"][0]
    d_sup = def_volume["distance"][1]
    
    
    for atom_sub in l_atom_sub : 
        if atom_sub["name"] == "N01" : 
            atomN1 = atom_sub
        if atom_sub["name"] == "C01" : 
            atomC1 = atom_sub
        
    count = 0
    for x_test in [atomC1["x"] + x * 0.3 for x in range (-100,100)] : 
        for y_test in [atomC1["y"] + y * 0.3 for y in range (-100,100)] : 
            for z_test in [atomC1["z"] + z * 0.3 for z in range (-100,100)] :  
                atom_test = structure.genericAtom(x_test, y_test, z_test)
                distance = calcul.distanceTwoatoms(atom_test, atomC1)
                l_angleC1 = calcul.anglePrimaryAmineCalculVol(atomC1, atomN1, atom_test)
                
#                 print distance, l_angleC1
                if distance < d_sup and distance > d_inf: 
                    if l_angleC1[0] > angle_inf and l_angleC1[0] < angle_sup :
                        count = count + 1
                        atom_test["count"] = count
                        atom_test["resSeq"] = count
                        writePDBfile.coordinateStructure(atom_test, "HETATM", filout)
                
                
                            
    filout.close()
    WriteParameter (pr_volume + subs + ".param", subs, def_volume, count) 
  
    

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

