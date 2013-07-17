import structure
import loadFile
import retrieveAtom
import writePDBfile
import repertory
import tool
import writeFile
import runScriptR


from numpy import *
from copy import deepcopy



def groupAtomCoord (list_atom1):
    
    l_out = []
    for atom in list_atom1 : 
        l_out.append ([float(atom["x"]), float(atom["y"]), float(atom["z"])])
    
    return l_out


# def transformMatrix2List (mat_cord):
#     
#     l_out = []
#     group_atom_rotated = array(mat_cord) # tranform array
#     for rowmat in group_atom_rotated : 
#         d_atom = {}
#         d_atom["x"] = rowmat[0]
#         d_atom["y"] = rowmat[1]
#         d_atom["z"] = rowmat[2]
#         d_atom["element"] = "C"
#         d_atom["resSeq"] = "LIZ"
#         
#         l_out.append(d_atom)
#     return l_out
        

def rigid_transform_3D(A, B):
    
    print len (A), len (B)
    
    assert len(A) == len(B)
    if len(A) != len(B) :
        print "Not same number of atoms" 
        return None

    N = A.shape[0] # total points

    centroid_A = sum(A,axis=0) / float(N)
    centroid_B = sum(B,axis=0) / float(N)


#     # centre the points
    AA = A - tile(centroid_A, (N, 1))
#     print tile(centroid_A, (N, 1)), "11111111"
    BB = B - tile(centroid_B, (N, 1))
# 
#     # dot is matrix multiplication for array
    H = transpose(AA) * BB
# 
    U, S, Vt = linalg.svd(H)
# 
    R = Vt.T * U.T
    t = -R*centroid_A.T + centroid_B.T
# 
#     # special reflection case
    if linalg.det(R) < 0:
        print "Reflection detected"
        Vt[2,:] *= -1
        R2 = Vt.T * U.T
        t2 = -R2*centroid_A.T + centroid_B.T
    
    if "R2" in locals() : 
        A1 = applyTranformation (R, t, A )
        A2= applyTranformation (R2, t2, A )
        rmse1 = rmse(A1, B)
        rmse2 = rmse(B, A2)
        if rmse2 < rmse1 : 
            return R2, t2
        
    return R, t




def applyTranformation (matrix_rotation, matrix_transloc, matrix_points = None , l_atom_in = []):
    
    
    # case with list of point in input
    if matrix_points == None : 
        matrix_points = mat(array(groupAtomCoord (l_atom_in)))
    
    # A2 new coord
    n = len(matrix_points) # number of points
    A2 = (matrix_rotation*matrix_points.T) + tile(matrix_transloc, (1, n))
    A2 = A2.T
    
    # return only matrix if list atom in input -> empty

    if l_atom_in == [] : 
        return A2
    else :
        return transformListPoint (l_atom_in, A2)



def transformListPoint (l_atom_in, matrix_coord_before_transloc) :
    
    l_out = deepcopy(l_atom_in)
    nb_atom = len (matrix_coord_before_transloc)
    
    i = 0
    while i < nb_atom : 
        l_out[i]["x"] = matrix_coord_before_transloc[i].tolist()[0][0]
        l_out[i]["y"] = matrix_coord_before_transloc[i].tolist()[0][1]
        l_out[i]["z"] = matrix_coord_before_transloc[i].tolist()[0][2]
    
        i = i + 1
    
    
    return l_out
     
    
    
    

def rmse (m_points1, m_points2):
    
    # condition same size matrix
    n = len (m_points1)
    err = m_points2-m_points1
    err = multiply(err, err)
    err = sum(err)
    rmse = sqrt(err/n)
    
    return rmse
    







def globalNeighbor (atom_interest_close, substruct, p_dir_result) : 
    
    p_dir_result = repertory.imposeNeighbors (p_dir_result)
    l_at_ref =  structure.substructureCoord(substruct)
    l_superimpose = []
    
    for at_central in atom_interest_close[substruct] : 
        PDB_ID = at_central["PDB"]
        serial_at_central = at_central["serial"]
        name_ligand =  at_central["resName"] 
        
        print PDB_ID, serial_at_central
        # all atom ligand
        l_at_lig = loadFile.ligandInPDB(PDB_ID, name_ligand)
#         for at_ligand in l_at_lig : 
#             print at_ligand
        l_at_subs = retrieveAtom.substructure (substruct, serial_at_central, l_at_lig)
        
        
        
        v_atom_ref = mat(array(groupAtomCoord(l_at_ref[0:3])))
        v_atom_central = mat(array(groupAtomCoord(l_at_subs[0:3])))


        rotation, translocation =  rigid_transform_3D(v_atom_central, v_atom_ref)
        
        v_atom_rotated = applyTranformation(rotation, translocation, v_atom_central)
        l_atom_rotated = applyTranformation(rotation, translocation, l_atom_in=l_at_subs)

        print rmse(v_atom_central, v_atom_ref), "RMSE 1"
        print rmse(v_atom_ref, v_atom_rotated), "RMSE 2"


#         print v_atom_rotated
#         print "************compare**********"
#         print l_at_subs
#         print l_atom_rotated
#         print "/////////////////////////////"
        
        l_atom_neighbors = at_central["neighbors"]
        try : 
            l_atom_neighbor_rotated = applyTranformation(rotation, translocation, l_atom_in=l_atom_neighbors)
            l_superimpose = l_superimpose + l_atom_rotated + l_atom_neighbor_rotated
        except : 
            continue
    
    # color with b factor
    tool.colorAtomType (l_superimpose)
    
    
    p_file_coord = writeFile.coordinates3D (l_superimpose, p_dir_result + substruct + "_neigbor.coord", substruct) 
    writePDBfile.coordinateSection(p_dir_result + substruct + "_neigbor.pdb", l_superimpose , "HETATM", "Superimpose neighbors " + substruct)
    runScriptR.plot3D (p_file_coord)
    
        
            
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    



