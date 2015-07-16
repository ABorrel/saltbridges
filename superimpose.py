import structure
import loadFile
import retrieveAtom
import writePDBfile
import pathManage
import tool
import writeFile
import runScriptR
import statistic


from numpy import *
from copy import deepcopy



def groupAtomCoord (l_atom):
    
    print "--- l-19", l_atom
    
    l_out = []
    for atom in l_atom :
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
    
    
    if len(A) != len(B) :
        print "Not same number of atoms" 
        return None, None
    assert len(A) == len(B)

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
    







def globalNeighbor (atom_interest_close, subs, p_dir_result) : 
    
    p_dir_result = pathManage.imposeNeighbors (p_dir_result)
    # extract from ideal position
    l_at_ref =  structure.substructureCoord(subs)
    l_superimpose_neighbor = []
    l_superimpose_subs = []
    
    l_RMSE = []
    
    for at_central in atom_interest_close[subs] : 
        PDB_ID = at_central["PDB"]
        serial_at_central = at_central["serial"]
        name_ligand =  at_central["resName"] 
        
        print PDB_ID, serial_at_central
        # all atom ligand
        l_at_lig = loadFile.ligandInPDB(PDB_ID, name_ligand)
#         for at_ligand in l_at_lig : 
#             print at_ligand
        l_at_subs = retrieveAtom.substructure (subs, at_central, l_at_lig)
        
        
        v_atom_ref = mat(array(groupAtomCoord(l_at_ref)))
        
        v_atom_central = mat(array(groupAtomCoord(l_at_subs)))

        rotation, translocation =  rigid_transform_3D(v_atom_central, v_atom_ref)
        if rotation == None or translocation == None : 
            continue
        
        v_atom_rotated = applyTranformation(rotation, translocation, v_atom_central)
        l_subs_rotated = applyTranformation(rotation, translocation, l_atom_in = l_at_subs)

        RMSE_rot = rmse(v_atom_ref, v_atom_rotated)
        print rmse(v_atom_central, v_atom_ref), "RMSE 1"
        print rmse(v_atom_ref, v_atom_rotated), "RMSE 2"

        l_RMSE.append (str(RMSE_rot))

#         print v_atom_rotated
#         print "************compare**********"
#         print l_at_subs
#         print l_subs_rotated
#         print "/////////////////////////////"
        
        l_atom_neighbors = at_central["neighbors"]
        try : 
            l_atom_neighbor_rotated = applyTranformation(rotation, translocation, l_atom_in=l_atom_neighbors)
            l_superimpose_neighbor = l_superimpose_neighbor + l_atom_neighbor_rotated
            l_superimpose_subs = l_superimpose_subs + l_subs_rotated
        except : 
            continue
    
    # color with b factor
    tool.colorAtomType (l_superimpose_neighbor)
    
    # write gif
    pr_init_gif = p_dir_result + "/gif/" + subs + "/"
    pathManage.CreatePathDir(pr_init_gif)
    p_file_coord = writeFile.coordinates3D (l_superimpose_neighbor + l_superimpose_subs, pr_init_gif + subs + "_neigbor.coord", subs) 
    runScriptR.plot3D (p_file_coord)
    
    # write one PDB by atom close type 
    pr_init_PDB = p_dir_result + "/PDB/" + subs + "/" 
    pathManage.CreatePathDir(pr_init_PDB)
    file_RMSE = open (pr_init_PDB + "RMSE", "w")
    file_RMSE.write ("\n".join(l_RMSE) + "\n")
    file_RMSE.close ()
    writeFile.coordinates3DPDBbyNeighborType (l_superimpose_neighbor, l_superimpose_subs, subs, pr_init_PDB)
    
    
    
        
            
        
    
    
def SuperimposeFirstNeighbors (st_atom, pr_result):
    
    
    d_nb_neighbor = structure.nbNeighbor ()
    
    
    for subs in st_atom.keys () : 
        if subs == "global" : 
            continue
        
        l_at_ref =  structure.substructureCoord(subs)
        d_atom_superiposed = {}
        
        pr_superimpose = pr_result +  subs + "/"
        pathManage.CreatePathDir(pr_superimpose)
        
        nb_ind = d_nb_neighbor[subs] # number of considered neighbors
        
        for at_central in st_atom[subs] : 
            PDB_ID = at_central["PDB"]
            name_ligand =  at_central["resName"] 
        
            # all atom ligand
            l_at_lig = loadFile.ligandInPDB(PDB_ID, name_ligand)
            l_at_subs = retrieveAtom.substructure (subs, at_central, l_at_lig)
        
        
        
            v_atom_ref = mat(array(groupAtomCoord(l_at_ref[0:3])))
            v_atom_central = mat(array(groupAtomCoord(l_at_subs[0:3])))


            rotation, translocation =  rigid_transform_3D(v_atom_central, v_atom_ref)
            if rotation == None or translocation == None : 
                continue
        
            v_atom_rotated = applyTranformation(rotation, translocation, v_atom_central)
            l_atom_rotated = applyTranformation(rotation, translocation, l_atom_in=l_at_subs)

            print rmse(v_atom_central, v_atom_ref), "RMSE 1"
            print rmse(v_atom_ref, v_atom_rotated), "RMSE 2"


            #         print v_atom_rotated
            #         print "************compare**********"
            #         print l_at_subs
            #         print l_atom_rotated
        
        
        
            
            l_neighbor = deepcopy(at_central["neighbors"])
            if len (l_neighbor) < nb_ind : 
                continue
            else : 
                l_combination = []
                l_atom_neighbors = []
                for i_neighbors in range (1, nb_ind + 1) : 
                    atom_class, d_atom = statistic.searchMoreClose (l_neighbor) 
                    l_combination.append (atom_class)
                    l_atom_neighbors.append (d_atom)
                    
                
                l_combination.sort ()
                k = "_".join (l_combination)
                
                
                if not k in d_atom_superiposed.keys () : 
                    d_atom_superiposed[k] = []
                
                l_atom_neighbor_rotated = applyTranformation(rotation, translocation, l_atom_in = l_atom_neighbors)
                d_atom_superiposed[k] = d_atom_superiposed[k] + l_atom_rotated + l_atom_neighbor_rotated
                
                
        for k in d_atom_superiposed.keys () : 
            # write gif
            pr_init_gif = pr_superimpose + "gif/" 
            pathManage.CreatePathDir(pr_init_gif)
            p_file_coord = writeFile.coordinates3D (d_atom_superiposed[k], pr_init_gif + k + ".coord", subs) 
            runScriptR.plot3D (p_file_coord, option = "global")
    
            # write one PDB by atom close type 
            pr_init_PDB = pr_superimpose + "/PDB/"
            pathManage.CreatePathDir(pr_init_PDB)
            writeFile.coordinates3DPDB (d_atom_superiposed[k], subs, pr_init_PDB + k + ".pdb" )   
                
                
