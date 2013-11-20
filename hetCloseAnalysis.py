from copy import deepcopy


# remove in neighbor structure compound with iron close
def removeNeighborIron (atom_interest_close):
    
    
    if type (atom_interest_close) is list : 
        ld_out = []
        nb_central_atom = len (atom_interest_close )
        i = 0
        while i < nb_central_atom : 
            l_neighbors = atom_interest_close[i]["neighbors"]
            if closeIron (l_neighbors) == 1 :
                ld_out.append (deepcopy(atom_interest_close[i])) 
                del atom_interest_close[i]
                nb_central_atom = nb_central_atom - 1
                continue
            else : 
                i = i + 1 
                
    else :
        ld_out = {}
        for sub_struct in atom_interest_close.keys() : 
            nb_central_atom = len (atom_interest_close[sub_struct] )
            i = 0
            while i < nb_central_atom : 
                l_neighbors = atom_interest_close[sub_struct][i]["neighbors"]
                if closeIron (l_neighbors) == 1 :
                    if not sub_struct in ld_out.keys () : 
                        ld_out [sub_struct] = []
                    ld_out [sub_struct].append (deepcopy(atom_interest_close[sub_struct][i]))
                    del atom_interest_close[sub_struct][i]
                    nb_central_atom = nb_central_atom - 1
                    continue
                else : 
                    i = i + 1
                    
    return ld_out
                    
    
    
def closeIron (l_neighbors) : 
    # peut etre changer pour linstant je suprime si pas dans la liste des aa
    l_ions = ["ZN", "FE2", "ZN", "CU", "SO4", "CA","NI", "LU", "TB", "NA"]
    listAminoAcid = ["ILE", "LEU", "LYS", "PHE", "TYR", "VAL", "SER", "MET", "ARG", "TRP", "PRO", "GLY", "GLU", "ASN", "HIS", "ALA", "ASP", "GLN", "THR", "CYS", "HOH"]
    d = 10
    nb_neigbor = len (l_neighbors)
    if nb_neigbor == 0 : 
        return 0
    i = 0
    while i < nb_neigbor : 
        if l_neighbors[i]["distance"] < d : 
            res_out = l_neighbors[i]["resName"]
            res_seq = l_neighbors[i]["resSeq"]
            d = l_neighbors[i]["distance"]
        i = i + 1 
    
    if not res_out in listAminoAcid :
        print res_out
        
        
        return 1
    else : 
        return 0



