from copy import deepcopy


# remove in neighbor structure compound with iron close
def removeNeighborIron (st_atom, p_summary):
    
    filout = open (p_summary,"w")
    
    ld_out = {}
    for sub_struct in st_atom.keys() : 
        nb_central_atom = len (st_atom[sub_struct] )
        if nb_central_atom == 0 : 
            continue
        i = 0
        while i < nb_central_atom : 
            l_neighbors = st_atom[sub_struct][i]["neighbors"]
            ion_close = closeIron (l_neighbors)
            if ion_close != 0 :
                if not sub_struct in ld_out.keys () : 
                    ld_out [sub_struct] = []
                ld_out [sub_struct].append (deepcopy(st_atom[sub_struct][i]))
                filout.write (str(st_atom[sub_struct][i]["PDB"]) + "\t" + str (sub_struct) + "\t" + str(ion_close["resName"]) + "\t" + str (ion_close["distance"]) + "\n")
                del st_atom[sub_struct][i]
                nb_central_atom = nb_central_atom - 1
                continue
            else : 
                i = i + 1
    
    filout.close ()              
    return ld_out
                    
    
    
def closeIron (l_neighbors) : 
    # peut etre changer pour linstant je suprime si pas dans la liste des aa
    l_ions = ["ZN", "FE2", "ZN", "CU", "SO4", "CA","NI", "LU", "TB", "NA"]
    listAminoAcid = ["ILE", "LEU", "LYS", "PHE", "TYR", "VAL", "SER", "MET", "ARG", "TRP", "PRO", "GLY", "GLU", "ASN", "HIS", "ALA", "ASP", "GLN", "THR", "CYS", "HOH"]
    nb_neigbor = len (l_neighbors)
    if nb_neigbor == 0 : 
        return 0
    i = 0
    while i < nb_neigbor : 
        if not l_neighbors[i]["resName"] in listAminoAcid :
            return deepcopy(l_neighbors[i])
        else:
            i = i + 1
    
    return 0



