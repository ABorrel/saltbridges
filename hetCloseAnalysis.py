from copy import deepcopy


# remove in neighbor structure compound with iron close
def removeNeighborIron (st_atom, p_summary):
    
    l_ions = ['CR', 'AL', 'MG', 'MN', 'ZN', 'CA', 'FE', 'CU', 'CD', 'NI', 'PD', 'CO', 'BA', 'HG', 'BE', 'LI', 'NA','CS', 'K']
    l_DNA = ["A", "T", "C", "G"]
    
    filout = open (p_summary,"w")
    d_count = {}
    
    ld_out = {}
    for subs in st_atom.keys() : 
        if not subs in d_count.keys ():
            d_count[subs] = {}
            d_count[subs]["Ions"] = 0
            d_count[subs]["DNA"] = 0
            d_count[subs]["other"] = 0
            
        nb_central_atom = len (st_atom[subs] )
        if nb_central_atom == 0 : 
            continue
        i = 0
        while i < nb_central_atom : 
            l_neighbors = st_atom[subs][i]["neighbors"]
            ion_close = closeIon (l_neighbors)
           
            if ion_close != 0 :
                if ion_close["resName"] in l_ions : 
                    d_count[subs]["Ions"] = d_count[subs]["Ions"] + 1
                elif ion_close["resName"] in l_DNA : 
                    d_count[subs]["DNA"] = d_count[subs]["DNA"] + 1
                else : 
                    d_count[subs]["other"] = d_count[subs]["other"] + 1
                
                if not subs in ld_out.keys () : 
                    ld_out [subs] = []
                ld_out [subs].append (deepcopy(st_atom[subs][i]))
                filout.write (str(st_atom[subs][i]["PDB"]) + "\t" + str (subs) + "\t" + str(ion_close["resName"]) + "\t" + str (ion_close["distance"]) + "\n")
                del st_atom[subs][i]
                nb_central_atom = nb_central_atom - 1
                continue
            else : 
                i = i + 1
    
    
    
    filout.close ()   
    
    filout_resume = open (p_summary[0:-4] + ".res","w")
    
    for sub in d_count.keys () : 
        filout_resume.write ("===== " + str (sub) + " =====\n")
        filout_resume.write ("DNA close: " + str (d_count[sub]["DNA"]) + "\n")
        filout_resume.write ("Ions close: " + str (d_count[sub]["Ions"]) + "\n")
        filout_resume.write ("Other ligands: " + str (d_count[sub]["other"]) + "\n\n")
    filout_resume.close ()
               
    return ld_out
                    
    
    
def closeIon (l_neighbors) : 
    # peut etre changer pour linstant je suprime si pas dans la liste des aa

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



