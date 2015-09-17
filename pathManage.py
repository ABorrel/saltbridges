from os import makedirs, listdir, path
from re import search

# globals()["repInit"] = "/home/student10/stage/"
# globals()["repInit"] = "/home/alexandre/intership/"
globals()["repInit"] = "/home/borrel/saltBridgesProject/"

repInit = "/home/borrel/saltBridgesProject/"
repPDB = "/home/borrel/PDB/"


def result(rep_add = ""):

    if rep_add != "" : 
        rep = repInit + "result/" + rep_add + "/"
    else : 
        rep = repInit + "result/"
         
    try: makedirs(rep, mode=0777)
    except: pass
    return rep

def withoutAtLeastOneSummary(dir_in):

    rep = dir_in + "withoutAtLeastOne/"
    try: makedirs(rep, mode=0777)
    except: pass
    return rep


def pathDitrectoryPDB ():
    
    return repPDB


def openPdbFile ():

    rep = repPDB
    return rep


def logFile () :

    rep = repInit + "log/"
    try: makedirs(rep, mode=0777)
    except: pass
    return rep

def aminoAcid(structure, dir_in):

    rep = dir_in + "/" + structure + "/aminoAcid/"
    try: makedirs(rep, mode=0777)
    except: pass
    return rep

def resultDistance(rep_int):

    rep = rep_int + "distanceType/"
    try: makedirs(rep, mode=0777)
    except: pass
    return rep


def resultInteraction (pr_int, name_in = ""):
    
    if name_in == "" : 
        pr_salt_bridges = pr_int + "CountInteraction/"
    else :
        pr_salt_bridges = pr_int + "CountInteraction/" + name_in + "/"
        
    if not path.isdir(pr_salt_bridges) : 
        makedirs(pr_salt_bridges, mode=0777)
    else : 
        pass
    return pr_salt_bridges



def lengthBond (rep_int):

    rep = rep_int + "lengthBond/"
    try: makedirs(rep, mode=0777)
    except: pass
    return rep


def resultSub(structure, dir_out):

    rep = dir_out + structure + "/"
    try: makedirs(rep, mode=0777)
    except: pass
    return rep



def scriptR():
    
    rep = "./"
    return rep

def resultAminoAcidGlobal():

    rep = repInit + "result/AminoAcidGlobal/"
    try: makedirs(rep, mode=0777)
    except: pass
    return rep

def globalProportionAtom(dir_in) :

    rep = dir_in + "globalProportionAtom/"
    try: makedirs(rep, mode=0777)
    except: pass
    return rep

def globalProportionType(dir_in) :

    rep = dir_in + "globalProportionType/"
    try: makedirs(rep, mode=0777)
    except: pass
    return rep


def resultAngle(pr_int, name_in = ""):

    if name_in == "" : 
        pr_angle = pr_int + "Angle/"
    else :
        pr_angle = pr_int + "Angle/" + name_in + "/"
        
    if not path.isdir(pr_angle) : 
        makedirs(pr_angle, mode=0777)
    else : 
        pass
    return pr_angle


def countGlobalProx (pr_int, name_in = ""):

    if name_in == "" : 
        pr_count = pr_int + "countProx/"
    else :
        pr_count = pr_int + "countProx/" + name_in + "/"
        
    if not path.isdir(pr_count) : 
        makedirs(pr_count, mode=0777)
    else : 
        pass
    return pr_count


def countNeighbor(pr_int, name_in = ""):

    if name_in == "" : 
        pr_count_neighbor = pr_int + "neighbor/"
    else :
        pr_count_neighbor = pr_int + "neighbor/" + name_in + "/"
        
    if not path.isdir(pr_count_neighbor) : 
        makedirs(pr_count_neighbor, mode=0777)
    else : 
        pass
    return pr_count_neighbor



def parsingDataset(pr_int):
    
    rep = pr_int + "parsingDataset/"
    try: makedirs(rep, mode=0777)
    except: pass
    return rep



def retriveDataSetFile (directory_in) :
    """
    search in directory file with dataset
    return : list of files
    """
    l_out = []
    list_files = listdir(directory_in)
    for filin in list_files : 
        if search("dataset", filin) and search (".txt", filin) and not search ("RMN", filin) and not search ("NMR", filin) and not search ("OUT", filin) and not search ("~", filin): 
            if path.getsize(directory_in + filin) > 5 : 
                l_out.append (directory_in + filin)
    
    return l_out


def retrieveSummaryFile (strut, name_dataset) : 
    
    l_out = []
    path_dir = result(name_dataset)
    
    l_file_folder  = listdir(path_dir)
    
    for file_folder in l_file_folder : 
        if search("^[1-9]", file_folder) or search("OUT", file_folder) :
            path_retrieve = path_dir + file_folder + "/" + strut + "/" + "summary" + strut
            
            if path.isfile(path_retrieve) : 
#                 print path_retrieve, "OUUTTTT"
                l_out.append (path_retrieve)
                
    return l_out


def imposeNeighbors (p_dir_result) : 
    
    rep = p_dir_result + "superimpose/"
    try: makedirs(rep, mode=0777)
    except: pass
    return rep


def bondLength (p_dir_result) : 
    
    rep = p_dir_result + "bondLength/"
    try: makedirs(rep, mode=0777)
    except: pass
    return rep


def coplorIMD (p_dir_result) : 
    
    rep = p_dir_result + "coplarIMD/"
    try: makedirs(rep, mode=0777)
    except: pass
    return rep


def coplorGUA (p_dir_result) : 
    
    rep = p_dir_result + "coplarGUA/"
    try: makedirs(rep, mode=0777)
    except: pass
    return rep


def twoArea(pr_int, name_in = "") : 
    
    
    if name_in == "" : 
        pr_area = pr_int + "splitArea/"
    else :
        pr_area = pr_int + "splitArea/" + name_in + "/"
        
    if not path.isdir(pr_area) : 
        makedirs(pr_area, mode=0777)
    else : 
        pass
    return pr_area
    

def combination (pr_int, name_in = "") : 
    
    if name_in == "" : 
        pr_combinatoire = pr_int + "combination/"
    else :
        pr_combinatoire = pr_int + "combination/" + name_in + "/"
        
    if not path.isdir(pr_combinatoire) : 
        makedirs(pr_combinatoire, mode=0777)
    else : 
        pass
    return pr_combinatoire


def typeSubStructure (directory_in, element) : 
     
     
    if directory_in[-1] == "/" : 
        directory_in = directory_in[0:-1]
    else :
        pass
    dir_out = directory_in + "/" + element + "/"
    try: makedirs(dir_out, mode=0777)
    except: pass
    return dir_out
 
def CreatePathDir (pr_in):
     
    try : makedirs( pr_in, mode=0777 )
    except : pass
     
    return pr_in


def ResultAngleCriteria (pr_in, name_in) : 
    
    
    if name_in == "" : 
        pr_out = pr_in + "AngleSelect/"
    else :
        pr_out = pr_in + "AngleSelect/" + name_in + "/"
        
    if not path.isdir(pr_out) : 
        makedirs(pr_out, mode=0777)
    else : 
        pass
    return pr_out
    
    
    
    
    

