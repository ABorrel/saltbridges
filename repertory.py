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

def pdbechemFile ():

    rep = repInit + "PDBeChem/pdb/"
    try: makedirs(rep, mode=0777)
    except: pass
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


def lengthBond (rep_int):

    rep = rep_int + "lengthBond/"
    try: makedirs(rep, mode=0777)
    except: pass
    return rep


def resultStruct(structure, dir_out):

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
        pr_angle = pr_int + "countProx/"
    else :
        pr_angle = pr_int + "countProx/" + name_in + "/"
        
    if not path.isdir(pr_angle) : 
        makedirs(pr_angle, mode=0777)
    else : 
        pass
    return pr_angle






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
    list_out = []
    list_files = listdir(directory_in)
    for filin in list_files : 
        if search("dataset", filin) and not search (".stat", filin) and not search ("RMN", filin) and not search ("NMR", filin) and not search ("OUT", filin): 
            if path.getsize(directory_in + filin) > 5 : 
                list_out.append (directory_in + filin)
    
    return list_out


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

#############################Serine Protease###################3


def datasetSerineProtease ():

    rep = repInit + "dataset_SerineProtease/"
    return rep


def resultSerineProtease ():
    
    rep = repInit + "resultSerineProtein/"
    try: makedirs(rep, mode=0777)
    except: pass
    return rep

def datasetSerineProteaseType (rep):

    rep = repInit + "dataset_SerineProtease/" + rep + "/"
    return rep    


def typeSubStructure (directory_in, element) : 
    
    
    if directory_in[-1] == "/" : 
        directory_in = directory_in[0:-1]
    else :
        pass
    dir_out = directory_in + "/" + element + "/"
    try: makedirs(dir_out, mode=0777)
    except: pass
    return dir_out


