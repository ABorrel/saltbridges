from os import makedirs

#globals()["repInit"] = "/home/student10/stage/"
#globals()["repInit"] = "/home/alexandre/intership/"
globals()["repInit"] = "/home/alexandre/project/saltBridgesProject/"




def result():

    rep = repInit + "result/"
    try: makedirs(rep, mode=0777)
    except: pass
    return rep

def withoutAtLeastOneSummary():

    rep = result() + "withoutAtLeastOne/"
    try: makedirs(rep, mode=0777)
    except: pass
    return rep



def openPdbFile ():

    rep = repInit + "PDB/"
    try: makedirs(rep, mode=0777)
    except: pass
    return rep

def pdbechemFile ():

    rep = repInit + "PDBeChem/pdb/"
    try: makedirs(rep, mode=0777)
    except: pass
    return rep

def logFile () :

    rep = repInit + "result/log/"
    try: makedirs(rep, mode=0777)
    except: pass
    return rep

def aminoAcid(structure):

    rep = resultStruct(structure)
    rep = rep + "aminoAcid/"
    try: makedirs(rep, mode=0777)
    except: pass
    return rep

def resultDistance():

    rep = repInit + "result/result_distance/"
    try: makedirs(rep, mode=0777)
    except: pass
    return rep

def resultStruct(structure):

    rep = repInit + "result/" + structure +"/"
    try: makedirs(rep, mode=0777)
    except: pass
    return rep



def scriptR():
    rep =  "./../../"
    return rep

def resultAminoAcidGlobal():

    rep = repInit + "result/AminoAcidGlobal/"
    try: makedirs(rep, mode=0777)
    except: pass
    return rep

def globalProportionAtom() :

    rep = repInit + "result/globalProportionAtom/"
    try: makedirs(rep, mode=0777)
    except: pass
    return rep

def globalProportionType() :

    rep = repInit + "result/globalProportionType/"
    try: makedirs(rep, mode=0777)
    except: pass
    return rep


def resultAngle(type):

    rep = repInit + "result/angle/"
    try: makedirs(rep, mode=0777)
    except: pass
    
    rep = rep + type + "/"
    try: makedirs(rep, mode=0777)
    except: pass
    return rep

def parsingDataset():
    
    rep = repInit + "result/parsingDataset/"
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

    rep = repInit + "dataset_SerineProtease/"+rep + "/"
    return rep    
