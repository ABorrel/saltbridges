from re import search
from urllib import urlretrieve
from copy import deepcopy

import loadFile
import parsing
import searchPDB


def retrieveFirstSeq(pdb):
    """Retieve in pdb site sequence of pdb
    retrieve first sequence of file"""

    adresseSeq = ("http://www.pdb.org/pdb/files/fasta.txt?structureIdList=%s" % pdb)

    try:
        pathFileFasta = urlretrieve(adresseSeq)
        sequence = retrieveFirstSeqFileFasta(pathFileFasta)
    except:
        print "Impossible retrieve"
        return ""
    
    return sequence

    

def retrieveFirstSeqFileFasta(pathFileFasta):
    """Retrieve the first sequence of fasta file
    in : fadta file in temp repertory
    out : sequences of amino acid -> string format -> no condition -> empty string"""

    fileFasta = open(pathFileFasta[0], "r")
    linesRead = fileFasta.readlines()
    fileFasta.close()

    lenFile = len(linesRead)
    sequence = ""

    i = 1
    while i != lenFile:
        if search("^>", linesRead[i]):
            i = lenFile
        else:
            seqLine = linesRead[i].split("\n")[0]
            sequence = sequence + seqLine
            i = i + 1

    return sequence



def CheckComplexQuality(l_in, name_lig, limit_RX, limit_RFree, one_PDB_out, debug = 1):
    """Check if the PDB file is similar and preserve file with the best Quality
    remove file contain only DNA, RNA structure
    in : list of PDB files
    out : void -> change directly PDB list"""

    l_PDB = deepcopy(l_in)
    
    nb_PDB = len(l_PDB)
    l_RX = []

    i = 0
    while i < nb_PDB:
        
        # load PDB
        d_PDB_load = loadFile.ExtractInfoPDBID(l_PDB[i])
        
        # case PDB not found => often difference between list and database
        if d_PDB_load == {} :
            del l_PDB[i]
            nb_PDB = nb_PDB - 1
            continue  
        
        # check Header
        header_PDB = d_PDB_load["Header"]
        # Check DNA or RNA
        if search ("dna", header_PDB) or search ("rna", header_PDB):
            if debug == 1 : print "Exit => out dna", l_PDB[i]
            del l_PDB[i]
            nb_PDB = nb_PDB - 1
            continue
        
        # check quality
        RX = d_PDB_load["RX"]
        RFree = d_PDB_load["RFree"]
        if float (RX) > limit_RX or float (RFree) > limit_RFree : 
            if debug == 1 : print "Exit => Quality structure", l_PDB[i], RX, RFree
            del l_PDB[i]
            nb_PDB = nb_PDB - 1
            continue
         
        # check if the structure contain the substructure
        # -> because when we improve the resolution the ligand change and the search PDB -> ERROR
        l_interest_sub = searchPDB.interestStructure(d_PDB_load[name_lig][0])
        if l_interest_sub == [] : 
            if debug == 1 : print "Exit => No Sub found", l_PDB[i], RX, RFree
            del l_PDB[i]
            nb_PDB = nb_PDB - 1
            continue
            
        # Check false ligand, hooked to the protein
        l_atom_protein = d_PDB_load["protein"]
        l_atom_lig = d_PDB_load[name_lig][0]
        
        if parsing.checkLigandHooked (l_atom_protein, l_atom_lig) == 1:
            if debug == 1 : print "Exit => ligand = Modificated Amino acid", l_PDB[i]
            nb_PDB = nb_PDB - 1
            del l_PDB[i]
            continue
        
        l_RX.append (RX)
        i = i + 1
    
    
    
    if len (l_PDB) == 1 or len (l_PDB) == 0 :
        return l_PDB
    
    # case where we want only one PDB -> retrieve best quality
    if one_PDB_out == 1 : 
        return [l_PDB[l_PDB.index (min (l_PDB))]]
    
    return l_PDB




def SelectBestComplexRX (l_PDB) : 
    """
    Retrieve best PDB among list of PDB
    """
    
    l_RX = []
    l_Rfree = []
    
    for PDB_ID in l_PDB : 
        l_quality =  parsing.Quality(PDB_ID)
        l_RX.append (l_quality[0])
        l_Rfree.append (l_quality[1])
    
    # based on Resolution
    i_best_PDB = l_RX.index (max (l_RX))
    
    return [l_PDB[i_best_PDB]]

    
    

