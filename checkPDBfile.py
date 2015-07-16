from re import search
from urllib import urlretrieve
from copy import deepcopy

import loadFile
import parsing



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



######compare####### -> need blast or gobal alignment


# def compareSeqofPDB(listPDB, listSeq, lengthListPDB, i, j):
#     """Compare sequence of amino acid of list of pdb files pairwise sequence comparaison, preserve best Quality file and remove sequences and PDB in list
#     in : - list PDB file
#          - list of amino acid sequences
#          - number of PDB in list
#          - position files comparaison in list
#     out : - list pdb file
#           - list of amino acid sequences
#           - number of PDB file in list PDB after comparaison
#           - position files comparaison in list"""
# 
# 
#     if i == j:
#         return listPDB, listSeq, lengthListPDB, i, j
# 
#     if listSeq[i]["seq"] == "" or listSeq[j]["seq"] == "":
#         return listPDB, listSeq, lengthListPDB, i, j
# 
#     elif listSeq[i]["seq"] == listSeq[j]["seq"]:
#         if listSeq[i] < listSeq[j]["Quality"]:
#             del listPDB[j]
#             del listSeq[j]
#             j = j - 1
#             lengthListPDB = lengthListPDB - 1
#         else:
#             del listPDB[i]
#             del listSeq[i]
#             i = i - 1
#             lengthListPDB = lengthListPDB - 1
# 
#         return listPDB, listSeq, lengthListPDB, i, j
# 
#     elif incluedSequence(listSeq[i], listSeq[j]) == 1:
#         del listPDB[j]
#         del listSeq[j]
#         j = j - 1
#         lengthListPDB = lengthListPDB - 1
# 
#         return listPDB, listSeq, lengthListPDB, i, j
# 
#     elif incluedSequence(listSeq[j], listSeq[i]) == 1:
#         del listPDB[i]
#         del listSeq[i]
#         i = i - 1
#         lengthListPDB = lengthListPDB - 1
# 
#         return listPDB, listSeq, lengthListPDB, i, j
# 
#     else:
# 
#         return listPDB, listSeq, lengthListPDB, i, j
# 
# 
# def incluedSequence(seq1, seq2):
#     """Function which look if seq1 is inclued in seq2
#     in : - seq1
#          - seq2
#     seq 1 bigger seq 2
#     out : 0 or 1
#           0 -> no inclued
#           1 -> inclued"""
# 
#     try:
#         if search(seq1, seq2):
#             return 1
#         else:
#             return 0
# 
#     except:
#         return 0


def CheckComplexQuality(l_in, name_lig, limit_RX, limit_RFree, one_PDB_out, debug = 0):
    """Check if the PDB file is similar and preserve file with the best Quality
    remove file contain only DNA, RNA structure
    in : list of PDB files
    out : void -> change directly PDB list"""

    l_PDB = deepcopy(l_in)
    
    nb_PDB = len(l_PDB)
    l_RX = []

    i = 0
    while i < nb_PDB:
        
        # case of obsolete PDB 
        try : header_PDB = parsing.header(l_PDB[i])
        except : 
            i = i + 1
            continue
        
        # Check DNA or RNA
        if search ("dna", header_PDB) or search ("rna", header_PDB):
            if debug : print "out dna", nb_PDB
            del l_PDB[i]
            nb_PDB = nb_PDB - 1
            continue
        else :
            l_quality = parsing.Quality(l_PDB[i])
            if l_quality[0] > limit_RX or l_quality[1] > limit_RFree : 
                del l_PDB[i]
                nb_PDB = nb_PDB - 1
                continue
            else :
                l_atom_complex = loadFile.globalPDB(l_PDB[i])
                # only first ligand included in PDB
                l_atom_ligand = parsing.retrieveLigand (l_atom_complex, name_lig)[0] 
                if parsing.checkLigandHooked (l_atom_complex, l_atom_ligand) == 1:
                    nb_PDB = nb_PDB - 1
                    del l_PDB[i]
                else : 
                    l_RX.append (l_quality[0])
                    i = i + 1
    
    
    
    if len (l_PDB) == 1 or len (l_PDB) == 0 :
         
        return l_PDB
    
    # case where we want only one PDB -> retrieve best quality
    if one_PDB_out == 1 : 
        
        return [l_PDB[l_PDB.index (min (l_PDB))]]




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

    
    

