from re import search
from urllib import urlretrieve

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
#     """Compare sequence of amino acid of list of pdb files pairwise sequence comparaison, preserve best resolution file and remove sequences and PDB in list
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
#         if listSeq[i] < listSeq[j]["resolution"]:
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


def checkPDB(listPDB, nameLigand, limit_RX, limit_RFree):
    """Check if the PDB file is similar and preserve file with the best resolution
    remove file contain only DNA, RNA structure
    in : list of PDB files
    out : list of PDB files"""

    nb_PDB = len(listPDB)

    i = 0
    while i < nb_PDB:
        print i
        header_PDB = parsing.header(listPDB[i])
        l_quality = parsing.resolution(listPDB[i])
        print l_quality
        if l_quality[0] > limit_RX or l_quality[1] > limit_RFree : 
            del listPDB[i]
            nb_PDB = nb_PDB - 1
        elif search ("dna", header_PDB) or search ("rna", header_PDB):
            del listPDB[i]
            nb_PDB = nb_PDB - 1
            print "out dna", nb_PDB
        else :
            l_atom_complex = loadFile.globalPDB(listPDB[i])
            # only first ligand included in PDB
            l_atom_ligand = parsing.retrieveLigand (l_atom_complex, nameLigand)[0] 
            if parsing.checkLigandHooked (l_atom_complex, l_atom_ligand) == 1:
                nb_PDB = nb_PDB - 1
                del listPDB[i]
            else : 
                i = i + 1

    return



def selectBestPDBamongList (l_PDB) : 
    
    l_RX = []
    l_Rfree = []
    
    for PDB_ID in l_PDB : 
        l_quality =  parsing.resolution(PDB_ID)
        l_RX.append (l_quality[0])
        l_Rfree.append (l_quality[1])
    
    i_best_PDB = l_RX.index (max (l_RX))
    
    return [l_PDB[i_best_PDB]]

    
    

