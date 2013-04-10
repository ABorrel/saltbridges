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



######compare#######


def compareSeqofPDB(listPDB, listSeq, lengthListPDB, i, j):
    """Compare sequence of amino acid of list of pdb files pairwise sequence comparaison, preserve best resolution file and remove sequences and PDB in list
    in : - list PDB file
         - list of amino acid sequences
         - number of PDB in list
         - position files comparaison in list
    out : - list pdb file
          - list of amino acid sequences
          - number of PDB file in list PDB after comparaison
          - position files comparaison in list"""


    if i == j:
        return listPDB, listSeq, lengthListPDB, i, j

    if listSeq[i]["seq"] == "" or listSeq[j]["seq"] == "":
        return listPDB, listSeq, lengthListPDB, i, j

    elif listSeq[i]["seq"] == listSeq[j]["seq"]:
        if listSeq[i] < listSeq[j]["resolution"]:
            del listPDB[j]
            del listSeq[j]
            j = j - 1
            lengthListPDB = lengthListPDB - 1
        else:
            del listPDB[i]
            del listSeq[i]
            i = i - 1
            lengthListPDB = lengthListPDB - 1

        return listPDB, listSeq, lengthListPDB, i, j

    elif incluedSequence(listSeq[i], listSeq[j]) == 1:
        del listPDB[j]
        del listSeq[j]
        j = j - 1
        lengthListPDB = lengthListPDB - 1

        return listPDB, listSeq, lengthListPDB, i, j

    elif incluedSequence(listSeq[j], listSeq[i]) == 1:
        del listPDB[i]
        del listSeq[i]
        i = i - 1
        lengthListPDB = lengthListPDB - 1

        return listPDB, listSeq, lengthListPDB, i, j

    else:

        return listPDB, listSeq, lengthListPDB, i, j


def incluedSequence(seq1, seq2):
    """Function which look if seq1 is inclued in seq2
    in : - seq1
         - seq2
    seq 1 bigger seq 2
    out : 0 or 1
          0 -> no inclued
          1 -> inclued"""

    try:
        if search(seq1, seq2):
            return 1
        else:
            return 0

    except:
        return 0


def checkPDB(listPDB, nameLigand):
    """Check if the PDB file is similar and preserve file with the best resolution
    remove file contain only DNA, RNA structure
    in : list of PDB files
    out : list of PDB files"""

    listSeq = []
    lengthListPDB = len(listPDB)

    if lengthListPDB == 1:
        if ligandHooked(nameLigand, listPDB[0]) == 1:
            del listPDB[0]
            return
        
        PDBType = parsing.methods(listPDB[0])
        if PDBType == "dna" or PDBType == "rna" or PDBType == "dna-rna":
            print "out dna"
            del listPDB[0]
            return
        else:
            return
    else:
        
        i = 0
        while i < lengthListPDB:
            if ligandHooked(nameLigand, listPDB[i]) == 1:
                lengthListPDB = lengthListPDB - 1
                del listPDB[i]
                continue

            PDBType = parsing.methods(listPDB[i])

            if PDBType == "dna" or PDBType == "rna" or PDBType == "dna-rna":
                del listPDB[i]
                lengthListPDB = lengthListPDB - 1
                print "out dna", lengthListPDB
            else:
                Sequence = {}
                Sequence["name"] = listPDB[i]
                Sequence["seq"] = retrieveFirstSeq(listPDB[i])
                Sequence["resolution"] = parsing.resolution(listPDB[i])
                listSeq.append(Sequence)
                i = i + 1
    i = 0
    while i < lengthListPDB:
        j = i + 1
        while j < lengthListPDB:
            tempi = i
            listPDB, listSeq, lengthListPDB, i, j = compareSeqofPDB(listPDB, listSeq, lengthListPDB, i, j)
            j = j + 1
            if tempi != i:
                break
        i = i + 1

    return




def ligandHooked(nameLigand, filePDB):
    """Check if ligand is hook in protein
    in : ligand -> list of atom (dictionnary)
    out : 0 -> no hook  1 -> hook"""

    ligand = loadFile.ligandInPDBConnectMatrixLigand(filePDB, nameLigand)

    listSerial = []
    for atom in ligand:
        listSerial.append(atom["serial"])

    for atom in ligand:
        for connect in atom["connect"]:
            if not connect in listSerial:
                print "Ligand Hooked", filePDB
                return 1
    return 0
