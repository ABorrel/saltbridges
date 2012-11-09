from urllib import urlretrieve
from os import system

def retrieveFirstSeq(pdbSerial, repertoryFasta):
    """Retieve in pdb site sequence of pdb -> fasta format
    in: repertory out fasta, PDB serial
    out: empty, fill repertory fasta"""

    adresseSeq = ("http://www.pdb.org/pdb/files/fasta.txt?structureIdList=%s" % pdbSerial)

    try:
        pathFileFasta = urlretrieve(adresseSeq)
        print pathFileFasta
        
        nameFasta = pdbSerial + ".fasta"
        
        cmd = "mv " + pathFileFasta[0] + " " + repertoryFasta + nameFasta
        print cmd
        system (cmd)
        
    except:
        print "Impossible retrieve"
