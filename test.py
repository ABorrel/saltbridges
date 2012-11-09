import tool
import checkPDBfile
import parsing
import retrieveAtom
import loadFile
import calcul
import structure
import searchPDB
import repertory



def testDistance():
    '''calculate the distance between two atom'''

    fileListeLigands = open ("list_ligands", "r")
    lectfileListeLigands = fileListeLigands.readlines()
    fileListeLigands.close()


    nameLigands = manageligandsliste(lectfileListeLigands, 0)
    print nameLigands
    File = "./pdb/" + nameLigands + ".pdb"
    ATOM = parsePDB(File)

    print ATOM
    for i in ATOM:
        for j in i['CONECT']:
            distance = distanceTwoPoints(i, ATOM[int (j) - 1])
            print "distance", i['compound'], ATOM[int(j) - 1]['compound'], distance


def testcoplanar():
    file = "/home/student10/stage/PDBeChem/pdb/10C.pdb"
    ATOM = parsing.parsePDBechem(file)
    print ATOM

    listN = substructure.searchN(ATOM)

    for N in listN:
        tool.coplanar(ATOM[int(N) - 1], ATOM)

    print tool.coplanar(ATOM[1], ATOM)

def testCompareSeq(pdbFile1, pdbFile2) :

    list = []
    out = checkPDBfile.compareSeqofPDB(pdbFile1, pdbFile2, list)

    print out
    #print h1["name"],h2["name"]
    #print h1["resolution"],h2["resolution"]


def testRetrieveligand (pdbName, nameLigands) :

    coords = coordinatePDB.loadLigandInPDB(pdbName, nameLigands)
    print "in"
    print coords



def testAngleTertiary():

    ligand = loadFile.globalPDB("2v3d", "NBV")
    print ligand

    nitrogen = retrieveAtom.serial(7930, ligand)
    print nitrogen
    counterIon = retrieveAtom.serial(1842, ligand)

    print calcul.angleTertiaryAmine(nitrogen, counterIon, ligand)


def testAngleSecondary():

    ligand = loadFile.globalPDB("2j0k", "4ST")
    print ligand

    nitrogen = retrieveAtom.serial(9789, ligand)
    print nitrogen
    counterIon = retrieveAtom.serial(3499, ligand)

    print calcul.angleSecondaryAmine(nitrogen, counterIon, ligand)


def testSearch():
    
    ligandInPDB = loadFile.resultLigandPDB(repertory.result() + "resultLigandInPDB")
    resultFilterPDB = structure.resolutionFilter()

    nbLigand = len(ligandInPDB.keys())
    listligand = ligandInPDB.keys()
    
    i = 8065
    while (i < 8066):
        nameLigand = listligand[i]
        print nameLigand, i
        #logFile.write(str(nameLigand) + " " + str(i) + "\n")

        atomLigand = loadFile.ligandInPDBConnectMatrixLigand(ligandInPDB[nameLigand][0], nameLigand)
        listStruct = searchPDB.interestStructure(atomLigand)
        print listStruct
        i = i + 1
    
def repertoryTest():
    print repertory.result()
    
