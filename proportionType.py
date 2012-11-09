import structure


def globalNeighbors (listAtom, count):
    """Count numbers of neighbor with atom type
    in: list atom global
    out: count with keys -> numbers of neighbors"""

    for atom in listAtom:
        nbNeighbor = numberNeigthbor(atom["neighbors"])
        for neighbor in atom["neighbors"]:
            #print count
            count["allNumberNeighbors"][neighbor["classification"]] = count["allNumberNeighbors"][neighbor["classification"]] + 1
            if not nbNeighbor in count.keys():
                count[nbNeighbor] = structure.countClassificationAtoms()

            if neighbor["classification"] in count[nbNeighbor].keys():
                count[nbNeighbor][neighbor["classification"]] = count[nbNeighbor][neighbor["classification"]] + 1

            else:
                count[nbNeighbor]["others"] = count[nbNeighbor]["others"] + 1




def amine(listAmine, count):
    """Count numbers of neighbor with atom type
    in: list atom global
    out: count with keys -> numbers of neighbors"""
    
    for type in listAmine.keys():
        for nitrogen in listAmine[type]:
            nbNeighbor = numberNeigthbor(nitrogen["neighbors"])
            
            for neighbor in nitrogen["neighbors"]:
                count[type]["allNumberNeighbors"][neighbor["classification"]] = count[type]["allNumberNeighbors"][neighbor["classification"]] + 1
                
                if not nbNeighbor in count[type].keys():
                    count[type][nbNeighbor] = structure.countClassificationAtoms()
                if not nbNeighbor in count["GlobalAmine"].keys():
                    count["GlobalAmine"][nbNeighbor] = structure.countClassificationAtoms()


                if neighbor["classification"] in count[type][nbNeighbor].keys():
                    count[type][nbNeighbor][neighbor["classification"]] = count[type][nbNeighbor][neighbor["classification"]] + 1
                    count["GlobalAmine"][nbNeighbor][neighbor["classification"]] = count["GlobalAmine"][nbNeighbor][neighbor["classification"]] + 1

                else:
                    count[type][nbNeighbor]["others"] = count[type][nbNeighbor]["others"] + 1
                    count["GlobalAmine"][nbNeighbor]["others"] = count["GlobalAmine"][nbNeighbor]["others"] + 1


def numberNeigthbor(listNeighbor):
    """Count atom neighbor
    in: Structure of global atom
    out: count in string"""

    count = 0
    for neighbor in listNeighbor:
        if not neighbor == []:
            count = count + 1
    return str(count)

