import structure


def globalNeighbors (listAtom, count):
    """Search for each atom the neighbor and retrieve the number of neighbor and the number of each element
    in: list Atom
    out: append the value in count structure"""

    for atom in listAtom:
        nbNeighbor = numberNeigthbor(atom["neighbors"])
        for neighbor in atom["neighbors"]:
            if not nbNeighbor in count.keys():
                count[nbNeighbor] = structure.countElements()

            if neighbor["element"] in count[nbNeighbor].keys():
                count[nbNeighbor][neighbor["element"]] = count[nbNeighbor][neighbor["element"]] + 1

            else:
                count[nbNeighbor]["others"] = count[nbNeighbor]["others"] + 1




def amine(listAmine, count):
    """For amine structure search the neighbor and count the different element
    in: structure where the neighbor is count
    out: append in count structure"""
    
    for type in listAmine.keys():
        for nitrogen in listAmine[type]:
            nbNeighbor = numberNeigthbor(nitrogen["neighbors"])
            for neighbor in nitrogen["neighbors"]:
                if not nbNeighbor in count[type].keys():
                    count[type][nbNeighbor] = structure.countElements()
                if not nbNeighbor in count["GlobalAmine"].keys():
                    count["GlobalAmine"][nbNeighbor] = structure.countElements()


                if neighbor["element"] in count[type][nbNeighbor].keys():
                    count[type][nbNeighbor][neighbor["element"]] = count[type][nbNeighbor][neighbor["element"]] + 1
                    count["GlobalAmine"][nbNeighbor][neighbor["element"]] = count["GlobalAmine"][nbNeighbor][neighbor["element"]] + 1

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

