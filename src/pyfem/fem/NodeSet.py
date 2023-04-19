import re

from numpy import array

from pyfem.utils.parser import getType
from pyfem.utils.itemList import itemList
from pyfem.utils.logger import getLogger

logger = getLogger()


# -------------------------------------------------------------------------------
#
# -------------------------------------------------------------------------------

class NodeSet(itemList):

    def __init__(self):
        self.rank = -1
        self.groups = {}

    # -------------------------------------------------------------------------------
    #
    # -------------------------------------------------------------------------------

    def getNodeCoords(self, nodeIDs):
        return array(self.get(nodeIDs))

    # -------------------------------------------------------------------------------
    #
    # -------------------------------------------------------------------------------

    def readFromFile(self, fname):

        logger.info("Reading nodes ................")

        fin = open(fname, 'r')

        line = fin.readline()

        while line:
            if line.replace(" ", "").startswith('<Nodes>'):
                self.readNodalCoords(fin)

            if line.replace(" ", "").startswith('gmsh'):
                ln = line.replace('\n', '').replace('\t', '').replace(' ', '').replace('\r', '').replace(';', '')
                ln = ln.split('=', 1)
                self.readGmshFile(ln[1][1:-1])
                break

            line = fin.readline()

        fin = open(fname, 'r')

        line = fin.readline()

        while line:
            if line.replace(" ", "").startswith('<NodeGroup'):
                if 'name' in line:
                    label = line.split('=')[1].replace('\n', '').replace('>', '').replace(' ', '').replace('\"',
                                                                                                           '').replace(
                        '\'', '')
                    self.readNodegroup(fin, label)

            line = fin.readline()

        for key in self.groups:
            self.groups[key] = list(set(self.groups[key]))

    # -------------------------------------------------------------------------------
    #
    # -------------------------------------------------------------------------------

    def readGmshFile(self, fname):

        import meshio

        mesh = meshio.read(fname, file_format="gmsh")

        obj3d = ["pris", "pyra", "hexa", "wedg", "tetr"]

        self.rank = 2

        for key in mesh.cell_sets_dict:
            for typ in mesh.cell_sets_dict[key]:
                if (typ[:4] in obj3d):
                    self.rank = 3

        for nodeID, p in enumerate(mesh.points):
            self.add(nodeID, p[:self.rank])

        for key in mesh.cell_sets_dict:
            if key == "gmsh:bounding_entities":
                pass
            else:
                for typ in mesh.cell_sets_dict[key]:
                    for idx in mesh.cell_sets_dict[key][typ]:
                        iNodes = mesh.cells_dict[typ][idx]
                        for nodeID in iNodes:
                            self.addToGroup(key, nodeID)

    # -------------------------------------------------------------------------------
    #
    # -------------------------------------------------------------------------------

    def addToGroup(self, modelType, ID):

        if modelType not in self.groups:
            self.groups[modelType] = [int(ID)]
        else:
            self.groups[modelType].append(int(ID))

    # -------------------------------------------------------------------------------
    #
    # -------------------------------------------------------------------------------

    def __repr__(self):
        msg = "Number of nodes ............ %6d\n" % len(self)

        if len(self.groups) > 0:
            msg += "  Number of  groups .......... %6d\n" % len(self.groups)
            msg += "  -----------------------------------\n"
            msg += "    name                       #nodes\n"
            msg += "    ---------------------------------\n"

            for name in self.groups:
                msg += "    %-16s           %6d \n" % (name, len(self.groups[name]))

        return msg

    # -------------------------------------------------------------------------------
    #
    # -------------------------------------------------------------------------------

    def readNodalCoords(self, fin):

        while True:
            line = fin.readline()

            if line.replace(" ", "").startswith('</Nodes>'):
                return

            line = re.sub('\s{2,}', ' ', line)
            a = line.split(';')

            for a in a[:-1]:
                b = a.strip().split(' ')

                if b[0].startswith("//") or b[0].startswith("#"):
                    break
                if len(b) > 1 and type(eval(b[0])) == int:
                    if self.rank == -1:
                        self.rank = len(b) - 1

                    self.add(eval(b[0]), [eval(crd) for crd in b[1:]])

                # -------------------------------------------------------------------------------

    #
    # -------------------------------------------------------------------------------

    def readNodegroup(self, fin, key):

        while True:
            line = fin.readline()

            if line.replace(" ", "").startswith('</NodeGro'):
                return

            a = line.split()

            for b in a:
                if getType(b) == int:
                    self.addToGroup(key, b)
