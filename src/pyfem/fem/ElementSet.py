import re

from pyfem.utils.data_structures import SolverStatus
from pyfem.utils.item_list import ItemList
from pyfem.utils.logger import getLogger

logger = getLogger()


class ElementSet(ItemList):

    def __init__(self, nodes, props):

        ItemList.__init__(self)

        self.nodes = nodes
        self.props = props
        self.solverStat = SolverStatus()
        self.groups = {}



    def __iter__(self):

        elements = []

        for groupName in self.iterGroupNames():
            for element in self.iterElementGroup(groupName):
                elements.append(element)

        return iter(elements)



    def __repr__(self):
        str_ = "Number of elements ......... %6d\n" % len(self)

        if len(self.groups) > 0:
            str_ += "  Number of  groups .......... %6d\n" % len(self.groups)
            str_ += "  -----------------------------------\n"
            str_ += "    name                       #elems\n"
            str_ += "    ---------------------------------\n"
            for name in self.groups:
                str_ += "    %-16s           %6d\n" % (name, len(self.groups[name]))

        return str_



    def getDofTypes(self):

        dofTypes = []

        for element in self:
            for dof_type in element.dofTypes:
                if dof_type not in dofTypes:
                    dofTypes.append(dof_type)

        return dofTypes



    def read_from_file(self, fname):

        logger.info("Reading elements .............")

        fin = open(fname)

        while True:
            line = fin.readline()

            if line.startswith('<Elements>') == True:
                while True:
                    line = fin.readline()

                    if line.startswith('</Elements>') == True:
                        return

                    line = re.sub('\s{2,}', ' ', line)
                    a = line.split(';')

                    for a0 in a[:-1]:
                        b = a0.strip().split(' ')

                        if b[0].startswith("//") or b[0].startswith("#"):
                            break
                        if len(b) > 1 and type(eval(b[0])) == int:
                            self.add(eval(b[0]), eval(b[1]), [eval(node_id) for node_id in b[2:]])

            elif line.startswith('gmsh') == True:
                ln = line.replace('\n', '').replace('\t', '').replace(' ', '').replace('\r', '').replace(';', '')
                ln = ln.split('=', 1)
                self.read_gmsh_file(ln[1][1:-1])
                return



    def read_gmsh_file(self, fname):

        import meshio

        mesh = meshio.read(fname, file_format="gmsh")

        elemID = 0

        for key in mesh.cell_sets_dict:
            for typ in mesh.cell_sets_dict[key]:
                for idx in mesh.cell_sets_dict[key][typ]:
                    iNodes = mesh.cells_dict[typ][idx]
                    self.add(elemID, key, iNodes.tolist())
                    elemID = elemID + 1

    # -------------------------------------------------------------------------------
    #  add element
    # -------------------------------------------------------------------------------

    def add(self, ID, modelName, elementNodes):

        # Check if the model exists

        if hasattr(self.props, modelName):

            modelProps = getattr(self.props, modelName)

            # Check if the model has a type
            if not hasattr(modelProps, 'type'):
                raise RuntimeError('Missing type for model ' + modelName)

            model_type = getattr(modelProps, 'type')

            modelProps.rank = self.nodes.rank
            modelProps.solverStat = self.solverStat

            element = getattr(__import__('pyfem.elements.' + model_type, globals(), locals(), model_type, 0), model_type)

            # Create the element

            elem = element(elementNodes, modelProps)

            #  Check if the node ids are valid:

            for node_id in elem.getNodes():
                if not node_id in self.nodes:
                    raise RuntimeError('Node ID ' + str(node_id) + ' does not exist')

            #  Add the element to the element set:

            ItemList.add(self, ID, elem)

            #  Add the element to the correct group:

            self.add_to_group(modelName, ID)



    def add_to_group(self, model_type, ID):

        if model_type not in self.groups:
            self.groups[model_type] = [ID]
        else:
            self.groups[model_type].append(ID)



    def addGroup(self, groupName, groupIDs):
        self.groups[groupName] = groupIDs



    def iterGroupNames(self):
        return self.groups



    def iterElementGroup(self, groupName):
        if groupName == "All":
            return iter(self)
        elif isinstance(groupName, list):
            elems = []
            for name in groupName:
                elems += self.get(self.groups[name])
            return iter(elems)
        else:
            return iter(self.get(self.groups[groupName]))



    def elementGroupCount(self, groupName):
        if groupName == "All":
            return len(self)
        elif isinstance(groupName, list):
            length = 0;
            for name in groupName:
                length += len(self.groups[name])
            return length
        else:
            return len(self.groups[groupName])

    #
    #
    #

    def getFamilyIDs(self):

        familyIDs = []
        fam = ["CONTINUUM", "INTERFACE", "SURFACE", "BEAM", "SHELL"]

        for elem in self:
            familyIDs.append(fam.index(elem.family))

        return familyIDs



    def commitHistory(self):

        for element in list(self.values()):
            element.commitHistory()
