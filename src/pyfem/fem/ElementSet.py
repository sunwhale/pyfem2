import re
import os

from pyfem.utils.data_structures import SolverStatus, Properties
from pyfem.utils.IntegerIdDict import IntegerIdDict
from pyfem.utils.logger import get_logger
from pyfem.fem.NodeSet import NodeSet

logger = get_logger()

ELEMENTS_START = "<Elements>"
ELEMENTS_END = "</Elements>"
GROUP_START = "<ElementGroups"
GROUP_END = "</ElementGroups>"
GMSH_START = "gmsh"
COMMENT_STARTERS = ["//", "#"]


class ElementSet(IntegerIdDict):

    def __init__(self, nodes: NodeSet, props: Properties):

        IntegerIdDict.__init__(self)

        self.nodes = nodes
        self.props = props
        self.solver_status = SolverStatus()
        self.groups = {}

    def __iter__(self):
        for group_name in self.iter_group_names():
            for element in self.iter_element_group(group_name):
                yield element

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

    def get_dof_types(self):
        dof_types = {dof_type for element in self for dof_type in element.dof_types}
        return list(dof_types)

    def read_from_file(self, file_name: str) -> None:

        logger.info("Reading elements .............")

        # with open(file_name, 'r') as f:
        #     for line in f:
        #         line = line.strip().replace(' ', '').replace('\n', '').replace('\t', '').replace('\r', '')
        #         if line.startswith(GMSH_START):  # If the line starts with GMSH_START, read a Gmsh file
        #             clean_line = line.replace(';', '').replace('\"', '').replace('\'', '')
        #             gmsh_file_name = clean_line.split('=')[1]
        #             self.read_gmsh_file(gmsh_file_name)
        #         elif line.startswith(ELEMENTS_START):  # If the line starts with NODES_START, read the node coordinates
        #             self.read_node_coords(f)
        #         elif line.startswith(GROUP_START):  # If the line starts with GROUP_START, read a node group
        #             if 'name' in line:
        #                 clean_line = line.replace('>', '').replace('\"', '').replace('\'', '')
        #                 label = clean_line.split('=')[1]
        #                 self.read_node_group(f, label)
        #
        #     # Remove duplicate entries in each group
        #     for key in self.groups:
        #         self.groups[key] = list(set(self.groups[key]))

        fin = open(file_name)

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
                            self.add_item_by_id(eval(b[0]), eval(b[1]), [eval(node_id) for node_id in b[2:]])

            elif line.startswith('gmsh') == True:
                ln = line.replace('\n', '').replace('\t', '').replace(' ', '').replace('\r', '').replace(';', '')
                ln = ln.split('=', 1)
                self.read_gmsh_file(ln[1][1:-1])
                return

    def read_gmsh_file(self, file_name):

        import meshio

        mesh = meshio.read(file_name, file_format="gmsh")

        elemID = 0

        for key in mesh.cell_sets_dict:
            for typ in mesh.cell_sets_dict[key]:
                for idx in mesh.cell_sets_dict[key][typ]:
                    iNodes = mesh.cells_dict[typ][idx]
                    self.add_item_by_id(elemID, key, iNodes.tolist())
                    elemID = elemID + 1


    def add_item_by_id(self, ID, modelName, elementNodes):

        # Check if the model exists

        if hasattr(self.props, modelName):

            modelProps = getattr(self.props, modelName)

            # Check if the model has a type
            if not hasattr(modelProps, 'type'):
                raise RuntimeError('Missing type for model ' + modelName)

            model_type = getattr(modelProps, 'type')

            modelProps.rank = self.nodes.rank
            modelProps.solver_status = self.solver_status

            element = getattr(__import__('pyfem.elements.' + model_type, globals(), locals(), model_type, 0),
                              model_type)

            # Create the element

            elem = element(elementNodes, modelProps)

            #  Check if the node ids are valid:

            for node_id in elem.getNodes():
                if not node_id in self.nodes:
                    raise RuntimeError('Node ID ' + str(node_id) + ' does not exist')

            #  Add the element to the element set:

            IntegerIdDict.add_item_by_id(self, ID, elem)

            #  Add the element to the correct group:

            self.add_to_group(modelName, ID)

    def add_to_group(self, model_type, ID):

        if model_type not in self.groups:
            self.groups[model_type] = [ID]
        else:
            self.groups[model_type].append(ID)

    def add_group(self, groupName, groupIDs):
        self.groups[groupName] = groupIDs

    def iter_group_names(self):
        return self.groups

    def iter_element_group(self, groupName):
        if groupName == "All":
            return iter(self)
        elif isinstance(groupName, list):
            elems = []
            for name in groupName:
                elems += self.get_items_by_ids(self.groups[name])
            return iter(elems)
        else:
            return iter(self.get_items_by_ids(self.groups[groupName]))

    def element_group_count(self, groupName):
        if groupName == "All":
            return len(self)
        elif isinstance(groupName, list):
            length = 0;
            for name in groupName:
                length += len(self.groups[name])
            return length
        else:
            return len(self.groups[groupName])

    def get_family_ids(self):

        familyIDs = []
        fam = ["CONTINUUM", "INTERFACE", "SURFACE", "BEAM", "SHELL"]

        for elem in self:
            familyIDs.append(fam.index(elem.family))

        return familyIDs

    def commit_history(self):

        for element in list(self.values()):
            element.commit_history()


if __name__ == "__main__":
    from pyfem.utils.parser import file_parser
    os.chdir('F:\\Github\\pyfem\\examples\\gmsh\\rectangle')
    props = file_parser('rectangle.pro')
    nset = NodeSet()
    nset.read_from_file('rectangle.dat')
    elset = ElementSet(nset, props)
    elset.read_from_file('rectangle.dat')
    print(elset.items())

    os.chdir('F:\\Github\\pyfem\\examples\\mesh')
    props = file_parser('PatchTest8_3D.pro')
    nset = NodeSet()
    nset.read_from_file('PatchTest8_3D.dat')
    elset = ElementSet(nset, props)
    elset.read_from_file('PatchTest8_3D.dat')
    print(elset.items())