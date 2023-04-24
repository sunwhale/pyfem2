import os
import re
from typing import Dict, List, Union, TextIO, Iterator

import meshio

from pyfem.fem.NodeSet import NodeSet
from pyfem.utils.IntegerIdDict import IntegerIdDict
from pyfem.utils.data_structures import SolverStatus, Properties
from pyfem.utils.logger import get_logger

logger = get_logger()

ELEMENTS_START = "<Elements>"
ELEMENTS_END = "</Elements>"
GMSH_START = "gmsh"
COMMENT_STARTERS = ["//", "#"]


class ElementSet(IntegerIdDict):

    def __init__(self, nodes: NodeSet, props: Properties):
        super().__init__()
        self.nodes = nodes
        self.props = props
        self.solver_status = SolverStatus()
        self.groups = {}

    def __iter__(self) -> iter:
        for group_name in self.iter_group_names():
            for element in self.iter_element_group(group_name):
                yield element

    def __repr__(self) -> str:
        str_ = "Number of elements ......... %6d\n" % len(self)

        if len(self.groups) > 0:
            str_ += "  Number of  groups .......... %6d\n" % len(self.groups)
            str_ += "  -----------------------------------\n"
            str_ += "    name                       #elems\n"
            str_ += "    ---------------------------------\n"
            for name in self.groups:
                str_ += "    %-16s           %6d\n" % (name, len(self.groups[name]))

        return str_

    def get_dof_types(self) -> List[str]:
        dof_types = {dof_type for element in self for dof_type in element.dof_types}
        return list(dof_types)

    def read_from_file(self, file_name: str) -> None:

        logger.info("Reading elements .............")

        with open(file_name, 'r') as f:
            for line in f:
                line = line.strip().replace(' ', '').replace('\n', '').replace('\t', '').replace('\r', '')
                if line.startswith(GMSH_START):  # If the line starts with GMSH_START, read a Gmsh file
                    clean_line = line.replace(';', '').replace('\"', '').replace('\'', '')
                    gmsh_file_name = clean_line.split('=')[1]
                    self.read_gmsh_file(gmsh_file_name)
                elif line.startswith(ELEMENTS_START):  # If the line starts with NODES_START, read the node coordinates
                    self.read_element_connectivity(f)

    def read_element_connectivity(self, f: TextIO) -> None:
        while True:
            line = f.readline().strip().replace('\"', '').replace('\'', '')  # Remove leading/trailing white space

            if line.replace(' ', '').startswith(ELEMENTS_END):
                # If the line starts with NODES_END, return from the function
                return

            items = [x for x in re.split(r';|\s', line) if x]  # Split the line into items

            first_item = items[0] if items else ''
            if any(first_item.startswith(comment) for comment in COMMENT_STARTERS):
                # If the line is a comment, skip it
                continue

            if items and items[0].isdigit():
                # If the first item is a digit, add the element types and connectivity to the ElementSet object
                self.add_item_by_element_id(int(items[0]), str(items[1]), [int(node_id) for node_id in items[2:]])

    def read_gmsh_file(self, file_name: str) -> None:
        mesh = meshio.read(file_name, file_format="gmsh")

        element_id = 0

        for cell_set in mesh.cell_sets_dict:
            for mesh_type in mesh.cell_sets_dict[cell_set]:
                for id_ in mesh.cell_sets_dict[cell_set][mesh_type]:
                    cell_nodes = mesh.cells_dict[mesh_type][id_]
                    self.add_item_by_element_id(element_id, cell_set, cell_nodes.tolist())
                    element_id += 1

    def add_item_by_element_id(self, id_: int, model_name: str, element_nodes: List[int]) -> None:
        if hasattr(self.props, model_name):
            model_props = getattr(self.props, model_name, None)

            if not model_props:  # Check if the model exists
                raise RuntimeError(f"Model {model_name} does not exist")

            model_type = getattr(model_props, 'type', None)
            if not model_type:
                raise RuntimeError(f"Missing type {model_type} for model {model_name}")

            model_props.rank = self.nodes.rank
            model_props.solver_status = self.solver_status

            # Import the element module and create the element
            element_module = __import__('pyfem.elements.' + model_type, globals(), locals(), model_type, 0)
            element = getattr(element_module, model_type)
            elem = element(element_nodes, model_props)  # Create the element

            # Check if the node ids are valid
            invalid_node_ids = [node_id for node_id in elem.getNodes() if node_id not in self.nodes]
            if invalid_node_ids:
                raise RuntimeError(f"Invalid node IDs: {invalid_node_ids}")

            self.add_item_by_id(id_, elem)  # Add the element to the element set
            self.add_to_group(model_name, id_)  # Add the element to the correct group

    def add_to_group(self, model_type: str, id_: int) -> None:
        if model_type not in self.groups:
            self.groups[model_type] = [id_]
        else:
            self.groups[model_type].append(id_)

    def add_group(self, group_name: str, group_ids: List[int]) -> None:
        self.groups[group_name] = group_ids

    def iter_group_names(self) -> Iterator[str]:
        return iter(self.groups)

    def iter_element_group(self, group_name: Union[str, List[str]]) -> Iterator[object]:
        if group_name == "All":
            return iter(self)
        elif isinstance(group_name, list):
            elems = [item for name in group_name for item in self.get_items_by_ids(self.groups[name])]
            return iter(elems)
        else:
            return iter(self.get_items_by_ids(self.groups[group_name]))

    def element_group_count(self, group_name: str) -> int:
        if group_name == "All":
            return len(self)
        elif isinstance(group_name, list):
            length = 0
            for name in group_name:
                length += len(self.groups[name])
            return length
        else:
            return len(self.groups[group_name])

    def get_family_ids(self) -> List[int]:

        families = ['CONTINUUM', 'INTERFACE', 'SURFACE', 'BEAM', 'SHELL']
        family_ids = [families.index(elem.family) for elem in self]

        return family_ids

    def commit_history(self) -> None:

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
    print(elset.groups)

    os.chdir('F:\\Github\\pyfem\\examples\\mesh')
    props = file_parser('PatchTest8_3D.pro')
    nset = NodeSet()
    nset.read_from_file('PatchTest8_3D.dat')
    elset = ElementSet(nset, props)
    elset.read_from_file('PatchTest8_3D.dat')
    print(elset)

    for e in elset.iter_element_group(['ContElem','ContElem2']):
        print(type(e))
    # print(IntegerIdDict.add_item_by_id)
