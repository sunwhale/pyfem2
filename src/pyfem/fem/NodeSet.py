import os
import re
from typing import Dict, List, Union, TextIO

import meshio
import numpy as np

from pyfem.utils.IntegerIdDict import IntegerIdDict
from pyfem.utils.logger import get_logger

logger = get_logger()

NODES_START = "<Nodes>"
NODES_END = "</Nodes>"
GROUP_START = "<NodeGroup"
GROUP_END = "</NodeGroup>"
GMSH_START = "gmsh"
COMMENT_STARTERS = ["//", "#"]


class NodeSet(IntegerIdDict):

    def __init__(self):
        super().__init__()
        self.rank: int = -1
        self.groups: Dict[str, List[int]] = {}

    def __repr__(self) -> str:
        """
        String representation of the object.

        Returns:
            str_ (str): String representation of the nodes.
        """
        str_ = "Number of nodes ............ %6d\n" % len(self)

        if len(self.groups) > 0:
            str_ += "  Number of  groups .......... %6d\n" % len(self.groups)
            str_ += "  -----------------------------------\n"
            str_ += "    name                       #nodes\n"
            str_ += "    ---------------------------------\n"

            for name in self.groups:
                str_ += "    %-16s           %6d \n" % (name, len(self.groups[name]))

        return str_

    def get_node_coords(self, node_ids: Union[int, List[int]]) -> np.ndarray[float]:
        """
        Given node ids, return the coordinates of the nodes in a numpy array.

        Parameters:
            node_ids (Union[int, List[int]]): The node ids for which to get the coordinates.

        Returns:
            np.ndarray[float]: An array of node coordinates.
        """
        return np.array(self.get_items_by_ids(node_ids))

    def read_from_file(self, file_name: str) -> None:
        """
        Read the NodeSet object from a file.

        Parameters:
            file_name (str): The name of the file to read.

        Returns:
            None
        """
        logger.info("Reading nodes ................")

        with open(file_name, 'r') as f:
            for line in f:
                line = line.strip().replace(' ', '').replace('\n', '').replace('\t', '').replace('\r', '')
                if line.startswith(GMSH_START):  # If the line starts with GMSH_START, read a Gmsh file
                    clean_line = line.replace(';', '').replace('\"', '').replace('\'', '')
                    gmsh_file_name = clean_line.split('=')[1]
                    self.read_gmsh_file(gmsh_file_name)
                elif line.startswith(NODES_START):  # If the line starts with NODES_START, read the node coordinates
                    self.read_node_coords(f)
                elif line.startswith(GROUP_START):  # If the line starts with GROUP_START, read a node group
                    if 'name' in line:
                        clean_line = line.replace('>', '').replace('\"', '').replace('\'', '')
                        label = clean_line.split('=')[1]
                        self.read_node_group(f, label)

            # Remove duplicate entries in each group
            for key in self.groups:
                self.groups[key] = list(set(self.groups[key]))

    def read_gmsh_file(self, file_name: str) -> None:
        """
        Read a Gmsh file and add the nodes and groups to the NodeSet object.

        Parameters:
            file_name (str): The name of the Gmsh file to read.

        Returns:
            None
        """
        mesh = meshio.read(file_name, file_format='gmsh')
        obj3d = ['pris', 'pyra', 'hexa', 'wedg', 'tetr']  # A list of 3D mesh types
        self.rank = 2  # Default rank is 2 (for 2D meshes)
        for cell_set in mesh.cell_sets_dict:
            for mesh_type in mesh.cell_sets_dict[cell_set]:
                if mesh_type[:4] in obj3d:
                    self.rank = 3  # If any 3D mesh type is found, set the rank to 3
        for node_id, coords in enumerate(mesh.points):
            self.add_item_by_id(node_id, coords[:self.rank])
            # Add the node with id and coordinates to the NodeSet object
        for cell_set in mesh.cell_sets_dict:
            if cell_set == 'gmsh:bounding_entities':
                pass
            else:
                for mesh_type in mesh.cell_sets_dict[cell_set]:
                    for id_ in mesh.cell_sets_dict[cell_set][mesh_type]:
                        cell_nodes = mesh.cells_dict[mesh_type][id_]
                        for node_id in cell_nodes:
                            self.add_to_group_by_id(cell_set, node_id)  # Add each node to the appropriate group

    def add_to_group_by_id(self, cell_set: str, node_id: int) -> None:
        """
        Add a node to a group.

        Parameters:
            cell_set (str): The name of the group.
            node_id (int): The id of the node to add.

        Returns:
            None
        """
        if cell_set not in self.groups:
            self.groups[cell_set] = [int(node_id)]
        else:
            self.groups[cell_set].append(int(node_id))

    def read_node_coords(self, f: TextIO) -> None:
        """
        Read the node coordinates from the file object f.

        Parameters:
            f (TextIO): The file object to read.

        Returns:
            None
        """
        while True:
            line = f.readline().strip()  # Remove leading/trailing white space
            if line.replace(' ', '').startswith(NODES_END):
                # If the line starts with NODES_END, return from the function
                return

            items = [x for x in re.split(r';|\s', line) if x]  # Split the line into items

            first_item = items[0] if items else ''
            if any(first_item.startswith(comment) for comment in COMMENT_STARTERS):
                # If the line is a comment, skip it
                continue

            if items and items[0].isdigit():
                # If the first item is a digit, add the node coordinates to the NodeSet object
                if self.rank == -1:
                    self.rank = len(items) - 1
                self.add_item_by_id(int(items[0]), [float(crd) for crd in items[1:]])

    def read_node_group(self, f: TextIO, key: str) -> None:
        """
        Read a node group from the file object f with the given key.

        Parameters:
            f (TextIO): The file object to read.
            key (str): The name of the group to add the node to.

        Returns:
            None
        """
        while True:
            line = f.readline()
            if line.replace(' ', '').startswith(GROUP_END):
                # If the line starts with GROUP_END, return from the function
                return
            items = line.split()
            for item in items:
                if item.isdigit():  # If the item can be translated to an integer, add it to the node group
                    self.add_to_group_by_id(key, int(item))


if __name__ == "__main__":
    nset = NodeSet()
    os.chdir('F:\\Github\\pyfem\\examples\\gmsh\\rectangle')
    nset.read_from_file('rectangle.dat')
    print(nset)
    nset = NodeSet()
    os.chdir('F:\\Github\\pyfem\\examples\\mesh')
    nset.read_from_file('PatchTest8_3D.dat')
    print(nset)
