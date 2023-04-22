import os
import re
from typing import Dict, List, Union, TextIO

import meshio
import numpy as np

from pyfem.utils.item_list import ItemList
from pyfem.utils.logger import getLogger
from pyfem.utils.parser import get_type

logger = getLogger()

NODES_START = "<Nodes>"
NODES_END = "</Nodes>"
GROUP_START = "<NodeGroups>"
GROUP_END = "</NodeGroups>"
GMSH_START = "gmsh"
COMMENT_STARTERS = ["//", "#"]


class NodeSet(ItemList):

    def __init__(self):
        super().__init__()
        self.rank: int = -1
        self.groups: Dict[str, List[int]] = {}

    def __repr__(self) -> str:
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
        return np.array(self.get(node_ids))

    def read_from_file(self, file_name: str) -> None:
        logger.info("Reading nodes ................")

        with open(file_name, 'r') as f:
            for line in f:
                line = line.strip().replace(' ', '').replace('\n', '').replace('\t', '').replace('\r', '')
                if line.startswith(NODES_START):
                    self.read_node_coords(f)
                elif line.startswith(GMSH_START):
                    clean_line = line.replace(';', '').replace('\"', '').replace('\'', '')
                    gmsh_file_name = clean_line.split('=')[1]
                    self.read_gmsh_file(gmsh_file_name)
                elif line.startswith(GROUP_START):
                    if 'name' in line:
                        clean_line = line.replace('>', '').replace('\"', '').replace('\'', '')
                        label = clean_line.split('=')[1]
                        self.read_node_group(f, label)

            for key in self.groups:
                self.groups[key] = list(set(self.groups[key]))

    def read_gmsh_file(self, file_name: str) -> None:
        mesh = meshio.read(file_name, file_format='gmsh')
        obj3d = ['pris', 'pyra', 'hexa', 'wedg', 'tetr']
        self.rank = 2
        for cell_set in mesh.cell_sets_dict:
            for mesh_type in mesh.cell_sets_dict[cell_set]:
                if mesh_type[:4] in obj3d:
                    self.rank = 3
        for node_id, coords in enumerate(mesh.points):
            self.add(node_id, coords[:self.rank])
        for cell_set in mesh.cell_sets_dict:
            if cell_set == 'gmsh:bounding_entities':
                pass
            else:
                for mesh_type in mesh.cell_sets_dict[cell_set]:
                    for id_ in mesh.cell_sets_dict[cell_set][mesh_type]:
                        cell_nodes = mesh.cells_dict[mesh_type][id_]
                        for node_id in cell_nodes:
                            self.add_to_group(cell_set, node_id)

    def add_to_group(self, cell_set: str, node_id: int) -> None:
        if cell_set not in self.groups:
            self.groups[cell_set] = [int(node_id)]
        else:
            self.groups[cell_set].append(int(node_id))

    def read_node_coords(self, f: TextIO) -> None:
        while True:
            line = f.readline().strip()
            if line.replace(' ', '').startswith(NODES_END):
                return

            # items = [x.strip() for x in re.split(r';|\s', line) if x]
            #
            # print(items)
            # for item in items:
            #     if any(item.startswith(comment) for comment in COMMENT_STARTERS):
            #         break
            #     parts = item.split()
            #     print(parts)
            #
            #     if len(parts) > 1 and isinstance(eval(parts[0]), int):
            #         if self.rank == -1:
            #             self.rank = len(parts) - 1
            #         self.add(eval(parts[0]), [eval(crd) for crd in parts[1:]])

            line = re.sub('\s{2,}', ' ', line)
            a = line.split(';')
            for a in a[:-1]:
                b = a.strip().split(' ')
                for comment in COMMENT_STARTERS:
                    if line.startswith(comment):
                        break
                if len(b) > 1 and type(eval(b[0])) == int:
                    if self.rank == -1:
                        self.rank = len(b) - 1
                    self.add(eval(b[0]), [eval(crd) for crd in b[1:]])

    def read_node_group(self, f: TextIO, key: str) -> None:
        while True:
            line = f.readline()

            if line.replace(" ", "").startswith('</NodeGro'):
                return

            a = line.split()

            for b in a:
                if get_type(b) == int:
                    self.add_to_group(key, b)


if __name__ == "__main__":
    nset = NodeSet()
    # os.chdir('F:\\Github\\pyfem\\examples\\gmsh\\rectangle')
    # nset.read_from_file('rectangle.dat')
    os.chdir('F:\\Github\\pyfem\\examples\\mesh')
    nset.read_from_file('PatchTest8_3D.dat')
    print(nset.__dict__)
    print(nset.items())
    # print(type(nset.get_node_coords(1)))
