from typing import List, Union

import numpy as np
from numpy import array, zeros, where
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

from pyfem.fem.Constraint import Constraint
from pyfem.fem.ElementSet import ElementSet
from pyfem.utils.IntegerIdDict import IntegerIdDict
from pyfem.utils.logger import get_logger
from pyfem.utils.parser import read_node_table, NodeTable

logger = get_logger()


class DofSpace:
    def __init__(self, elements: ElementSet) -> None:
        self.constraint = None
        self.dof_types = elements.get_dof_types()
        self.dofs = array(list(range(len(elements.nodes) * len(self.dof_types)))).reshape(
            (len(elements.nodes), len(self.dof_types)))
        self.nodes = elements.nodes
        self.id_map = IntegerIdDict()
        for index, id_ in enumerate(elements.nodes):
            self.id_map.add_item_by_id(id_, index)
        self.all_constrained_dofs = []
        self.number_of_dofs = self.get_number_of_dofs()

    def __str__(self) -> str:
        return str(self.dofs)

    def get_number_of_dofs(self) -> int:
        return len(self.dofs.flatten())

    def set_constrain_factor(self, factor, load_case='All_'):
        if load_case == 'All_':
            for name in self.constraint.constrained_factors.keys():
                self.constraint.constrained_factors[name] = factor
        else:
            self.constraint.constrained_factors[load_case] = factor

    def read_from_file(self, file_name: str) -> None:
        logger.info("Reading constraints ..........")
        node_table = read_node_table(file_name, "NodeConstraints", self.nodes)
        self.constraint = self.create_constrain(node_table)

    def create_constrain(self, node_tables: Union[List[NodeTable], None] = None) -> Constraint:
        constraint = Constraint(self.get_number_of_dofs())
        if node_tables is None:
            label = "main"
            constraint.constrained_dofs[label] = []
            constraint.constrained_values[label] = []
            constraint.constrained_factors[label] = 1.0
            return constraint

        for node_table in node_tables:
            label = node_table.sub_label
            constraint.constrained_dofs[label] = []
            constraint.constrained_values[label] = []
            constraint.constrained_factors[label] = 1.0

            for item in node_table.data:
                dof_type = item[0]
                node_id = item[1]
                val = item[2]

                if node_id not in self.nodes:
                    raise RuntimeError('Node ID ' + str(node_id) + ' does not exist')

                index = self.id_map.get_items_by_ids(node_id)

                if dof_type not in self.dof_types:
                    raise RuntimeError('DOF type "' + dof_type + '" does not exist')

                if len(item) == 3:
                    dof_id = self.dofs[index, self.dof_types.index(dof_type)]
                    constraint.add_constraint(dof_id, val, label)

                else:
                    slave_node_id = item[4]
                    slave_dof_type = item[3]
                    factor = item[5]

                    if not slave_node_id[0] in self.nodes:
                        raise RuntimeError('Node ID ' + str(slave_node_id) + ' does not exist')

                    slave_index = self.id_map.get_items_by_ids(slave_node_id)

                    if slave_dof_type not in self.dof_types:
                        raise RuntimeError('DOF type "' + slave_dof_type + '" does not exist')

                    slave_dof = self.dofs[slave_index, self.dof_types.index(slave_dof_type)]

                    dof_id = self.dofs[index, self.dof_types.index(dof_type)]

                    constraint.add_constraint(dof_id, [val, slave_dof, factor], label)

        constraint.check_constraints(self, node_tables)

        constraint.flush()

        return constraint

    def get_dof_ids_by_type(self, node_ids: Union[int, List[int]], dof_type: str) -> Union[int, List[int]]:
        """
        Get dof_id (or list of dof_ids) for a given dof_type and node (or list of nodes).
        :param node_ids:
        :param dof_type:
        :return:
        """
        return self.dofs[self.id_map.get_items_by_ids(node_ids), self.dof_types.index(dof_type)]

    def get_dof_ids_by_types(self, node_ids: List[int], dof_types: List[str]) -> List[int]:
        """
        Get list of dof_ids for given list of dof_types and list of nodes.
        :param node_ids:
        :param dof_types:
        :return:
        """
        dofs_ = []
        for node_id in node_ids:
            for dof_type in dof_types:
                dofs_.append(self.dofs[self.id_map.get_items_by_ids(node_id), self.dof_types.index(dof_type)])
        return dofs_

    def get_dof_name_by_id(self, dof_id: int) -> str:
        """
        get the dof name as a string. For example 'u[0]'.
        :param dof_id:
        :return:
        """
        return self.get_dof_type_name(dof_id) + '[' + str(self.get_node_id_by_dof_id(dof_id)) + ']'

    def get_node_id_by_dof_id(self, dof_id: int) -> int:
        return self.nodes.get_id_by_index(int(where(self.dofs == dof_id)[0]))

    def get_dof_type_id(self, dof_id: int) -> int:
        return int(where(self.dofs == dof_id)[1])

    def get_dof_type_name(self, dof_id: int) -> str:
        return self.dof_types[self.get_dof_type_id(dof_id)]

    # def get_items_by_node_ids(self, node_ids: Union[int, List[int]]) -> List[int]:
    #     """
    #     The function is not used.
    #     """
    #     return self.dofs[self.id_map.get_items_by_ids(node_ids)].flatten()

    # def copy_constrain(self, dof_types: list = None) -> Constraint:
    #     """
    #     The function is not used.
    #     """
    #     new_constrain = deepcopy(self.constraint)
    #
    #     if type(dof_types) is str:
    #         dof_types = [dof_types]
    #
    #     for dof_type in dof_types:
    #         for dof_id in self.dofs[:, self.dof_types.index(dof_type)]:
    #             for label in new_constrain.constrained_factors.keys():
    #                 new_constrain.add_constraint(dof_id, 0.0, label)
    #
    #     new_constrain.flush()
    #
    #     return new_constrain

    def solve(self, A: coo_matrix, rhs: np.ndarray, constraint: Constraint = None) -> np.ndarray:
        """
        Solves the system Ax = rhs using the internal constraint matrix.
        Returns the total solution vector x.

        :param A:
        :param rhs:
        :param constraint:
        :return x:
        """

        if constraint is None:
            constraint = self.constraint

        if len(A.shape) == 2:

            a = zeros(self.number_of_dofs)

            constraint.add_constrained_values(a)

            A_constrained = constraint.C.transpose() * (A * constraint.C)

            rhs_constrained = constraint.C.transpose() * (rhs - A * a)

            x_constrained = spsolve(A_constrained, rhs_constrained)

            x = constraint.C * x_constrained

            constraint.add_constrained_values(x)

        elif len(A.shape) == 1:

            x = rhs / A

            constraint.set_constrained_values(x)

        else:
            raise ValueError

        return x

    # def eigen_solve(self, A: coo_matrix, B: coo_matrix, count: int = 5) -> Tuple[np.ndarray]:
    #     """
    #     Calculates the first count eigenvalues and eigenvectors of a system with ( A lambda B ) x
    #     :param A:
    #     :param B:
    #     :param count:
    #     :return:
    #     """
    #
    #     A_constrained = dot(dot(self.constraint.C.transpose(), A), self.constraint.C)
    #     B_constrained = dot(dot(self.constraint.C.transpose(), B), self.constraint.C)
    #
    #     eigen_values, eigen_vectors = eigsh(A_constrained, count, B_constrained, sigma=0., which='LM')
    #
    #     x = zeros(shape=(self.number_of_dofs, count))
    #
    #     for i, psi in enumerate(eigen_vectors.transpose()):
    #         x[:, i] = self.constraint.C * psi
    #
    #     return eigen_values, x

    # def norm(self, r: np.ndarray, constraint: Constraint = None) -> np.ndarray:
    #     """
    #     Calculates the norm of vector r excluding the constrained dofs
    #     :param r:
    #     :param constraint:
    #     :return:
    #     """
    #     if constraint is None:
    #         constraint = self.constraint
    #     return scipy.linalg.norm(constraint.C.transpose() * r)

    # def mask_prescribed(self, a, val: float = 0.0, constraint: Constraint = None):
    #     """
    #     Replaced the prescribed dofs by val
    #     :param a:
    #     :param val:
    #     :param constraint:
    #     :return:
    #     """
    #     if constraint is None:
    #         constraint = self.constraint
    #
    #     a[constraint.constrained_dofs['None']] = val
    #
    #     return a


if __name__ == "__main__":
    from pyfem.utils.parser import file_parser
    from pyfem.fem.NodeSet import NodeSet
    import time
    import os

    t1 = time.time()

    # os.chdir('F:\\Github\\pyfem\\examples\\gmsh\\rectangle')
    # props = file_parser('rectangle.pro')
    # nset = NodeSet()
    # nset.read_from_file('rectangle.dat')
    # elset = ElementSet(nset, props)
    # elset.read_from_file('rectangle.dat')
    # print(elset.elements_count_in_group('All'))

    os.chdir('F:\\Github\\pyfem\\examples\\mesh')
    props = file_parser('PatchTest8_3D.pro')
    nset = NodeSet()
    nset.read_from_file('PatchTest8_3D.dat')
    elset = ElementSet(nset, props)
    elset.read_from_file('PatchTest8_3D.dat')
    dofs = DofSpace(elset)
    dofs.read_from_file('PatchTest8_3D.dat')

    # for e in elset.iter_element_group(['ContElem','ContElem2']):
    #     print(type(e))
    # print(IntegerIdDict.add_item_by_id)

    t2 = time.time()

    total = t2 - t1
    print("Time elapsed = ", total, " [s].\n")
    print(type(dofs.constraint.constrained_dofs['None']))
    # print(dofs.constraint.constrained_dofs)
    # print(dofs.constraint.constrained_values)
    # print(dofs.constraint.constrained_factors)
    # print(dofs.get_dof_ids_by_types([0], ['u', 'v']))
