from copy import deepcopy
from typing import List, Union, TextIO, Iterator
import scipy.linalg
from numpy import array, dot, zeros, where
from scipy.sparse.linalg import eigsh
from scipy.sparse.linalg import spsolve

from pyfem.fem.Constrain import Constrain
from pyfem.fem.ElementSet import ElementSet
from pyfem.utils.IntegerIdDict import IntegerIdDict
from pyfem.utils.logger import get_logger
from pyfem.utils.parser import read_node_table, NodeTable

logger = get_logger()


class DofSpace:
    def __init__(self, elements: ElementSet) -> None:
        self.constrain = None
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
            for name in self.constrain.constrained_factors.keys():
                self.constrain.constrained_factors[name] = factor
        else:
            self.constrain.constrained_factors[load_case] = factor

    def read_from_file(self, file_name: str) -> None:
        logger.info("Reading constraints ..........")
        node_table = read_node_table(file_name, "NodeConstraints", self.nodes)
        self.constrain = self.create_constrain(node_table)

    def create_constrain(self, node_tables: Union[List[NodeTable], None] = None) -> Constrain:
        constrain = Constrain(self.get_number_of_dofs())
        if node_tables is None:
            label = "main"
            constrain.constrained_dofs[label] = []
            constrain.constrained_values[label] = []
            constrain.constrained_factors[label] = 1.0
            return constrain

        for node_table in node_tables:
            label = node_table.sub_label

            constrain.constrained_dofs[label] = []
            constrain.constrained_values[label] = []
            constrain.constrained_factors[label] = 1.0

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
                    constrain.add_constraint(dof_id, val, label)

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

                    constrain.add_constraint(dof_id, [val, slave_dof, factor], label)

        constrain.check_constraints(self, node_tables)

        constrain.flush()

        return constrain

    def get_dof_ids_by_type(self, node_ids: Union[int, List[int]], dof_type: str) -> Union[int, List[int]]:
        """
        Get dof_id (or list of dof_ids) for a given dof_type and node (or list of nodes)
        :param node_ids:
        :param dof_type:
        :return:
        """
        return self.dofs[self.id_map.get_items_by_ids(node_ids), self.dof_types.index(dof_type)]

    def get_dof_ids_by_types(self, node_ids: List[int], dof_types: List[str]) -> List[int]:
        """
        Get list of dof_ids for given list of dof_types and list of nodes
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
        Returns the dof_id as a string. For example 'u[14]'
        """
        return self.getTypeName(dof_id) + '[' + str(self.getNodeID(dof_id)) + ']'

    def getNodeID(self, dof_id):

        '''
        Returns the node ID of dof_id
        '''

        return self.nodes.get_id_by_index(int(where(self.dofs == dof_id)[0]))

    def get_type(self, dof_id):

        '''
        Returns the type of dof_id
        '''

        return int(where(self.dofs == dof_id)[1])

    def getTypeName(self, dof_id):

        '''
        Returns the name of the dof_type
        '''

        return self.dof_types[self.get_type(dof_id)]

    def get_items_by_ids(self, node_ids):

        '''Returns all dof_ids for a list of nodes'''

        return self.dofs[self.id_map.get_items_by_ids(node_ids)].flatten()

    def copyConstrain(self, dof_types: list = None):

        '''
        
        '''

        newCons = deepcopy(self.constrain)

        if type(dof_types) is str:
            dof_types = [dof_types]

        for dof_type in dof_types:
            for iDof in self.dofs[:, self.dof_types.index(dof_type)]:
                for label in newCons.constrained_factors.keys():
                    newCons.add_constraint(iDof, 0.0, label)

        newCons.flush()

        return newCons

    def solve(self, A, b, constrainer=None):

        '''Solves the system Ax = b using the internal constraints matrix.
           Returns the total solution vector x.'''

        if constrainer is None:
            constrainer = self.constrain

        if len(A.shape) == 2:

            a = zeros(self.number_of_dofs)

            constrainer.addConstrainedValues(a)

            A_constrained = constrainer.C.transpose() * (A * constrainer.C)

            b_constrained = constrainer.C.transpose() * (b - A * a)

            x_constrained = spsolve(A_constrained, b_constrained)

            x = constrainer.C * x_constrained

            constrainer.addConstrainedValues(x)

        elif len(A.shape) == 1:
            x = b / A

            constrainer.setConstrainedValues(x)

        return x

    def eigensolve(self, A, B, count=5):

        '''Calculates the first count eigenvalues and eigenvectors of a
           system with ( A lambda B ) x '''

        A_constrained = dot(dot(self.constrain.C.transpose(), A), self.constrain.C)
        B_constrained = dot(dot(self.constrain.C.transpose(), B), self.constrain.C)

        eigvals, eigvecs = eigsh(A_constrained, count, B_constrained, sigma=0., which='LM')

        x = zeros(shape=(len(self), count))

        for i, psi in enumerate(eigvecs.transpose()):
            x[:, i] = self.constrain.C * psi

        return eigvals, x

    def norm(self, r, constrainer=None):

        '''
        Calculates the norm of vector r excluding the constrained dofs
        '''

        if constrainer is None:
            constrainer = self.constrain

        return scipy.linalg.norm(constrainer.C.transpose() * r)

    def maskPrescribed(self, a, val=0.0, constrainer=None):

        '''
        Replaced the prescribed dofs by val
        '''

        if constrainer is None:
            constrainer = self.constrain

        a[constrainer.constrained_dofs["None"]] = val

        return a


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
    print(dofs.get_dof_name_by_id(0))
    print(dofs.constrain.constrained_dofs)
    print(dofs.constrain.constrained_values)
    print(dofs.constrain.constrained_factors)
    print(dofs.get_dof_ids_by_types([0], ['u', 'v']))
