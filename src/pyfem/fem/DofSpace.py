from copy import deepcopy

import scipy.linalg
from numpy import array, dot, zeros, where
from scipy.sparse.linalg import eigsh
from scipy.sparse.linalg import spsolve

from pyfem.fem.Constrainer import Constrainer
from pyfem.utils.IntegerIdDict import IntegerIdDict
from pyfem.utils.logger import get_logger
from pyfem.utils.parser import read_node_table

logger = get_logger()


class DofSpace:
    def __init__(self, elements):
        self.dof_types = elements.get_dof_types()
        self.dofs = array(list(range(len(elements.nodes) * len(self.dof_types)))).reshape(
            (len(elements.nodes), len(self.dof_types)))
        self.nodes = elements.nodes

        # Create the ID map
        self.IDmap = IntegerIdDict()
        for ind, ID in enumerate(elements.nodes):
            self.IDmap.add_item_by_id(ID, ind)

        self.allConstrainedDofs = []

    def __str__(self):
        return str(self.dofs)

    def __len__(self):
        return len(self.dofs.flatten())

    def setConstrainFactor(self, fac, loadCase="All_"):

        if loadCase == "All_":
            for name in self.cons.constrainedFac.keys():
                self.cons.constrainedFac[name] = fac
        else:
            self.cons.constrainedFac[loadCase] = fac

    def read_from_file(self, fname):

        logger.info("Reading constraints ..........")

        NodeTable = read_node_table(fname, "NodeConstraints", self.nodes)

        self.cons = self.createConstrainer(NodeTable)

    def createConstrainer(self, nodeTables=None):

        cons = Constrainer(len(self))

        if nodeTables == None:
            label = "main"
            cons.constrainedDofs[label] = []
            cons.constrainedVals[label] = []
            cons.constrainedFac[label] = 1.0

            self.cons = cons
            return cons

        for NodeTable in nodeTables:

            label = NodeTable.sub_label

            cons.constrainedDofs[label] = []
            cons.constrainedVals[label] = []
            cons.constrainedFac[label] = 1.0

            for item in NodeTable.data:

                node_id = item[1]
                dof_type = item[0]
                val = item[2]

                if not node_id in self.nodes:
                    raise RuntimeError('Node ID ' + str(node_id) + ' does not exist')

                ind = self.IDmap.get_items_by_ids(node_id)

                if dof_type not in self.dof_types:
                    raise RuntimeError('DOF type "' + dof_type + '" does not exist')

                if len(item) == 3:
                    dofID = self.dofs[ind, self.dof_types.index(dof_type)]

                    cons.addConstraint(dofID, val, label)
                else:
                    slaveNodeID = item[4]
                    slaveDofType = item[3]
                    factor = item[5]

                    if not slaveNodeID[0] in self.nodes:
                        raise RuntimeError('Node ID ' + str(slaveNodeID) + ' does not exist')

                    slaveInd = self.IDmap.get_items_by_ids(slaveNodeID)

                    if slaveDofType not in self.dof_types:
                        raise RuntimeError('DOF type "' + slaveDofType + '" does not exist')

                    slaveDof = self.dofs[slaveInd, self.dof_types.index(slaveDofType)]

                    dofID = self.dofs[ind, self.dof_types.index(dof_type)]

                    cons.addConstraint(dofID, [val, slaveDof, factor], label)

                    # Check for all tyings whether master of slave is not slave itself
        cons.checkConstraints(self, nodeTables)

        cons.flush()

        return cons

    def getForType(self, node_ids, dof_type):

        '''
        Returns all dofIDs for given dof_type for a list of nodes
        '''

        return self.dofs[self.IDmap.get_items_by_ids(node_ids), self.dof_types.index(dof_type)]

    def getForTypes(self, node_ids, dof_types):

        '''
        Returns all dofIDs for given list of dof_type for a list of nodes
        '''

        dofs = []

        for node in node_ids:
            for dof_type in dof_types:
                dofs.append(self.dofs[self.IDmap.get_items_by_ids(node), self.dof_types.index(dof_type)])

        return dofs

    def getDofName(self, dofID):

        '''
        Returns the dofID as a string. For example 'u[14]'
        '''

        return self.getTypeName(dofID) + '[' + str(self.getNodeID(dofID)) + ']'

    def getNodeID(self, dofID):

        '''
        Returns the node ID of dofID
        '''

        return self.nodes.get_id_by_index(int(where(self.dofs == dofID)[0]))

    def get_type(self, dofID):

        '''
        Returns the type of dofID
        '''

        return int(where(self.dofs == dofID)[1])

    def getTypeName(self, dofID):

        '''
        Returns the name of the dof_type
        '''

        return self.dof_types[self.get_type(dofID)]

    def get_items_by_ids(self, node_ids):

        '''Returns all dofIDs for a list of nodes'''

        return self.dofs[self.IDmap.get_items_by_ids(node_ids)].flatten()

    def copyConstrainer(self, dof_types: list = None):

        '''
        
        '''

        newCons = deepcopy(self.cons)

        if type(dof_types) is str:
            dof_types = [dof_types]

        for dof_type in dof_types:
            for iDof in self.dofs[:, self.dof_types.index(dof_type)]:
                for label in newCons.constrainedFac.keys():
                    newCons.addConstraint(iDof, 0.0, label)

        newCons.flush()

        return newCons

    def solve(self, A, b, constrainer=None):

        '''Solves the system Ax = b using the internal constraints matrix.
           Returns the total solution vector x.'''

        if constrainer is None:
            constrainer = self.cons

        if len(A.shape) == 2:

            a = zeros(len(self))

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

        A_constrained = dot(dot(self.cons.C.transpose(), A), self.cons.C)
        B_constrained = dot(dot(self.cons.C.transpose(), B), self.cons.C)

        eigvals, eigvecs = eigsh(A_constrained, count, B_constrained, sigma=0., which='LM')

        x = zeros(shape=(len(self), count))

        for i, psi in enumerate(eigvecs.transpose()):
            x[:, i] = self.cons.C * psi

        return eigvals, x

    def norm(self, r, constrainer=None):

        '''
        Calculates the norm of vector r excluding the constrained dofs
        '''

        if constrainer is None:
            constrainer = self.cons

        return scipy.linalg.norm(constrainer.C.transpose() * r)

    def maskPrescribed(self, a, val=0.0, constrainer=None):

        '''
        Replaced the prescribed dofs by val
        '''

        if constrainer is None:
            constrainer = self.cons

        a[constrainer.constrainedDofs["None"]] = val

        return a
